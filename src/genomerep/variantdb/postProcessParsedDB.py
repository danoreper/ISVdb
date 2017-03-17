import logging.config
import yaml
print("wakkea!")
import ipdb

logging.config.dictConfig(yaml.load(open("../config/logging.config",'r')))
logger=logging.getLogger("mylog")
logger.debug("logging!")

def createIndex(con, tablename, cols):
    query = "CREATE INDEX " + "".join(cols) + " on " + tablename + "(" + ", ".join(cols) +")";
    con.query(query)

def buildStrainsTable(con, engine):
    createIndex(con, "founder_genotype", ["strain_name"])
    logger.debug("created founder_genotype_index")
    return con.query('''create TABLE strain 
                        (strain_id smallint unsigned NOT NULL AUTO_INCREMENT PRIMARY KEY, 
                        UNIQUE INDEX ix_strain(strain_name) USING HASH)  
                        ENGINE = ''' + engine + '''  
                      select distinct strain_name from founder_genotype;''')

def refactorFounderTable(con, engine):
    con.query("drop table if exists acop;")
    con.query('''create  table acop ENGINE = ''' + engine + '''
             select variant_id, strain_id, allele_index_1, allele_index_2, prob, is_max 
             FROM founder_genotype 
             inner join strain on founder_genotype.strain_name=strain.strain_name ;''')
    con.query("drop   table founder_genotype;")
    con.query("rename table acop TO founder_genotype;")
    
    createIndex(con, "founder_genotype", ["variant_id", "allele_index_1","allele_index_2","strain_id","prob"])
    createIndex(con, "founder_genotype", ["strain_id", "variant_id","allele_index_1","allele_index_2","prob"])
    createIndex(con, "founder_genotype", ["strain_id", "allele_index_1","allele_index_2","prob"])
    createIndex(con, "founder_genotype", ["allele_index_1", "allele_index_2","strain_id","variant_id","prob"])
    createIndex(con, "founder_genotype", ["allele_index_2","strain_id","variant_id","prob"])
        
    con.query("ALTER TABLE founder_genotype ADD CONSTRAINT fk_strain_id FOREIGN KEY(strain_id) REFERENCES strain(strain_id);")


def postProcessFounderMeta(con, engine):
    con.query("CREATE INDEX ix_strain_name on founder_meta(strain_name) USING HASH;")
    con.query("drop table if exists acop;")
    con.query('''create  table acop ENGINE = ''' + engine + '''
             select strain.strain_id, founder_meta.*   
             FROM founder_meta 
             inner join strain on founder_meta.strain_name=strain.strain_name ;''')
    con.query("alter table acop drop column strain_name;")
    con.query("drop table founder_meta;")
    con.query("rename table acop TO founder_meta;")
    con.query("CREATE INDEX ix_strain_id_variant_id on founder_meta(strain_id, variant_id);")
    con.query("CREATE INDEX ix_variant_id_strain_id on founder_meta(variant_id,strain_id);")


def postProcessInfo(con, allowedTranscripts,engine):
    con.query("CREATE INDEX ix_gene_name on info(gene_name) USING HASH;")
    con.query("CREATE INDEX ix_transcript_name on info(transcript_name) USING HASH;")
    con.query("CREATE INDEX ix_consequence_name on info(consequence_name) USING HASH;")
#get rid of transcripts and consequences that are not in the allowed transcripts list
    
    if len(allowedTranscripts)>1:
        con.query("drop table if exists allowed_transcripts")
        con.query("CREATE TABLE allowed_transcripts(transcript_name varchar(100)) engine = " + engine + ";")
        valueString = ",".join(["('" + str(x) + "')" for x in allowedTranscripts])
        con.query("INSERT into allowed_transcripts(transcript_name) values " + valueString + ";")
        con.query("create index ix on allowed_transcripts(transcript_name);")
        con.query("drop table if exists info2;")
        con.query("create table info2 select info.* from info, allowed_transcripts where info.transcript_name=allowed_transcripts.transcript_name;") #         x=raw_input("enter")
        con.query("drop table info;")
        con.query("drop table allowed_transcripts;")
        con.query("rename table info2 to info")
    
    con.query("create TABLE gene        (gene_id        int NOT NULL AUTO_INCREMENT PRIMARY KEY, UNIQUE INDEX ix_gene_name(gene_name)               USING HASH) ENGINE = " + engine + " select distinct gene_name        from info;")
    con.query("create TABLE transcript  (transcript_id  int NOT NULL AUTO_INCREMENT PRIMARY KEY, UNIQUE INDEX ix_transcript_name(transcript_name)   USING HASH) ENGINE = " + engine + " select distinct transcript_name  from info;")
    con.query("create TABLE consequence (consequence_id int NOT NULL AUTO_INCREMENT PRIMARY KEY, UNIQUE INDEX ix_consequence_name(consequence_name) USING HASH) ENGINE = " + engine + " select distinct consequence_name from info;")
    con.query('''create table acop ENGINE = ''' + engine + ''' select info.variant_id as variant_id, allele_index, gene_id, transcript_id, consequence_id 
    FROM info 
    inner join gene on info.gene_name=gene.gene_name 
    inner join transcript on info.transcript_name = transcript.transcript_name
    inner join consequence on info.consequence_name = consequence.consequence_name
    inner join allele on info.allele=allele.allele AND info.variant_id = allele.variant_id;''')
    con.query("drop  table info;")
    con.query("rename table acop TO consequence_info;")
    con.query("CREATE INDEX ix_allele_index on consequence_info(allele_index);")
    con.query("CREATE INDEX ix_transcript_id on consequence_info(transcript_id) USING HASH;")
    con.query("CREATE INDEX ix_consequence_id on consequence_info(consequence_id) USING HASH;")
    con.query("CREATE INDEX ix_variant_id_consequence_id on consequence_info(variant_id, consequence_id) USING HASH;")
    con.query("CREATE INDEX ix_variant_id_allele_index_consequence_id on consequence_info(variant_id, allele_index, consequence_id) USING HASH;")
    con.query("CREATE INDEX ix_toolong on consequence_info(variant_id, transcript_id, allele_index, consequence_id) USING HASH;")
    con.query("CREATE INDEX ix_gene_id_consequence_id on consequence_info(gene_id, consequence_id) USING HASH;")
    con.query("CREATE INDEX ix_gene_id_allele_index on consequence_info(gene_id, allele_index) USING HASH;")
    con.query("ALTER TABLE consequence_info ADD CONSTRAINT fk_allele_index   FOREIGN KEY(allele_index)   REFERENCES allele(allele_index);")
    con.query("ALTER TABLE consequence_info ADD CONSTRAINT fk_gene_id        FOREIGN KEY(gene_id)        REFERENCES gene(gene_id);")
    con.query("ALTER TABLE consequence_info ADD CONSTRAINT fk_transcript_id  FOREIGN KEY(transcript_id)  REFERENCES transcript(transcript_id);")
    con.query("ALTER TABLE consequence_info ADD CONSTRAINT fk_consequence_id FOREIGN KEY(consequence_id) REFERENCES consequence(consequence_id);")


def add_B6_referenceData(con):
    con.query("INSERT INTO consequence(consequence_name) VALUES('reference');")
    con.commit()
    con.query("INSERT INTO strain(strain_name) VALUES ('C57BL6J');")
    con.commit()
    #allele index 1 and 2 for b6 is 0. The probability is 1. It is always the max
    con.query('''INSERT INTO founder_genotype(strain_id, variant_id, allele_index_1, allele_index_2, prob, is_max) 
            select DISTINCT strain.strain_id, variant.variant_id, 0,0, 1,1 from variant, strain 
            where strain.strain_name='C57BL6J'; ''')
    con.commit()
    #add allele index and its consequence into the var info table
    con.query('''INSERT INTO consequence_info(variant_id, allele_index, gene_id, transcript_id, consequence_id) 
            select DISTINCT variant_id, 0, gene_id, transcript_id, consequence.consequence_id 
            from consequence_info, consequence where consequence.consequence_name='reference'; ''')
    con.commit()

def postprocess(con, allowedTranscripts, engine):
    logger.debug("postproc")
    logger.info("indexing variants")
    createIndex(con, "variant", ["variant_id","chrom","pos"])
    createIndex(con, "variant", ["variant_id","pos","chrom"])
    createIndex(con, "variant", ["pos","chrom","variant_id"])
    createIndex(con, "variant", ["pos","variant_id","chrom"])
    createIndex(con, "variant", ["chrom","variant_id","pos"])
    createIndex(con, "variant", ["chrom","pos","variant_id"])

    logger.info("building strains table")
    buildStrainsTable(con, engine)
    
    logger.info("refactoring founder table")
    refactorFounderTable(con, engine)
    
    logger.info("postproc meta table")
    postProcessFounderMeta(con, engine)
    
    logger.info("postproc info table")
    postProcessInfo(con, allowedTranscripts, engine)
    #add B6 information to database, even though it isnt explicitly in the vcf
    logger.info("adding b6 info")
    add_B6_referenceData(con)
    logger.info("done all post proc")
    con.close()

