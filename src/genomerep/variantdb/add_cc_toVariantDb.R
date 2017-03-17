# TODO: Add comment
# 
# Author: doreper
###############################################################################
library(RMySQL)
library(IRanges)
library(data.table)
library(reshape2)
library("Biostrings")
source("./stringutils.R")
source("./genomerep/cc_founderProbs.R")
source("./loadParams.R")


rootTempFileDir = prop$tmpdir
engine = prop$variantdb$engine
dir.create(rootTempFileDir)

add_unobserved_cclines_to_strains <- function(con, probsFiles) 
{
    
    
    strainMapping = data.table(dbGetQuery(con, "SELECT * from strain;"), key="strain_name")
    
    oldStrains = unlist(lapply(strsplit(basename(probsFiles),split=".csv"), "[",1))
    newStrains = setdiff(oldStrains, strainMapping$strain_name)
    maxId = max(dbGetQuery(con, "SELECT strain_id from strain;"))
    values = paste(paste0("(","'",newStrains,"'", ", ", (1+maxId):(length(newStrains)+maxId),")"),collapse=", ") 
    dbSendQuery(con, paste0("INSERT IGNORE INTO strain (strain_name, strain_id) VALUES ", values,";"))
}

buildInbredLineFounders <- function(con, engine, rilHaplotypeProbsDir, founders, karyotype) 
{
    print("building inbred line founders")
    dbSendQuery(con, "DROP TABLE IF EXISTS diplotype;")
    
    dbSendQuery(con, paste0("CREATE TABLE if not exists diplotype (
							strain_id smallint unsigned NOT NULL,
							variant_id int unsigned NOT NULL, 
							founder_id_1 smallint unsigned NOT NULL,
							founder_id_2 smallint unsigned NOT NULL,
							prob float unsigned not NULL
##,
							
							## FOREIGN KEY(variant_id) REFERENCES variant(variant_id),
							## FOREIGN KEY(strain_id)	REFERENCES strain(strain_id),
							## FOREIGN KEY(founder_id_1) REFERENCES strain(strain_id),
							## FOREIGN KEY(founder_id_1) REFERENCES strain(strain_id)
      ) ENGINE = ",engine,";"))      
	
    variantFrame = dbGetQuery(con, "SELECT * from variant;")
    probsFiles = stringutils$getFilesWithExtension(rilHaplotypeProbsDir, ".csv")
    add_unobserved_cclines_to_strains(con = con, probsFiles)
    print("added unobserved strains")
    strainMapping = data.table(dbGetQuery(con, "SELECT * from strain;"), key="strain_name")
   
    jobBounds = getJobBounds(variantFrame = variantFrame)
    
    cc_limit = prop$variantdb$cc_limit
    if(is.na(cc_limit))
    {
        cc_limit = length(probsFiles)
    }
    
    for(probFile in probsFiles[1:cc_limit])
    {
        strainNameAndDescent = buildFounderProbs$parseProbFilename(probFile)
                
        strainName = paste(strainNameAndDescent$strain, strainNameAndDescent$descent, sep="_")
        probData = buildFounderProbs$build(probFile, founders,karyotype)
        atempGt   = normalizePath(tempfile(paste0("diplos","_",strainName), rootTempFileDir, ".csv"))
      
        for(jb in 1:(length(jobBounds)-1))
        {
            jobStart = jobBounds[jb]
            jobEnd   = jobBounds[jb+1]-1
            print(paste0(jobStart, ":", jobEnd))
            subframe = variantFrame[jobStart:jobEnd,]
            gt = buildFounderProbs$hapProbPerPosition(subframe$variant_id, subframe$chrom, subframe$pos,probData)
            db_version_of_gt = convert_GT_to_DB_format(gt = gt, strainName, strainMapping = strainMapping)
            write.table(as.matrix(db_version_of_gt),   atempGt,   append=T, row.names=F, col.names=F, sep=",")
        }
                                        #once a full file for a given strain is assembled, load it into the database
        fileInsertQuery = paste0("LOAD DATA INFILE '", normalizePath(atempGt), "' INTO TABLE diplotype ",
                                 "FIELDS TERMINATED BY ',' ",
                                 "LINES TERMINATED BY '\n' ",
                                 ";")
        print(fileInsertQuery)
        x=dbSendQuery(con,fileInsertQuery)
        ##then delete the file... we dont want to double up our space by having this redundant version floating around
      unlink(atempGt)
  }
  
                                        #insert all the founder lines into the diplotype_founder table
  dbSendQuery(con, "insert into diplotype 
					select distinct strain_id, variant_id, strain_id as founder_id_1, strain_id as founder_id_2, 1 as prob from founder_genotype;")
  
  
##  createIndex(con, "diplotype", c("founder_id_2"))
##  createIndex(con, "diplotype", c("prob"))
  
  createIndex(con, "diplotype", c("founder_id_1","founder_id_2"))
  createIndex(con, "diplotype", c("strain_id",  "variant_id", "founder_id_1", "founder_id_2", "prob"))
  createIndex(con, "diplotype", c("variant_id", "strain_id",  "founder_id_1", "founder_id_2", "prob"))
  
  createIndex(con, "diplotype", c("variant_id", "founder_id_1", "founder_id_2","prob"))
  createIndex(con, "diplotype", c("strain_id",  "founder_id_1", "founder_id_2","prob"))
  
  
  addForeignKey(con, "diplotype", "strain_id",  "strain", "strain_id")
  addForeignKey(con, "diplotype", "variant_id", "variant", "variant_id")
  addForeignKey(con, "diplotype", "founder_id_1", "founder", "strain_id")
  addForeignKey(con, "diplotype", "founder_id_2", "founder", "strain_id")
}

getJobBounds <- function(variantFrame) 
{
	numVariants = nrow(variantFrame)
	by.increment = 1000000
#divide up the jobs per vcf folder into equal sizes. The end of one job is the (start - 1) of the next job. 
#Cap it off by the last job as necessary given rounding
        jobBounds = seq(1, numVariants, by=by.increment)
	if(length(jobBounds)<(numVariants+1))
	{
		jobBounds = c(jobBounds, numVariants+1)
	}
	return(jobBounds)
}

convert_GT_to_DB_format <- function(gt, strainName, strainMapping) 
{
    melted = data.table::melt(gt,id=c("variant_id"), variable.name = "diplotype", value.name="prob")
    melted = melted[melted$prob>=prop$variantdb$truncate_cc_prob,]
    rm(gt)
    
    chars = as.character(melted$diplotype)
    id1 = strainMapping[sub("\\..*","",chars ),strain_id]
    id2 = strainMapping[sub(".*\\.","",chars ),strain_id]
    melted$founder_id_1 = pmin(id1, id2)
    melted$founder_id_2 = pmax(id1, id2)
    
    melted$strain_id    = strainMapping[J(strain_name=strainName),strain_id]
    melted = data.frame(melted)
    melted = melted[,c("strain_id", "variant_id", "founder_id_1", "founder_id_2", "prob")]
    
    return(melted)
}

#for every pair of founder alleles, build up what the genotype per strain will be. Ensure that redundant genotypes are
#summed over. Creates a table entitled inbred_alleles_no_max ( the next step of processing will add a field indicating 
#whether a given genotype is the most likely one.
buildInbredAlleles_nomax <- function(con)
{

    ##build a table of the alleles for inbred line, assuming only one of the founder1 alleles is
    ##used, and only one of the founder 2 alleles is used. We will call this 4 times to get all 4 possible duplotypes given the founders
    buildInbredAllelesQuery = function(f1_allele_index, f2_allele_index, name)
    {
        q = paste0("select * from (
						select 
						il.strain_id,
						il.variant_id,  
						least(   f1.allele_index_",f1_allele_index, ", f2.allele_index_", f2_allele_index, ") as allele_index_1,
						greatest(f1.allele_index_",f1_allele_index, ", f2.allele_index_", f2_allele_index, ") as allele_index_2,
						il.prob*f1.prob*f2.prob/4 as prob 
						from 
						diplotype il 
						inner join founder_genotype f1 on il.variant_id=f1.variant_id and il.founder_id_1=f1.strain_id
						inner join founder_genotype f2 on il.variant_id=f2.variant_id and il.founder_id_2=f2.strain_id
        ) as ", name)
        return(q)
    }

    buildInbredAllelesQuery = function(f1_allele_index, f2_allele_index, name=NULL)
    {
        q = ""
        if(!is.null(name))
        {
            q = paste0(q, " select * from (")
        }
        
	q = paste0(q, "select 
		      il.strain_id,
                      il.variant_id,  
		      least(   f1.allele_index_",f1_allele_index, ", f2.allele_index_", f2_allele_index, ") as allele_index_1,
                      greatest(f1.allele_index_",f1_allele_index, ", f2.allele_index_", f2_allele_index, ") as allele_index_2,
		      il.prob*f1.prob*f2.prob/4 as prob 
		      from 
		      diplotype il 
	              inner join founder_genotype f1 on il.variant_id=f1.variant_id and il.founder_id_1=f1.strain_id
		      inner join founder_genotype f2 on il.variant_id=f2.variant_id and il.founder_id_2=f2.strain_id")
        
        if(!is.null(name))
        {
            q = paste0(q, ") as ", name)
        }
        return(q)
    }
    
    print("about to run query")
    dbSendQuery(con, "Drop table if exists inbred_alleles_no_max;")
    print("dropped")
    ##take the union of all possible paris of founder alleles- 4 different combinations

    singleQuery = F

    if(singleQuery)
    {
        
        query = paste0("create table inbred_alleles_no_max",
                       " ENGINE = ",engine, " ", 
                       
                       buildInbredAllelesQuery(1,1, "t1"),
                       " UNION ALL ",
                       buildInbredAllelesQuery(1,2, "t2"),
                       " UNION ALL ",
                       buildInbredAllelesQuery(2,1, "t3"),
                       " UNION ALL ",
                       buildInbredAllelesQuery(2,2, "t4"),";")
        print(query)
        dbSendQuery(con, query)
        print("made query")
    } else {
    
        ## ##take the union of all possible paris of founder alleles- 4 different combinations
        ## #Doing this all in a single union all call for some reasons hangs the mariadb. Thus breaking it up into 4 tables, taking the union one table at a time.
        
        dbSendQuery(con, statement = paste("create table inbred_alleles_no_max", "ENGINE = ", engine, buildInbredAllelesQuery(1,1), ";"))

        dbSendQuery(con, statement = paste("insert into inbred_alleles_no_max",  buildInbredAllelesQuery(1,2), ";"))
        dbSendQuery(con, statement = paste("insert into inbred_alleles_no_max",  buildInbredAllelesQuery(2,1), ";"))
        dbSendQuery(con, statement = paste("insert into inbred_alleles_no_max",  buildInbredAllelesQuery(2,2), ";"))
    }
        
	#sum redundant diplotypes before returning the table
	#we unfortunately have to do this in a separate steps so as to first create an index on the derived allele_index columns
	#before applying a potentially slow sum over groupe
    createIndex(con, "inbred_alleles_no_max", c("strain_id", "variant_id", "allele_index_1", "allele_index_2","prob"))
	
    dbSendQuery(con, "Drop table if exists temp_table;")
    query = paste0("create table temp_table(strain_id smallint, variant_id int, allele_index_1 tinyint, allele_index_2 tinyint, prob float)", 
                       " ENGINE = ",engine, " ",
                       "select strain_id, variant_id, allele_index_1, allele_index_2, sum(prob) as prob ", 
                       "from inbred_alleles_no_max ",
                       "group by strain_id, variant_id, allele_index_1, allele_index_2;")
    
    dbSendQuery(con, query)
    dbSendQuery(con, "drop table if exists inbred_alleles_no_max;")
    dbSendQuery(con, "rename table temp_table TO inbred_alleles_no_max;")
    createIndex(con, "inbred_alleles_no_max", c("strain_id", "variant_id", "allele_index_1", "allele_index_2", "prob"))
}


#Creates a table inbred_alleles_only_max consisting of only the max likelihood genotype
#records per strain-variant in inbred_alleles_no_max
buildInbredAlleles_onlyMax <- function(con)
{
    dbSendQuery(con, "Drop table if exists inbred_alleles_only_max;")
    query = paste0("create table inbred_alleles_only_max(
			is_max tinyint, 
			strain_id smallint, 
			variant_id int,
			allele_index_1 tinyint, 
			allele_index_2 tinyint, 
			prob float)", 
                   
                   " ENGINE = ",engine, " ",
                   "select 
			1 as is_max, 
			t1.strain_id, 
			t1.variant_id, 
			t1.allele_index_1, 
			t1.allele_index_2, 
			t1.prob 
			
			from inbred_alleles_no_max t1 
			
			left join inbred_alleles_no_max t2 on 
			t1.strain_id = t2.strain_id and 
			t1.variant_id = t2.variant_id and 
			(t1.prob < t2.prob OR
				(t1.prob=t2.prob and concat_ws(',', t1.allele_index_1, t1.allele_index_2)<
								 	 concat_ws(',', t2.allele_index_1, t2.allele_index_2))
			) where t2.strain_id IS NULL;")
	dbSendQuery(con, query)
}

mergeMaxAndNoMax <- function(con) 
{
	createIndex(con, "inbred_alleles_only_max", 
			c("strain_id", "variant_id", "allele_index_1", "allele_index_2"))

	dbSendQuery(con, "Drop table if exists genotype;")
	query = paste0("create table genotype(strain_id smallint unsigned not null,
                							 variant_id int unsigned not null, 
											 allele_index_1 tinyint unsigned,
											 allele_index_2 tinyint unsigned,
											 prob float unsigned not null,
											 is_max tinyint unsigned not null)",
			" ENGINE = ",engine, " ",
			
			"select 
			t1.strain_id, 			
			t1.variant_id, 
			t1.allele_index_1,
			t1.allele_index_2, 
			t1.prob, 
			ifnull(t2.is_max,0) as is_max 
			
			from inbred_alleles_no_max t1 
			
			left join inbred_alleles_only_max t2 on 
			t1.strain_id     = t2.strain_id and
			t1.variant_id    = t2.variant_id and 
			t1.allele_index_1= t2.allele_index_1 and
			t1.allele_index_2= t2.allele_index_2;")

	dbSendQuery(con, query)
	dbSendQuery(con, "drop table if exists inbred_alleles_only_max;")
	dbSendQuery(con, "drop table if exists inbred_alleles_no_max;")
	
	createIndex(con, "genotype", c("strain_id",  "variant_id",     "allele_index_1", "allele_index_2","is_max","prob"))
	createIndex(con, "genotype", c("variant_id", "strain_id",      "allele_index_1", "allele_index_2","is_max","prob"))
	createIndex(con, "genotype", c("variant_id", "allele_index_1", "allele_index_2", "is_max",  "prob", "strain_id"))
	createIndex(con,  "genotype", c("strain_id", "allele_index_1", "allele_index_2", "is_max","prob","variant_id"))

	## createIndex(con, "genotype", c("strain_id",  "variant_id",  "allele_index_1", "allele_index_2","prob"))
	## createIndex(con, "genotype", c("variant_id", "strain_id",   "allele_index_1", "allele_index_2","prob"))
	## createIndex(con, "genotype", c("variant_id", "allele_index_1", "allele_index_2","prob","strain_id"))
	## createIndex(con,  "genotype", c("strain_id", "allele_index_1", "allele_index_2","prob","variant_id"))
	
	addForeignKey(con, "genotype", "strain_id",  "strain", "strain_id")
	addForeignKey(con, "genotype", "variant_id", "variant", "variant_id")
	addForeignKey(con, "genotype", "allele_index_1", "allele", "allele_index")
	addForeignKey(con, "genotype", "allele_index_2", "allele", "allele_index")
	addForeignKey(con, "genotype", "allele_index_1", "varinfo", "allele_index")
	addForeignKey(con, "genotype", "allele_index_2", "varinfo", "allele_index")
	
}

buildInbredLineFixed <- function(con, engine) 
{
	dbSendQuery(con, "DROP table if exists inbred_line_fixed;")
	query = paste0("create table inbred_line_fixed(strain_id smallint not null, variant_id int not null, fixed tinyint not null)", 
			" ENGINE = ",engine," ", 
			"select strain_id, variant_id, (sum(ilf.founder_id_1!=ilf.founder_id_2))=0 as fixed ",
			"from diplotype ilf ",
			"group by ilf.strain_id, ilf.variant_id;")
	dbSendQuery(con, query)
	
	createIndex(con, "inbred_line_fixed",c("strain_id","variant_id","fixed"))
	createIndex(con, "inbred_line_fixed",c("variant_id","strain_id","fixed"))
	createIndex(con, "inbred_line_fixed",c("strain_id", "fixed"))
	createIndex(con, "inbred_line_fixed",c("variant_id","fixed"))
	createIndex(con, "inbred_line_fixed",c("fixed"))
}

buildHaplotype <- function(con, hap_prefix, tablename)
{
                                        #	
    
    dbSendQuery(con, paste0("drop table if exists ", tablename,"_sampling ;"))
    makeHalfQuery = function(hap_prefix, tablename, halfIndex)
    {
        query = paste0("Select * FROM (Select strain_id, variant_id, ",hap_prefix, "_",halfIndex," as ", hap_prefix, ",", " prob/2 as prob ", "from ", tablename, ") as t", halfIndex) 
        return(query)
    }
    
    dbSendQuery(con, "drop table if exists temp;")
    query = paste0("create table temp", "(strain_id smallint unsigned not null, variant int unsigned not null,", hap_prefix, " tinyint unsigned, prob float unsigned not null) ",
                   " ENGINE = ",engine, " ",
                   makeHalfQuery(hap_prefix, tablename,1),
                   " UNION ALL ",
                   makeHalfQuery(hap_prefix, tablename,2),
                   ";")
    
    dbSendQuery(con, query)
    createIndex(con, "temp", c("strain_id","variant_id",hap_prefix))
    query = paste0("create table ", tablename,"_sampling (prob float unsigned not null)",
                   " ENGINE = ",engine, " ",
                   "select strain_id, variant_id, ", hap_prefix,", sum(prob) as prob from temp group by strain_id, variant_id,", hap_prefix,";")
    print(query)
    dbSendQuery(con, query)
    
    createIndex(con, paste0(tablename,"_sampling "), c("strain_id","variant_id", hap_prefix))
    createIndex(con, paste0(tablename,"_sampling "), c("variant_id","strain_id", hap_prefix))
    createIndex(con, paste0(tablename,"_sampling "), c("variant_id", hap_prefix))
    createIndex(con, paste0(tablename,"_sampling "), c("strain_id", hap_prefix))
    
                                        #TODO create other indexes as necessary
    
                                        #TODO add a max probability haplotype field
    dbSendQuery(con, "drop table if exists temp")
    return(query)
}


#dbSendQuery(con, "CREATE INDEX ix_variant_id on inbred_line(variant_id);");
createIndex <- function(con, tableName, colsToIndex)
{
	astring = paste0("CREATE INDEX ix", paste(colsToIndex, collapse=""), " ON ", tableName, "(", paste(colsToIndex, collapse=","),")") 
	print(astring)
	dbSendQuery(con, astring)
	dbCommit(con)
	return(astring)
}

addForeignKey <- function(con, tablename, col, tablenameref, colref)
{
#	
    query = paste0("ALTER TABLE ", tablename, " ADD CONSTRAINT fk_",col, " FOREIGN KEY(",col,") REFERENCES ",tablenameref, "(",colref,");")
    print(query)
    ##dbSendQuery(con, query)
    print("doneAddingForeign")
}

createIndexShort <- function(con, tableName, colsToIndex, shortname)
{
    astring = paste0("CREATE INDEX ix", shortname, " ON ", tableName, "(", paste(colsToIndex, collapse=","),")") 
    print(astring)
    dbSendQuery(con, astring)
    dbCommit(con)
    return(astring)
}

rebuild_ccdb <- function(db, rilHaplotypeProbsDir, founders, karyotype, host, user, password) 
{
    ##engine = "memory"
    print("connecting to db")
    con = dbConnect(MySQL(), user=user, password =password, dbname=db, host=host)
    dbSendQuery(con, 'set @@max_heap_table_size=4294967296;')
    print("about to begin builing il founder")
    pracma::tic()	
    buildInbredLineFounders(con = con, engine = engine, rilHaplotypeProbsDir = rilHaplotypeProbsDir, founders = founders, karyotype = karyotype)
    print("added all inbred line founders to db")
    pracma::toc()
    pracma::tic()
    buildInbredAlleles_nomax(con = con)
    print("joined alleles")
    pracma::toc()
    
    pracma::tic()
    buildInbredAlleles_onlyMax(con = con)
    print("got max diplotypes")
    pracma::toc() 
    pracma::tic()
    mergeMaxAndNoMax(con = con)
    print("merged in max boolean")
    
    pracma::toc()
    
    pracma::tic()
    hap_prefix = "allele_index"
    tablename  = "genotype"
    buildHaplotype(con, hap_prefix, tablename)
    print("built allele_hap")
    pracma::toc()
    
    pracma::tic()
    
##    buildInbredLineFixed(con = con, engine = engine)
	
    pracma::toc()

    hap_prefix = "founder_id"
    tablename = "diplotype"
    buildHaplotype(con, hap_prefix, tablename)
    pracma::toc()

    dbDisconnect(con)
    print("added all cc to db!")
}
