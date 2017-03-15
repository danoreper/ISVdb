library("RMySQL")
# Here requires loadParams.R
# ../config/defaultvdb.yaml

source("./loadParams.R")

##TODO pass in
db         = "exon_1410"

message(db)

tmpdir     = prop$tmpdir
dbhost     = prop$variantdb$host
dbuser     = prop$variantdb$user
dbpassword = prop$variantdb$password
dbengine   = prop$variantdb$engine

con = dbConnect(MySQL(), user = dbuser, password = dbpassword, dbname=db, host=dbhost)
dbSendQuery(con, 'set @@max_heap_table_size=4294967296;')

build_genotype_view <- function(con)
{
    tablename = "genotype_view"
    query = paste("Drop table if exists ", tablename)
    dbSendQuery(con, query)
    
    query = paste("create table ", tablename,
                  "select",
                  "v1.variant_id as variant_id, ",
                  "v1.chrom as chrom, ",
                  "v1.pos as pos, ",
                  "strain.strain_name as strain_name, ",
##                  "a1.allele_index as allele_index_1",
##                  "a2.allele_index as allele_index_2",
                  "a1.allele as allele_1, ",
                  "a2.allele as allele_2, ",
                  "prob, ", 
                  "is_max, ", 
                  "c1.consequence_name as consequence_1, ",
                  "c2.consequence_name as consequence_2, ",
                  "g.gene_name as gene_name, ",
                  "t.transcript_name as transcript_name ",
                  
                  "from (select * from variant where variant_id>0 and variant_id<=1000)as v1",

                  "inner join genotype ila ",
                  "ignore index (ixstrain_idvariant_idallele_index_1allele_index_2is_maxprob, ixvariant_idstrain_idallele_index_1allele_index_2is_maxprob) ",
                  "on v1.variant_id=ila.variant_id ",
                  "inner join allele a1 on a1.variant_id = ila.variant_id and a1.allele_index=ila.allele_index_1 ",
                  "inner join allele a2 on a2.variant_id = ila.variant_id and a2.allele_index=ila.allele_index_2 ",
                  "inner join strain on ila.strain_id=strain.strain_id ",
                  #           "inner join allele ref ignore index(allele_index) ",
##                  "on ref.variant_id=v1.variant_id and ref.allele_index=0 ",


                  "inner join varinfo vinfo1 on vinfo1.variant_id = v1.variant_id and vinfo1.allele_index=ila.allele_index_1",
                  "inner join consequence c1 on c1.consequence_id = vinfo1.consequence_id",
                  "inner join gene g on vinfo1.gene_id = g.gene_id",
                  "inner join transcript t on vinfo1.transcript_id = t.transcript_id",

                  "inner join varinfo vinfo2 on vinfo2.variant_id = v1.variant_id and vinfo2.allele_index=ila.allele_index_2 and vinfo1.transcript_id = vinfo2.transcript_id",
                  "inner join consequence c2 on c2.consequence_id = vinfo2.consequence_id ",
                  #"WHERE v1.variant_id<1000 ",
                  "ENGINE=",,
                  ";")
    
  message(query)
#    dbSendQuery(con, query)              
    
}
    
build_genotype_view(con)
