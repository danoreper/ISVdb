library("RMySQL")

##TODO pass in
source("./loadParams.R")

db         = "exon_1410"

tmpdir     = prop$tmpdir
dbhost     = prop$variantdb$host
dbuser     = prop$variantdb$user
dbpassword = prop$variantdb$password
dbengine   = prop$variantdb$engine

con = dbConnect(MySQL(), user = dbuser, password = dbpassword, dbname=db, host=dbhost)
dbSendQuery(con, 'set @@max_heap_table_size=4294967296;')

dbSendQuery(con, paste('SET default_storage_engine=',dbengine,";",sep=""))

source("./genomerep/build_view_query.R")

# Put the where, all tables contains variant to the where.

build_diplotype_sampling_view <- function(con)
{
    maxid<-as.numeric(dbGetQuery(con,'select max(variant_id) from variant;'))
    sep<-100000
    lowerbound<-1+c(0:floor(maxid/sep))*sep

    tablename = "diplotype_sampling_view"     #"genotype_view"
    query = paste("Drop table if exists ", tablename)
    dbSendQuery(con, query)

    for (i in c(1:length(lowerbound))){
        if (i==1){
		query=paste("create table ", tablename,select_query_diplotype_sampling(1,sep))
		dbSendQuery(con, query)
        }else{
		query=paste("insert into", tablename, select_query_diplotype_sampling(lowerbound[i],lowerbound[i]+sep-1))
		dbSendQuery(con, query)
	}
    }
    query = "CREATE INDEX pIndexDS on diplotype_sampling_view (chrom,pos)"
    dbSendQuery(con, query)
    query = "ALTER TABLE diplotype_sampling_view PARTITION BY KEY(strain_id) PARTITIONS 10";
    dbSendQuery(con, query)
}
system.time(build_diplotype_sampling_view(con))
#"CREATE INDEX pIndex on genotype_view (chrom,pos);"


##############################
build_genotype_sampling_view <- function(con)
{
    maxid<-as.numeric(dbGetQuery(con,'select max(variant_id) from variant;'))
    sep<-100000
    lowerbound<-1+c(0:floor(maxid/sep))*sep

    tablename = "genotype_sampling_view"     #"genotype_view"
    query = paste("Drop table if exists ", tablename)
    dbSendQuery(con, query)

    for (i in c(1:length(lowerbound))){
        if (i==1){
                query=paste("create table ", tablename,select_query_genotype_sampling(1,sep))
                dbSendQuery(con, query)
        }else{
                query=paste("insert into", tablename, select_query_genotype_sampling(lowerbound[i],lowerbound[i]+sep-1))
                dbSendQuery(con, query)
        }
    }
    query = "CREATE INDEX pIndexGS on genotype_sampling_view (chrom,pos)"
    dbSendQuery(con, query)
    query = "ALTER TABLE genotype_sampling_view PARTITION BY KEY(strain_id) PARTITIONS 10";
    dbSendQuery(con, query)
}
system.time(build_genotype_sampling_view(con))

build_diplotype_view <- function(con)
{
    maxid<-as.numeric(dbGetQuery(con,'select max(variant_id) from variant;'))
    sep<-100000
    lowerbound<-1+c(0:floor(maxid/sep))*sep

    tablename = "diplotype_view"     #"genotype_view"
    query = paste("Drop table if exists ", tablename)
    dbSendQuery(con, query)

    for (i in c(1:length(lowerbound))){
        if (i==1){
                query=paste("create table ", tablename,select_query_diplotype(1,sep))
                dbSendQuery(con, query)
        }else{
                query=paste("insert into", tablename, select_query_diplotype(lowerbound[i],lowerbound[i]+sep-1))
                dbSendQuery(con, query)
        }
    }
    query = "CREATE INDEX pIndexD on diplotype_view (chrom,pos)"
    dbSendQuery(con, query)
}
system.time(build_diplotype_view(con))

build_genotype_view <- function(con)
{
    maxid<-as.numeric(dbGetQuery(con,'select max(variant_id) from variant;'))
    sep<-100000
    lowerbound<-1+c(0:floor(maxid/sep))*sep

    tablename = "genotype_transcript_view"     #"genotype_view"
    query = paste("Drop table if exists ", tablename)
    dbSendQuery(con, query)

    for (i in c(1:length(lowerbound))){
        if (i==1){
                query=paste("create table ", tablename,select_query_genotype(1,sep))
                dbSendQuery(con, query)
        }else{
                query=paste("insert into", tablename, select_query_genotype(lowerbound[i],lowerbound[i]+sep-1))
                dbSendQuery(con, query)
        }
    }
    query = "CREATE INDEX pIndexGT on genotype_view (chrom,pos)"
    dbSendQuery(con, query)
}
system.time(build_genotype_view(con))


#> system.time(build_genotype_view(con))
#   user  system elapsed
#  0.005   0.000 357.650
