library(igraph)
############# Table Structure ###########
# variant V     variant_id   chrom  pos
# strian S strain_id strain_name
# diplotype D strain_id variant_id founder_id_1  founder_id_2
# genotype Gt strain_id variant_id allele_index_1 allele_index_2  prob  is_max
# consequence_info  Cf variant_id allele_index gene_id transcript_id consequence_id
# gene G gene_id  gene_name 
# allele A variant_id allele_index  allele
# transcript T transcript_id  transcript_name

fllnm <- c("variant","strian","diplotype","genotype","consequence_info","gene","allele","transcript")   # record the full name of tables
names(fllnm)<-c("V","S","D","Gt","Cf","G","A","T")

relation<-matrix(c(sort(c("V","Cf")),"V.variant_id","Cf.variant_id",
                   sort(c("V","D")),"V.variant_id","D.variant_id",
                   sort(c("V","Gt")),"V.variant_id","Gt.variant_id",
                   sort(c("V","A")),"V.variant_id","A.variant_id",
                   sort(c("S","D")),"S.strain_id","D.strain_id",
                   sort(c("S","Gt")),"S.strain_id","Gt.strain_id",
                   sort(c("D","Gt")),"D.strain_id","Gt.strain_id",
                   sort(c("D","Cf")),"D.variant_id","Cf.variant_id",
                   sort(c("D","A")),"D.variant_id","A.variant_id",
                   sort(c("Gt","Cf")),"Gt.variant_id","Cf.variant_id",
                   sort(c("Gt","Cf")),"Gt.allele_index_1","Cf.allele_index",
                   sort(c("Gt","Cf")),"Gt.allele_index_2","Cf.allele_index",
                   sort(c("Gt","A")),"Gt.variant_id","A.variant_id",
                   sort(c("Gt","A")),"Gt.allele_index_1","A.allele_indexd",
                   sort(c("Gt","A")),"Gt.allele_index_2","A.allele_index",
                   sort(c("Cf","G")),"Cf.gene_id","G.gene_id",
                   sort(c("Cf","A")),"Cf.variant_id","A.variant_id",
                   sort(c("Cf","A")),"Cf.allele_index","A.allele_index",
                   sort(c("Cf","T")),"Cf.transcript_id","Cf.transcript_id"
                 ),ncol=4,byrow=TRUE)
mygraph<-graph(c(t(relation[,c(1:2)])),directed = FALSE)
#plot(mygraph)
relation<-data.frame(relation)

SQLwriter <-function(VariantID=FALSE,Chr=TRUE,Pos=TRUE,Gene=TRUE,Allele=TRUE,Strain=TRUE) { 
  select <- c()
  if (VariantID){ select<-c(select,"V.variant_id")}
  if (Chr){ select<-c(select,"V.chrom")}
  if (Pos){ select<-c(select,"V.pos")}
  if (Gene){ select<-c(select,"G.gene_name")}
  if (Allele){ select<-c(select,"A.allele")}
  if (Strain){ select<-c(select,"S.strain_name")}
  select<-paste(select,collapse=",")
  select<-paste("SELECT",select,sep=" ")
  
  from<-"FROM"
  table_need<-c()
  # Find the table needed for the array
  if (VariantID){ table_need<-c(table_need,"V")}
  if (Chr){ table_need<-c(table_need,"V")}
  if (Pos){ table_need<-c(table_need,"V")}
  if (Gene){ table_need<-c(table_need,"G")}
  if (Allele){ table_need<-c(table_need,"A")}
  if (Strain){ table_need<-c(table_need,"S")}
  table_need<-unique(table_need)
  # Find the right way to join the tables
  # Get the original subgraph of the relation key
  # Remove unnessary keys
  
  
  where<-"WHERE"
  if (TRUE) {  }
  else{
  }
  return(paste(select,from,where,sep="\n")) 
}


