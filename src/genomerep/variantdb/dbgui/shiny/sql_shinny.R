input_ex1<-list(view=1,
                selectCollum1=c("variant_id","chrom","pos","strain_name","allele","prob","is_max","consequence","gene_name","transcript_name"),
                chr="X",
                strnum=31000000,
                endnum=32000000,
                prob_cutoff="0.5",
                ismax=1,
                homo="Homo",
                variantfile="",
                strainlistnew=c(1,2,3,4),
                variant_id_text="1,2",
                limnum=1000)

# 1797674
input_ex2<-list(view=3,
                selectCollum3=c("variant_id","chrom","pos","gene_name","strain_name","allele","prob"),
                chr="ALL",
                strnum=31000000,
                endnum=32000000,
                prob_cutoff="0",
                ismax=1,
                homo="ALL",
                variantfile="",
                strainp1new=c(1),
                strainp2new=c(2),
                limnum=1000)
# Don't de;et the code 

# 4.14 Version
generateQuery<-function(input){
  
  if (input$variantinputtype == 1){
    variant_list<- as.numeric(strsplit(input$variant_id_text,",")[[1]])
  }
  if (input$variantinputtype == 2){
    strnum<-input$strnum;endnum<-input$endnum;chr<-input$chr;
  }
  if (input$variantinputtype == 3){
    gene_list<-input$gene_list
  }
  
  limnum<-input$limnum;prob_cutoff<-input$prob_cutoff;
  ismax<-input$ismax;strainlist<-input$strainlistnew 
  homo<-input$homo;consequence_list<- input$consequence
  
  # For each kind of view, build the corresponding FROM, SELECT, WHERE
  if (input$view==1){ # genotype view 
    selectCollum<-input$selectCollum1
    WHERE_list<-c();
    
    FROM<-"genotype_transcript_view"

    if ("consequence" %in% selectCollum) { 
      selectCollum<-c(selectCollum[-which(selectCollum=="consequence")],c("consequence_1","consequence_2"))}
    if ("allele" %in% selectCollum) { 
      selectCollum<-c(selectCollum[-which(selectCollum=="allele")],c("allele_1","allele_2"))}
    
    if (homo=="Homo"){ WHERE_list<-c(WHERE_list,paste("allele_id_1 <=> allele_id_2")) }
    if (homo=="Hetero"){ WHERE_list<-c(WHERE_list,paste("NOT allele_id_1 <=> allele_id_2"))}
    
    if (input$variantinputtype == 1){
      WHERE_list<-c(WHERE_list,paste("(",paste(paste("variant_id=",variant_list,sep=""),collapse=" OR "),")",sep=""))}
    
    if (input$variantinputtype == 2){
      if (chr=="ALL"){}else{  WHERE_list<-c(WHERE_list,paste("chrom='",chr,"'",sep="")) }
      if (!is.na(strnum) ){WHERE_list<-c(WHERE_list,paste("pos>=",strnum,sep=""))}
      if (!is.na(endnum) ){WHERE_list<-c(WHERE_list,paste("pos<=",endnum,sep=""))}
    }
    
    if (input$variantinputtype == 3){
      WHERE_list<-c(WHERE_list,paste("(",paste(paste("gene_id=",gene_list,sep=""),collapse=" OR "),")",sep=""))
    }
    
    if (prob_cutoff>0){WHERE_list<-c(WHERE_list,paste("prob>",prob_cutoff,sep=""))}
    if (ismax) { WHERE_list<-c(WHERE_list, "is_max=1")}
    if (length(strainlist)>0){
      WHERE_list<-c(WHERE_list,paste("(",paste(paste("strain_id=",strainlist,sep=""),collapse=" OR "),")",sep="")) 
    }
    if (length(consequence_list)>0){
      WHERE_list<-c(WHERE_list,paste("(",paste(paste(" consequence_id_1= ",consequence_list," or consequence_id_2= ",
                                                     consequence_list,sep=""),collapse=" OR "),")",sep="")) 
    }
    
    if (length(WHERE_list)>0){
      WHERE<-paste(WHERE_list,collapse=" and ")
      SELECT<-paste(selectCollum,collapse=",")
      query<-paste("SELECT distinct ",SELECT,"FROM",FROM,"WHERE",WHERE,sep=" ")
    }else{
      SELECT<-paste(selectCollum,collapse=",")
      query<-paste("SELECT distinct ",SELECT,"FROM",FROM,sep=" ")
    }
  }
##########################################################################################  
  if (input$view==2){ # diplotype_view
    selectCollum<-input$selectCollum2
    WHERE_list<-c();
    FROM<-"diplotype_view"
    
    if ("founder_name" %in% selectCollum) { selectCollum<-c(selectCollum[-which(selectCollum=="founder_name")],c("founder_name_1","founder_name_2"))}
    
    if (homo=="Homo"){ WHERE_list<-c(WHERE_list,paste("founder_name_1 <=> founder_name_2")) }
    if (homo=="Hetero"){ WHERE_list<-c(WHERE_list,paste("NOT founder_name_1 <=> founder_name_2"))}
    
    if (chr=="ALL"){}else{  WHERE_list<-c(WHERE_list,paste("chrom='",chr,"'",sep="")) }
    if (!is.na(strnum) ){WHERE_list<-c(WHERE_list,paste("pos>=",strnum,sep=""))}
    if (!is.na(endnum) ){WHERE_list<-c(WHERE_list,paste("pos<=",endnum,sep=""))}
    if (prob_cutoff>0){WHERE_list<-c(WHERE_list,paste("prob>",prob_cutoff,sep=""))}
  #  if (ismax) { WHERE_list<-c(WHERE_list, "is_max=1")}
    if (length(strainlist)>0){
      WHERE_list<-c(WHERE_list,paste("(",paste(paste("strain_id=",strainlist,sep=""),collapse=" OR "),")",sep="")) # Change Later
    }
    
    if (length(WHERE_list)>0){
      WHERE<-paste(WHERE_list,collapse=" and ")
      SELECT<-paste(selectCollum,collapse=",")
      query<-paste("SELECT distinct ",SELECT,"FROM",FROM,"WHERE",WHERE,sep=" ")
    }else{
      SELECT<-paste(selectCollum,collapse=",")
      query<-paste("SELECT distinct ",SELECT,"FROM",FROM,sep=" ")
    }
    
  }
  
###################################################################################
  
  if (input$view==3){ # Genotype_sampling view
    selectCollum<-input$selectCollum3
    message(selectCollum)
    selectTrans<-c("s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos","s1.gene_name as gene_name",
                   "s1.strain_name as strain_name_s1, s2.strain_name as strain_name_s2",
                   "s1.allele as allele_s1, s2.allele as allele_s2","s1.prob*s2.prob as prob","s1.transcript_name as transcript_name")
    names(selectTrans)<-c("variant_id","chrom","pos","gene_name","strain_name","allele","prob","transcript_name")
    SELECT<-paste(selectTrans[selectCollum],collapse=",")
    FROM<- "genotype_sampling_view as s1 inner join genotype_sampling_view as s2 on s1.variant_id=s2.variant_id and s1.transcript_id=s2.transcript_id";
    
    strain_pair_matrix<-t(apply(merge(input$strainp1new,input$strainp2new),1,sort))
    strain_pair<-data.frame(p1=strain_pair_matrix[,1],p2=strain_pair_matrix[,2])
    message(strain_pair)
    WHERE_list<-c();
    
    if (homo=="Homo"){ WHERE_list<-c(WHERE_list,paste("s1.allele_index <=> s2.allele_index")) }
    if (homo=="Hetero"){ WHERE_list<-c(WHERE_list,paste("NOT s1.allele_index <=> s2.allele_index"))}
    if (prob_cutoff>0){WHERE_list<-c(WHERE_list,paste("s1.prob*s2.prob >",prob_cutoff,sep=""))}
    
    
    if (chr=="ALL"){}else{  
      WHERE_list<-c(WHERE_list,paste("s1.chrom='",chr,"'",sep="")) 
      WHERE_list<-c(WHERE_list,paste("s2.chrom='",chr,"'",sep=""))    }
    if (!is.na(strnum) ){
      WHERE_list<-c(WHERE_list,paste("s1.pos>=",strnum,sep=""))
      WHERE_list<-c(WHERE_list,paste("s2.pos>=",strnum,sep=""))   }
    if (!is.na(endnum) ){
      WHERE_list<-c(WHERE_list,paste("s1.pos<=",endnum,sep=""))
      WHERE_list<-c(WHERE_list,paste("s2.pos<=",endnum,sep=""))   }
    
    WHERE_list_strain<-c()
    for (i in c(1:length(strain_pair[[1]]))){
      WHERE_list_strain<-c(WHERE_list_strain,paste("(s1.strain_id=",strain_pair$p1[i]," AND s2.strain_id=",strain_pair$p2[i],")", sep=""))
    }
    WHERE_list<-c(WHERE_list,paste("(",paste(WHERE_list_strain,collapse = " OR "),")",sep=""))
    
    WHERE<-paste(WHERE_list,collapse=" and ")
    query<-paste("SELECT distinct ",SELECT,"FROM",FROM,"WHERE",WHERE,
                 "GROUP BY variant_id,s1.allele_index,s2.allele_index,s1.strain_id,s2.strain_id,s1.transcript_id",sep=" ")
  }
  
#############################################################################################  
  
  if (input$view==4){ # Genotype_sampling view
    selectCollum<-input$selectCollum4
    
    selectTrans<-c("s1.variant_id as variant_id","s1.chrom as chrom","s1.pos as pos","s1.gene_name as gene_name","s1.strain_name as strain_name_s1, s2.strain_name as strain_name_s2","s1.founder_name as founder_s1, s2.founder_name as founder_s2","s1.prob*s2.prob as prob")
    names(selectTrans)<-c("variant_id","chrom","pos","gene_name","strain_name","founder_name","prob")
    SELECT<-paste(selectTrans[selectCollum],collapse=",")
    FROM<- "diplotype_sampling_view as s1 inner join diplotype_sampling_view as s2 on s1.variant_id=s2.variant_id"
    
   # if (input$strainVersion){
  #    strain_pair_matrix<-t(apply(merge(input$strainp1old,input$strainp2old),1,sort))
  #  }else{
      strain_pair_matrix<-t(apply(merge(input$strainp1new,input$strainp2new),1,sort))
  #  }
    strain_pair<-data.frame(p1=strain_pair_matrix[,1],p2=strain_pair_matrix[,2])
    
    message(strain_pair)
    WHERE_list<-c();
    
    if (homo=="Homo"){ WHERE_list<-c(WHERE_list,paste("s1.founder_id <=> s2.founder_id")) }
    if (homo=="Hetero"){ WHERE_list<-c(WHERE_list,paste("NOT s1.founder_id <=> s2.founder_id"))}
    if (prob_cutoff>0){WHERE_list<-c(WHERE_list,paste("s1.prob*s2.prob >",prob_cutoff,sep=""))}
    
    if (chr=="ALL"){}else{  
      WHERE_list<-c(WHERE_list,paste("s1.chrom='",chr,"'",sep="")) 
      WHERE_list<-c(WHERE_list,paste("s2.chrom='",chr,"'",sep="")) 
    }
    if (!is.na(strnum) ){
      WHERE_list<-c(WHERE_list,paste("s1.pos>=",strnum,sep=""))
      WHERE_list<-c(WHERE_list,paste("s2.pos>=",strnum,sep=""))
    }
    if (!is.na(endnum) ){
      WHERE_list<-c(WHERE_list,paste("s1.pos<=",endnum,sep=""))
      WHERE_list<-c(WHERE_list,paste("s2.pos<=",endnum,sep=""))
    }
    
    WHERE_list_strain<-c()
    for (i in c(1:length(strain_pair[[1]]))){
      WHERE_list_strain<-c(WHERE_list_strain,paste("(s1.strain_id=",strain_pair$p1[i]," AND s2.strain_id=",strain_pair$p2[i],")", sep=""))
    }
    WHERE_list<-c(WHERE_list,paste("(",paste(WHERE_list_strain,collapse = " OR "),")",sep=""))
    
    WHERE<-paste(WHERE_list,collapse=" and ")
    query<-paste("SELECT ",SELECT,"FROM",FROM,"WHERE",WHERE,sep=" ")
    
  }
  message(query)
  query
}

#generateQuery(input_ex1)
#generateQuery(input_ex2)
