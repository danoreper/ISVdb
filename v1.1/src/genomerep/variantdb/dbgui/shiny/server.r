library(shiny)
library(pracma)
library(RMySQL)
con <- dbConnect(RMySQL::MySQL(), dbname = "exon_1410b", username='root',password='FrozenMixedVegetables',host="127.0.0.1",port=3308)

gene_list<-read.csv(file="genelist.csv",stringsAsFactors = FALSE)
gene_list_all<-as.character(c(gene_list$gene_id,gene_list$gene_id));names(gene_list_all)<-c(gene_list$gene_name,gene_list$mgi_symbol)
gene_list_all<-gene_list_all[-which(is.na(names(gene_list_all))==TRUE)]
gene_list<-as.list(gene_list_all)
gene_list_names<-names(gene_list)

source(file = "sql_shinny.R")

#source(file = "/Users/Yan2015/Dropbox/UNC/ValdarLab/Shinny/test/database.r")

strain_list<-read.csv(file="ccstrainmapid.csv")
strain_all<-as.character(strain_list$strain_id); names(strain_all)<-strain_list$strain_name
strain_list_new<-as.list(strain_all)

shinyServer(function(input, output, session) {
  total.start <- Sys.time()

  strain_pair<-data.frame()
  strain_pair_update<-function(){ 
    strain1<-sort(strain_new[as.numeric(input$strainp1new)])
    strain2<-sort(strain_new[as.numeric(input$strainp2new)])
    strain_pair_add<-t(apply(merge(strain1,strain2),1,sort))
    strain_pair<-rbind(strain_pair,strain_pair_add)
    strain_pair<-unique(strain_pair)
    out<-paste(apply( strain_pair[ , c(1,2) ] , 1 , paste , collapse = "-" ),sep="\n")
    out
  } 
#  output$stainpair <-renderText({ as.character(strain_pair_update()) })
  
  table_next <- eventReactive(input$submitvalue, {
    db.start<-Sys.time() # record time
    query<-generateQuery(input) # main function to generate Mysql query
    query<-paste(query,"LIMIT",input$limnum,sep=" ") # add output limit 
    out<-dbGetQuery(con,query) 
    db.end<-Sys.time()
    
    output$time1 <-renderText({ as.character(round(c(db.end-db.start),6)) }) # output time
    out
    })

  observeEvent(input$gene_search_action, {
    gene_input<-input$genesearch
    gene_match<-gene_list_names[grep(gene_input, gene_list_names)]
    if (length(gene_match)>40){
      gene_match<-gene_match[1:39]
    }
    updateSelectInput(session,"gene_list",choices=gene_list_all[gene_match])
  })
  
  observeEvent(input$strainlistnew_founders, {
    updateSelectInput(session,"strainlistnew",choices = strain_list_new , selected = c("1","2","3","4","5","6","7","8"))
  })
  
  
  output$table <- renderDataTable({table_next()}, options = list(pageLength = 20, searching=TRUE))

  output$downloadData <-downloadHandler(
    filename = function() { paste("output",'.txt', sep='') },
    content = function(file) {
      output$filename <-renderText({ as.character(file) })
      query<-generateQuery(input)
      sysquery<-paste("ssh ywcai@valdardb.its.unc.edu \"mysql --user='root' --password='FrozenMixedVegetables' 
                      --database='exon_1410' -e '",query,";' \" >",file, sep=" ")
      message(sysquery)
      system(sysquery)
    })
  
  total.end <- Sys.time()
  output$time2 <-renderText({ as.character(round(c(total.end-total.start),6)) })
})
