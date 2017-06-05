# Shinny
library(shiny)

# Strainlist
strain_list<-read.csv(file="ccstrainmapid.csv")
strain_all<-as.character(strain_list$strain_id); names(strain_all)<-strain_list$strain_name
strain_list_new<-as.list(strain_all)
# Consequence list
consequence_list<-read.csv(file="consequence.csv")
consequence_all<-as.character(consequence_list$consequence_id);names(consequence_all)<-consequence_list$consequence_name
consequence_list<-as.list(consequence_all)

gene_list<-read.csv(file="genelist.csv",stringsAsFactors = FALSE)
gene_list_all<-as.character(c(gene_list$gene_id,gene_list$gene_id));names(gene_list_all)<-c(gene_list$gene_name,gene_list$mgi_symbol)
gene_list_all<-gene_list_all[-which(is.na(names(gene_list_all))==TRUE)]
gene_list<-as.list(gene_list_all)
gene_list_default<-as.list(gene_list_all[1:10])

genotype_view_choice<-c("variant_id","chrom","pos","strain_name","allele","prob","is_max","consequence","gene_name","transcript_name")
diplotype_view_choice<-c("variant_id","chrom","pos","gene_name","strain_name","founder_name","prob")
genotype_samplying_view_choice<-c("variant_id","chrom","pos","gene_name","strain_name","allele","prob","transcript_name")
diplotype_samplying_view_choice<-c("variant_id","chrom","pos","gene_name","strain_name","founder_name","prob")

shinyUI(fluidPage(
  
  titlePanel("Query variants in the region"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("variantinputtype", label = h4("Choose variant by:"), 
                  choices = list("Variant_id" = 1, "Region" = 2, "Gene" = 3), 
                  selected = 1),
      
      conditionalPanel(condition = "input.variantinputtype == 1",
                       fileInput('variantfile', 'Upload variant_id file',accept=c('text/plain', '.csv')),
                       textInput("variant_id_text", label = h4("Input specific variant_id (comma splict)"), value = "")
      ),
      
      conditionalPanel(condition = "input.variantinputtype == 2",
                       fileInput('bedfile', 'Upload bed file',accept=c('text/plain', '.bed')),
                       fluidRow(
                         column(3,selectInput("chr", label = h5("Chrom"), choices = c(1:19,"X","ALL"),selected = "ALL")),
                         column(4,numericInput("strnum", label = h5("Start site"), value = NULL,step=1000000)),
                         column(4,numericInput("endnum", label = h5("End site"), value = NULL,step=1000000)))
      ),
      
      conditionalPanel(condition = "input.variantinputtype == 3",
                       fileInput('genefile', 'Upload gene list file',accept=c('text/plain', '.txt')),
                       textInput("genesearch", label="Search for gene name for ensembl ID", value=""),
                       actionButton("gene_search_action", "Search Gene"),                       
                       selectInput("gene_list", h5('Genes search result'),selectize=TRUE, 
                                   choices=gene_list_default, multiple = TRUE)
      ),
      
      selectInput("view", label = h4("Choose view type"), 
                  choices = list("Genotype view" = 1, "Diplotype view" = 2, "Genotype cross view" = 3, "Diplotype cross view"=4), 
                  selected = 1),

      conditionalPanel(condition = "input.view == 1", 
                       selectizeInput('selectCollum1', h5('Select included fields'),
                                      choices = genotype_view_choice, selected=genotype_view_choice, multiple = TRUE)),
      conditionalPanel(condition = "input.view == 2", 
                       selectizeInput('selectCollum2', h5('Select included fields'),
                                      choices =diplotype_view_choice, selected=diplotype_view_choice, multiple = TRUE)),
      conditionalPanel(condition = "input.view == 3", 
                       selectizeInput('selectCollum3', h5('Select included fields'),
                                      choices = genotype_samplying_view_choice, selected=genotype_samplying_view_choice, multiple = TRUE)),
      conditionalPanel(condition = "input.view == 4", 
                       selectizeInput('selectCollum4', h5('Select included fields'),
                                      choices = diplotype_samplying_view_choice, selected=diplotype_samplying_view_choice, multiple = TRUE)),
      
      fluidRow(
        column(3,numericInput("prob_cutoff", label = h5("Prob >="), value = 0.00,step=0.1)),
        column(4,selectInput("homo", label = h5("Heterozygous?"), choices = c("Homo","Hetero","All"), selected = "All")),
        conditionalPanel(condition = "input.view == 1 | input.view == 2", 
                         column(3,checkboxInput("ismax", label = h5("Ismax=True?"), value = FALSE)))
        ),
      conditionalPanel(condition = "input.view == 1", 
                         selectInput("consequence", label = h5("Consequence? (default is ALL)"), choices = consequence_list,multiple = TRUE )),

      conditionalPanel(condition = "input.view == 3 | input.view == 4", 
                      column(5,selectInput("strainp1new", h5('Strains'),selectize=TRUE, choices=strain_list_new, multiple = TRUE)),
                      column(5,selectInput("strainp2new", h5('Strains'),selectize=TRUE, choices=strain_list_new, multiple = TRUE))
     ),
     
      conditionalPanel(condition = "input.view == 1 | input.view == 2", 
                       selectInput('strainlistnew', h5('Strains (Default is all)'), choices=strain_list_new, multiple = TRUE),
                       actionButton("strainlistnew_founders", "Set founders")   ),
     
      fileInput('strainfile', 'Or Upload strain name file',accept=c('text/plain', '.csv')),
      
      actionButton("submitvalue", label ="Show Tables"),    
      downloadButton('downloadData', 'Download to file'),
  #    textInput("filename", label = h6("Outputfilename.format"), value = "out"),
      numericInput("limnum", label = h6("Outputlinelimit"), value = 1000),
      
      fluidRow(column(5, verbatimTextOutput("filename"))),
  
      h6("MySQL running time is:"),
      fluidRow(column(5, verbatimTextOutput("time1"))),
      h6("Gui running time is:"),
      fluidRow(column(5, verbatimTextOutput("time2"))),
      h6()
      ),
    
    mainPanel(
      h2("Table output:"),
      fluidRow(dataTableOutput("table"))
    )
  )
))