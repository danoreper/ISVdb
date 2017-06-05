# Shinny
library(shiny)

shinyUI(fluidPage(
  
  titlePanel("F to C transfer:"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("num", label = h3("Degree F:"), value = 100),
      hr(),
      h3("Degree C:"),
      fluidRow(column(5, verbatimTextOutput("value"))),
      
      helpText("F=32+1.8*C")
      ),
    
    mainPanel(
      h2("Donload the datafile:"),
      downloadButton('downloadData', 'Download'))
  )
))