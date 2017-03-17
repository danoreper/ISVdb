library(shiny)

shinyServer(function(input, output) {
  
  # You can access the value of the widget with input$num, e.g.
  output$value <- renderText({    round((input$num-32)/1.8,3)   })
  
  output$downloadData <-downloadHandler(
    filename = function() { paste('out', '.table', sep='') },
    content = function(file) {
      write.table(data.frame(F=input$num,C=round((input$num-32)/1.8,3)),
                  sep="\t",quote=FALSE,row.names=FALSE, file)
    }
  )
})