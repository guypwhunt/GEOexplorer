library(shiny)
library(rhandsontable)
ui <- fluidPage(
  titlePanel("Ttile"),
  sidebarLayout(
    sidebarPanel(
      actionButton("runButton","Change Dataframes")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("OldIrisTab", rHandsontableOutput('OldIris')),
        tabPanel("NewIrisTab", DT::dataTableOutput("NewIris"))
      ))))

server <- function(input,output,session)({
  values <- reactiveValues()
  output$OldIris <- renderRHandsontable({
    rhandsontable(iris)
  })

  observeEvent(input$runButton, {
    values$data <-  hot_to_r(input$OldIris)
  })

  output$NewIris <- DT::renderDataTable({
    values$data
  })

})
shinyApp(ui, server)
