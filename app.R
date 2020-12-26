# Load Packages
library(shiny)
library(DT)

# Load Data
geo2r_data <- readRDS("data/analysis-output.rds")

ui <- fluidPage(
  dataTableOutput('myTable')
)

server <- function(input, output, session){
  output$myTable <- renderDataTable({
    geo2r_data
  })
}

shinyApp(ui, server)