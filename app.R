# Load Packages
library(shiny)
library(DT)

# UPDATE HELPER FUNCTION TO RETURN DATA

# Load Data
geo2r_data <- readRDS("data/analysis-output.rds")

ui <- fluidPage(
  titlePanel("GEO2R Data Visualisation"),
  dataTableOutput('myTable')
)

server <- function(input, output, session){
  output$myTable <- renderDataTable({
    geo2r_data
  })
}

shinyApp(ui, server)