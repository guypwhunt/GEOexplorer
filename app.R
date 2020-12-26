# Load Packages
library(shiny)
library(DT)

# UPDATE HELPER FUNCTION TO RETURN DATA

# Load Data
geo2r_data <- readRDS("data/analysis-output.rds")

ui <- fluidPage(
  titlePanel("GEO2R Data Visualisation"),
  textInput("GEO accession code", "GEO accession code", "Please enter the GEO accession code"),
  dataTableOutput('myTable')
)

server <- function(input, output, session){
  output$myTable <- renderDataTable({
    geo2r_data
  })
}

shinyApp(ui, server)