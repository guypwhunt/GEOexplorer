# Load Packages
library(shiny)
library(DT)

# UPDATE HELPER FUNCTION TO RETURN DATA
# LOOK AT LESSON 6

# Load Data
geo2r_data <- readRDS("data/analysis-output.rds")

ui <- fluidPage(
  titlePanel("GEO2R Data Visualisation"),
  helpText("Select a GEO accession code to examine the gene expression data.

        Information will be collected from GEO2R."),
  textInput("GEO accession code", "GEO accession code", "Please enter the GEO accession code"),
  mainPanel(dataTableOutput('myTable'))
)

server <- function(input, output, session){
  output$myTable <- renderDataTable({
    geo2r_data
  })
}

shinyApp(ui, server)