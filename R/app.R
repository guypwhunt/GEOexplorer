#' A Function to Load the GEO2R Data Visualisation Shiny App
#'
#' This function loads the GEO2R Data Visualisation Shiny App. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of analysis and graphics to be performed/generated from the GEO2R gene expression data.
#' @keywords GEO
#' @export
#' @examples loadApp()
#' @import shiny plotly heatmaply shinyHeatmaply ggplot2 shinyBS shinyjs BiocManager DT
#' @author Guy Hunt
loadApp <- function() {
  # Load Packages
  library(shiny)
  library(plotly)
  library(heatmaply)
  library(shinyHeatmaply)
  library(ggplot2)
  library(shinyBS)
  library(shinyjs)
  library(DT)

  ui <- fluidPage(
    titlePanel("GEO Explorer"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of samples in a GEO Series to identify genes that are differentially expressed across experimental conditions. GEOexplorer extends GEO2R's functionalities by enabling a richer set of analysis and graphics to be performed/generated from the GEO2R gene expression data. The development of GEOexplorer was made possible because of the excellent code provided by GEO2R
(https://www.ncbi.nlm.nih.gov/geo/geo2r/)."),
    # Source the Side Bar UI Components
    sourceSideBarUi(),
    mainPanel(tabsetPanel(type = "tabs",
                          # Source the Dataset Information UI Components
                          sourceDatasetInformationUi(),
                          # Source the Exploratory Data Analysis UI Components
                          sourceExploratoryDataAnalysisUi(),
                          # Source the Differential Gene Expression UI Components
                          sourceDifferentialGeneExpressionAnalysisUi()
    )
    )
  )

  server <- function(input, output, session){
    # Source Server Components
    sourceServer(input, output, session)
  }

  shinyApp(ui, server)
}
