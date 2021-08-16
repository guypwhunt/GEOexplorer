#' A Function to Load the GEOexplorer Shiny App
#'
#' This function loads the GEOexplorer Shiny App.
#' The GEOexplorer Shiny App extends GEO2R's functionalities
#' by enabling a richer set of analysis and graphics to be
#' performed/generated from the gene expression data.
#' @keywords GEO
#' @export
#' @examples
#' app <- loadApp()
#' @import shinyHeatmaply Biobase
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @rawNamespace import(ggplot2, except = c(last_plot))
#' @importFrom grDevices palette
#' @importFrom graphics abline boxplot hist par text
#' @importFrom stats complete.cases cor model.matrix
#' na.omit prcomp princomp quantile
#' @importFrom vembedr embed_url
#' @author Guy Hunt
#' @return Large Shiny App
loadApp <- function() {
  ui <- fluidPage(
    titlePanel("GEOexplorer"),
    helpText(
      "GEO2R is an interactive web tool that allows
    users to compare two or more groups of experimental
    conditions in a GEO Series to identify genes that
    are differentially expressed across experimental
    conditions.
    GEOexplorer extends GEO2R's functionalities by enabling
    a richer set of analysis and graphics to be
    performed/generated from the GEO2R gene expression data.
    The development of GEOexplorer was made possible
    because of the excellent code provided by GEO2R
             (https://www.ncbi.nlm.nih.gov/geo/geo2r/)."
    ),

    # Source the Side Bar UI Components
    sourceSideBarUi(),
    mainPanel(
      tabsetPanel(
        type = "tabs",
        # Source the Dataset Information UI Components
        sourceDatasetInformationUi(),
        # Source the Exploratory Data Analysis UI Components
        sourceExploratoryDataAnalysisUi(),
        # Source the Differential Gene Expression UI Components
        sourceDifferentialGeneExpressionAnalysisUi()
      )
    )
  )

  server <- function(input, output, session) {
    # Source Server Components
    sourceServer(input, output, session)
  }

  shinyApp(ui, server)
}
