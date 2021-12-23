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
#' @importFrom shinybusy add_busy_spinner
#' @author Guy Hunt
#' @return Large Shiny App
loadApp <- function() {
  ui <- fluidPage(
    add_busy_spinner(spin = "fading-circle",
                     height = "100px",
                     width = "100px"),
    sourceUi()
  )

  server <- function(input, output, session) {
    # Source Server Components
    sourceServer(input, output, session)
  }

  shinyApp(ui, server)
}
