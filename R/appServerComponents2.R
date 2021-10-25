#' A Function to Return the Server Component
#'
#' A Function to Return the Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceServer()
#' @importFrom DT renderDataTable JS
#' @importFrom shinyBS addTooltip
#' @importFrom utils write.csv
#' @author Guy Hunt
#' @noRd
sourceServer2 <- function(input, output, session) {
  datasetInformationServer <- (

    output$table <- renderTable({

      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.

      req(input$file1)

      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {
          df <- read.csv(input$file1$datapath)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )

      return(df)

    })

  )
  return(datasetInformationServer)
}
