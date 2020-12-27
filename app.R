# Load Packages
library(shiny)
library(GEOquery)
library(limma)
library(umap)

ui <- fluidPage(
  titlePanel("GEO2R Data Visualisation"),
  helpText("Select a GEO accession code to examine the gene expression data.

        Information will be collected from GEO2R."),
  textInput("geo_accession_code", "GEO accession code", "GSE18384"),
  mainPanel(dataTableOutput('myTable'))
)

server <- function(input, output, session){
  
  data_input <- reactive({gset <- getGEO(input$geo_accession_code, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  ex <- exprs(gset)
  # log2 transform
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) }
  return(ex)})
  
  output$myTable <- renderDataTable({
    data_input()
  })
}

shinyApp(ui, server)