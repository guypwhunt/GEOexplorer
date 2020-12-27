  # Load Packages
  library(shiny)
  library(GEOquery)
  library(limma)
  library(umap)
  
  ui <- fluidPage(
    titlePanel("GEO2R Data Visualisation"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of Samples in a GEO Series in order to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of graphics to be generated from the GEO2R outputs."),
    helpText("Select a GEO accession code to examine the gene expression data."),
    textInput("geo_accession_code", "GEO accession code", "GSE18384"),
    helpText("Please enter your platform of interest, if the Series is associated with multiple platforms."),
    textInput("platform", "Platform", "GPL6246"),
    mainPanel(dataTableOutput('myTable'))
  )
  
  server <- function(input, output, session){
    
    data_input <- reactive({gset <- getGEO(input$geo_accession_code, GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep(input$platform, attr(gset, "names")) else idx <- 1
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