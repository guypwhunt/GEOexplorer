  # Load Packages
  library(shiny)
  library(GEOquery)
  library(limma)
  library(umap)
  library("maptools")
  
  ui <- fluidPage(
    titlePanel("GEO2R Data Visualisation"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of Samples in a GEO Series in order to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of graphics to be generated from the GEO2R outputs."),
    helpText("Select a GEO accession code to examine the gene expression data."),
    textInput("geo_accession_code", "GEO accession code", "GSE18384"),
    helpText("Please enter your platform of interest, if the Series is associated with multiple platforms."),
    textInput("platform", "Platform", "GPL6246"),
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Options", plotOutput("plot")),
                          tabPanel("Dataset", dataTableOutput('myTable')),
                          tabPanel("Box-and-Whisper Plot", plotOutput('boxPlot')),
                          tabPanel("Expression Density Plot", plotOutput('expressionDensity')),
                          tabPanel("Mean-Variance Plot", plotOutput('meanVariance')),
                          tabPanel("UMAP Plot", plotOutput('umap'))
    ))
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
    
    gset <- reactive({gset <- getGEO(input$geo_accession_code, GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep(input$platform, attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    return(gset)})
    
    # Data Set Plot
    output$myTable <- renderDataTable({
      data_input()
    })
    
    # Box-and-Whisker Plot
    output$boxPlot <- renderPlot({
      par(mar=c(7,4,2,1))
      title <- paste (input$geo_accession_code, "/", input$platform, sep ="")
      boxplot(data_input(), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
    })
    
    # Expression Value Distribution Plot
    output$expressionDensity <- renderPlot({
    par(mar=c(4,4,2,1))
    title <- paste (input$geo_accession_code, "/", input$platform, " value distribution", sep ="")
    plotDensities(data_input(), main=title, legend=F)
    })
    
    # Mean-Variance Plot
    output$meanVariance <- renderPlot({
    ex <- na.omit(data_input()) # eliminate rows with NAs
    plotSA(lmFit(ex), main= paste("Mean variance trend,", input$geo_accession_code))
    })
    
    # UMAP plot (multi-dimensional scaling)
    output$umap <- renderPlot({
    ex <- data_input()[!duplicated(data_input()), ]  # remove duplicates
    ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
    plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", pch=20, cex=1.5)
    pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
    })
  }
  
  shinyApp(ui, server)