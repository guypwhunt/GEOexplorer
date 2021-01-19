  # Load Packages
  library(shiny)
  library(GEOquery)
  library(limma)
  library(umap)
  library("maptools")
  library(ggplot2)
  source("geoIntegration.R")
  source("geo2rDataVisualisation.R")
  
  ui <- fluidPage(
    titlePanel("GEO2R Data Visualisation"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of Samples in a GEO Series in order to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of graphics to be generated from the GEO2R outputs."),
    sidebarPanel(
      helpText("Input a GEO accession code to examine the gene expression data."),
      textInput("geoAccessionCode", "GEO accession code", "GSE18384"),
      helpText("Input the platform of interest, if the series is associated with multiple platforms."),
      textInput("platform", "Platform", "GPL6246"),
      selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
                  choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Hochberg", "Holm", "Hommel", "None")),
      radioButtons("logTransformation",
                   label="Apply log transformation to the data:",
                   choices=list("Auto-Detect","Yes","No"),
                   selected="Auto-Detect"),
      # Need to get the radio button below working or delete
      radioButtons("limmaPrecisionWeights",
                   label="Apply limma precision weights (vooma):",
                   choices=list("Yes","No"),
                   selected="No"), 
      # Need to get the radio button below working or delete
      radioButtons("forceNormalization",
                   label="Force normalization:",
                   choices=list("Yes","No"),
                   selected="No"), 
      # Need to get the radio button below working or delete
      radioButtons("forceNormalization",
                   label="Category of Platform annotation to display on results:",
                   choices=list("Submitter supplied","NCBI generated"),
                   selected="No"), 
      br(),
      helpText("Plot displays"),
      sliderInput("integer", "Significance level cut-off:",
                  min = 0, max = 1,
                  value = 0.05),
      # add grouping functionality
    ),
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Dataset", dataTableOutput('myTable')),
                          tabPanel("GEO2R Data Visualization",
                          tabsetPanel(type = "tabs",
                          tabPanel("Box-and-Whisper Plot", plotOutput('boxPlot')),
                          tabPanel("Expression Density Plot", plotOutput('expressionDensity')),
                          tabPanel("Mean-Variance Plot", plotOutput('meanVariance')),
                          tabPanel("UMAP Plot", plotOutput('umap')))),
                          tabPanel("Next Generation Data Visualization")
    ))
  )
  
  server <- function(input, output, session){
    
    gsetData <- reactive({getGeoData(input$geoAccessionCode, input$platform)})
    
    dataInput <- reactive({extractGeoData(gsetData(), input$logTransformation)
    })
    
    
    # Data Set Plot
    output$myTable <- renderDataTable({
      dataInput()
    })
    
    # Box-and-Whisker Plot
    output$boxPlot <- renderPlot({
      boxAndWhiskerPlot(input$geoAccessionCode, input$platform, dataInput())
      })
    
    # Expression Value Distribution Plot
    output$expressionDensity <- renderPlot({
      expressionValueDistributionPlot(input$geoAccessionCode, input$platform, dataInput())
    })
    
    # Mean-Variance Plot
    output$meanVariance <- renderPlot({
      meanVariancePlot(input$geoAccessionCode, input$platform, dataInput())
    })
    
    # UMAP plot (multi-dimensional scaling)
    output$umap <- renderPlot({
      umapPlot(input$geoAccessionCode, input$platform, dataInput())
      })
  }
  
  shinyApp(ui, server)