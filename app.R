  # Load Packages
  library(shiny)
  library(GEOquery)
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  source("geoIntegration.R")
  source("geo2rDataVisualisation.R")
  source("analyticsFunctions.R")
  
  # Required for PCA data visualizations
  library(factoextra)
  
  
  ui <- fluidPage(
    titlePanel("GEO2R Data Visualisation"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of Samples in a GEO Series in order to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of graphics to be generated from the GEO2R outputs."),
    sidebarPanel(
      helpText("Input a GEO accession code to examine the gene expression data."),
      textInput("geoAccessionCode", "GEO accession code", "GSE18380"),
      helpText("Input the platform of interest, if the series is associated with multiple platforms."),
      textInput("platform", "Platform", "GPL4694"),
      radioButtons("logTransformation",
                   label="Apply log transformation to the data:",
                   choices=list("Auto-Detect","Yes","No"),
                   selected="Auto-Detect"),
      radioButtons("knnTransformation",
                   label="Apply k-nearest neighbors (KNN) algorithm to predict null data:",
                   choices=list("Yes","No"),
                   selected="No"),
      helpText("Rows with over 50% missing values are imputed using the overall mean per sample. Columns with over 80% will cause an error in the KNN computation."),
      br(),
      helpText("The Buttons Below here currently are not functional!!!"),
      # Need to get the radio button below working or delete
      selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
                  choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Hochberg", "Holm", "Hommel", "None")),
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
                          tabPanel("Next Generation Data Visualization",
                                   tabsetPanel(type = "tabs",
                                               tabPanel("Principal Component Analysis", tabsetPanel(type = "tabs",
                                                                                      tabPanel("Scree Plot", plotOutput('pcaScreePlot')),
                                                                                      tabPanel("Individuals Plot", plotOutput('pcaIndividualsPlot')),
                                                                                      tabPanel("Variables Plot", plotOutput('pcaVariablesPlot')),
                                                                                      tabPanel("Individuals and Variables Biplot", plotOutput('pcaBiplotPlot'))
                                                                                      )))
    )
    )
  )
  )
  
  server <- function(input, output, session){
    
    # Get GEO2R data
    gsetData <- reactive({getGeoData(input$geoAccessionCode, input$platform)})
    
    # Extract expression data
    expressionData <- reactive({extractExpressionData(gsetData())
    })
    
    # Apply log transformation to expression data if necessary
    dataInput <- reactive({logTransformExpressionData(expressionData(), input$logTransformation)
    })
    
    # Perform KNN transformation on log expression data if necessary
    knnDataInput <- reactive({knnDataTransformation(dataInput(), input$knnTransformation)
    })
    
    # Perform KNN transformation on log expression data for PCA input
    knnDataInputForPca <- reactive({knnDataTransformation(dataInput(), "Yes")
    })
    
    # Trigger this in the back ground
    # Move this into the analytics file
    # Perform PCA analysis on KNN transformation expression data
    pcaDataInput <- reactive({pca <- prcomp(knnDataTransformation(knnDataInputForPca(),"Yes"), scale = TRUE)
    return(pca)
    })
    
    # Data Set Plot
    output$myTable <- renderDataTable({
      knnDataInput()
    })
    
    # Box-and-Whisker Plot
    output$boxPlot <- renderPlot({
      boxAndWhiskerPlot(input$geoAccessionCode, input$platform, knnDataInput())
      })
    
    # Expression Value Distribution Plot
    output$expressionDensity <- renderPlot({
      expressionValueDistributionPlot(input$geoAccessionCode, input$platform, knnDataInput())
    })
    
    # Mean-Variance Plot
    output$meanVariance <- renderPlot({
      meanVariancePlot(input$geoAccessionCode, input$platform, knnDataInput())
    })
    
    # UMAP plot (multi-dimensional scaling)
    output$umap <- renderPlot({
      umapPlot(input$geoAccessionCode, input$platform, knnDataInput())
      })
    
    # These components need to be moved into a function
    # Principal component analysis scree plot
    output$pcaScreePlot <- renderPlot({
      fviz_eig(pcaDataInput())
    })
    
    # Principal component analysis individuals plot
    output$pcaIndividualsPlot <- renderPlot({
      fviz_pca_ind(pcaDataInput(),
                   col.ind = "cos2", # Color by the quality of representation
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   geom = "point",
                   repel = TRUE     # Avoid text overlapping
      )
    })
    
    # Principal component analysis variables plot
    output$pcaVariablesPlot <- renderPlot({
      fviz_pca_var(pcaDataInput(),
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   repel = TRUE     # Avoid text overlapping
      )
    })
    
    # Principal component analysis biplot of individuals and variables
    output$pcaBiplotPlot <- renderPlot({
      fviz_pca_biplot(pcaDataInput(), repel = TRUE,
                      col.var = "#2E9FDF", # Variables color
                      geom = "point",
                      col.ind = "#696969",  # Individuals color
      )
    })
    
  }
  
  shinyApp(ui, server)