  # Load Packages
  library(shiny)
  library(plotly)
  library(ggplot2)
  

  source("geoIntegrationFunctions/geoIntegrationFunctions.R")
  source("dataVisualizationFunctions/dataVisualizationFunctions.R")
  source("dataTransformationFunctions/dataTransformationFunctions.R")
  source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
  
  ui <- fluidPage(
    titlePanel("GEO2R Data Visualisation"),
    helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of Samples in a GEO Series to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of analysis and graphics to be performed/generated from the GEO2R gene expression data."),
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
      br()#,
      #helpText("The Buttons Below here currently are not functional!!!"),
      # Need to get the radio button below working or delete
      #selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
      #            choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Hochberg", "Holm", "Hommel", "None")),
      #radioButtons("limmaPrecisionWeights",
      #             label="Apply limma precision weights (vooma):",
      #             choices=list("Yes","No"),
      #             selected="No"), 
      # Need to get the radio button below working or delete
      #radioButtons("forceNormalization",
      #             label="Force normalization:",
      #             choices=list("Yes","No"),
      #             selected="No"), 
      # Need to get the radio button below working or delete
      #radioButtons("forceNormalization",
      #             label="Category of Platform annotation to display on results:",
      #             choices=list("Submitter supplied","NCBI generated"),
      #             selected="No"), 
      #br(),
      #helpText("Plot displays"),
      #sliderInput("integer", "Significance level cut-off:",
      #            min = 0, max = 1,
      #            value = 0.05),
      # add grouping functionality
    ),
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Dataset", dataTableOutput('myTable')),
                          tabPanel("Exploratory Data Analysis",
                                   tabsetPanel(type = "tabs",
                                   tabPanel("Original Visualizations",
                                            tabsetPanel(type = "tabs",
                                                        tabPanel("Box-and-Whisper Plot", br(), span("Generated using R boxplot. The plot below displays the distribution of the values of the genes in the dataset. Viewing the distribution can be useful for determining if the data in the dataset is suitable for differential expression analysis. Generally, median-centred values are indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('boxPlot')),
                                                        tabPanel("Expression Density Plot", br(), span("Generated using R limma (plotDensities). The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('expressionDensityPlot')),
                                                        tabPanel("Mean-Variance Plot", br(), span("Generated using R limma (plotSA, vooma). The plot below is used to check the mean-variance relationship of the expression data, after fitting a linear model. It can help show if there is a lot of variation in the data. Each point represents a gene. The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('meanVariancePlot')),
                                                        tabPanel("UMAP Plot", br(), span("Generated using R umap. Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique useful for visualizing how genes are related to each other. The number of nearest neighbors used in the calculation is indicated in the graph. The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('umapPlot')))),
                                   tabPanel("Interactive Visualizations",
                                            tabsetPanel(type = "tabs",
                                                        tabPanel("Box-and-Whisper Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. The quartiles are calculated using the linear method. Viewing the distribution can be useful for determining if the data in the dataset is suitable for differential expression analysis. Generally, median-centred values are indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveBoxAndWhiskerPlot')),
                                                        tabPanel("Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveDesnityPlot')),
                                                        tabPanel("3D Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveThreeDDesnityPlot')),
                                                        tabPanel("Mean-Variance Plot", br(), span("Generated using R limma and plotly. The plot below is used to check the mean-variance relationship of the expression data, after fitting a linear model. It can help show if there is a lot of variation in the data. Each point represents a gene. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveMeanVariancePlot')),
                                                        tabPanel("UMAP Plot", br(), span("Generated using R umap and plotly. Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique useful for visualizing how genes are related to each other. The number of nearest neighbors used in the calculation is indicated in the graph. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveUmapPlot'))
                                                        
                                                                                                    ))
                                                        
                                   )),
                          tabPanel("Principal Component Analysis",
                                   tabsetPanel(type = "tabs",
                                               tabPanel("Scree Plot", br(), span("Generated using R prcomp and visualised using R fviz_eig. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information. "), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues against the number of dimensions. The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('pcaScreePlot')),
                                               tabPanel("Individuals Plot", br(), span("Generated using R prcomp and visualised using R fviz_pca. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues for each individual (row) in the gene expression dataset for the top two principal components (Dim 1 and Dim 2). The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('pcaIndividualsPlot')),
                                               tabPanel("Variables Plot", br(), span("Generated using R prcomp and visualised using R fviz_pca. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues for each variable (column) in the gene expression dataset for the top two principal components (Dim 1 and Dim 2). The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('pcaVariablesPlot')),
                                               tabPanel("Individuals and Variables Biplot",  br(), span("Generated using R prcomp and visualised using R fviz_pca. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues for each variable (column) and individual (row) in the gene expression dataset for the top two principal components (Dim 1 and Dim 2). The plot shows data after log and KNN transformation if they were performed."), br(), plotOutput('pcaBiplotPlot'))
                                               )
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
    
    # Remove all incomplete rows
    naOmitInput <-reactive({naOmitTransformation(knnDataInput())
    })
    
    # Perform PCA analysis on KNN transformation expression data
    pcaDataInput <- reactive({pcaAnalysis(naOmitInput())
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
    output$expressionDensityPlot <- renderPlot({
      expressionValueDistributionPlot(input$geoAccessionCode, input$platform, knnDataInput())
    })
    
    # Mean-Variance Plot
    output$meanVariancePlot <- renderPlot({
      meanVariancePlot(input$geoAccessionCode, input$platform, naOmitInput())
    })
    
    # This was updated to mandatoryily use the KNN data, this may need to be reverted
    # UMAP plot (multi-dimensional scaling)
    output$umapPlot <- renderPlot({
      umapPlot(input$geoAccessionCode, input$platform, naOmitInput())
      })
    
    # Principal component analysis scree plot
    output$pcaScreePlot <- renderPlot({
      pcaScreePlot(pcaDataInput())
    })
    
    # Principal component analysis individuals plot
    output$pcaIndividualsPlot <- renderPlot({
      pcaIndividualsPlot(pcaDataInput())
    })
    
    # Principal component analysis variables plot
    output$pcaVariablesPlot <- renderPlot({
      pcaVariablesPlot(pcaDataInput())
    })
    
    # Principal component analysis biplot of individuals and variables
    output$pcaBiplotPlot <- renderPlot({
      pcaBiplotPlot(pcaDataInput())
    })
    
    # Interactive Box-and-Whisker Plot
    output$interactiveBoxAndWhiskerPlot <- renderPlotly({
      interactiveBoxAndWhiskerPlot(knnDataInput(), input$geoAccessionCode, input$platform)
    })
    
    # Interactive Density Plot
    output$interactiveDesnityPlot <- renderPlotly({
      interactiveDesnityPlot(naOmitInput(), input$geoAccessionCode, input$platform)
    })
    
    # 3D Interactive Density Plot
    output$interactiveThreeDDesnityPlot <- renderPlotly({
      interactiveThreeDDesnityPlot(naOmitInput(), input$geoAccessionCode, input$platform)
    })
    
    # Interactive UMAP Plot
    output$interactiveUmapPlot <- renderPlotly({
      interactiveUmapPlot(naOmitInput())
    })
    
    # Interactive Mean Variance Plot
    output$interactiveMeanVariancePlot <- renderPlotly({
      interactiveMeanVariancePlot(naOmitInput(),input$geoAccessionCode)
    })
  }
  
  shinyApp(ui, server)