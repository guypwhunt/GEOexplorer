# Load Packages
library(shiny)
library(plotly)
library(ggplot2)
library(shinyBS)
library(shinyjs)

source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
source("differentialGeneExpressionAnalysis/differentialGeneExpressionAnalysis.R")
source("interactiveDataVisualizationFunctions/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R")


ui <- fluidPage(
  titlePanel("GEO2R Data Visualisation"),
  helpText("GEO2R is an interactive web tool that allows users to compare two or more groups of samples in a GEO Series to identify genes that are differentially expressed across experimental conditions. GEO2R Data Visualisation extends GEO2R's functionalities by enabling a richer set of analysis and graphics to be performed/generated from the GEO2R gene expression data."),
  sidebarPanel(
    helpText("Input a GEO accession code to examine the gene expression data."),
    textInput("geoAccessionCode", "GEO accession code", "GSE18388"),
    helpText("Select the platform of interest."),
    selectInput("platform", "Platform",c()),
    radioButtons("logTransformation",
                 label="Apply log transformation to the data:",
                 choices=list("Auto-Detect","Yes","No"),
                 selected="Auto-Detect"),
    bsTooltip(id = "logTransformation", title = "The GEO database accepts a variety of data value types, including logged and unlogged data. Limma expects data values to be in log space. To address this, an auto-detect feature that checks the values of selected samples and automatically performs a log2 transformation on values determined not to be in log space. Alternatively, the user can select Yes to force log2 transformation, or No to override the auto-detect feature. The auto-detect feature only considers Sample values that have been assigned to a group, and applies the transformation in an all-or-none fashion", placement = "top", trigger = "hover"),
    uiOutput("logTransformationText"),
    br(),
    radioButtons("knnTransformation",
                 label="Apply k-nearest neighbors (KNN) algorithm to predict null data:",
                 choices=list("Yes","No"),
                 selected="No"),
    bsTooltip(id = "knnTransformation", title = "Rows with over 50% missing values are imputed using the overall mean per sample. Columns with over 80% will cause an error in the KNN computation.", placement = "top", trigger = "hover"),
    br(),
    radioButtons("platformAnnotation",
                 label="Category of Platform annotation to display on results:",
                 choices=list("Submitter supplied","NCBI generated"),
                 selected="NCBI generated")
  ),
  mainPanel(tabsetPanel(type = "tabs",
                        tabPanel("Dataset Information",
                                 tabsetPanel(type = "tabs",
                                            tabPanel("Experiment Information", br(), htmlOutput('experimentInfo')),
                                            tabPanel("Column Details", dataTableOutput('columnTable')),
                                            tabPanel("Dataset", dataTableOutput('table'))
                                            )),
                        tabPanel("Exploratory Data Analysis",
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Box-and-Whisper Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. The quartiles are calculated using the linear method. Viewing the distribution can be useful for determining if the data in the dataset is suitable for differential expression analysis. Generally, median-centred values are indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveBoxAndWhiskerPlot')),
                                             tabPanel("Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveDesnityPlot')),
                                             tabPanel("3D Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveThreeDDesnityPlot')),
                                             tabPanel("Mean-Variance Plot", br(), span("Generated using R limma and plotly. The plot below is used to check the mean-variance relationship of the expression data, after fitting a linear model. It can help show if there is a lot of variation in the data. Each point represents a gene. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveMeanVariancePlot')),
                                             tabPanel("UMAP Plot", br(), span("Generated using R umap and plotly. Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique useful for visualizing how genes are related to each other. The number of nearest neighbours used in the calculation is indicated in the graph. The plot shows data after log and KNN transformation if they were performed."), br(), br(), numericInput("knn", "Input the k-nearest neighbors value  to use:", 2, min = 2,step = 1), br(), plotlyOutput('interactiveUmapPlot')),
                                             tabPanel("PCA Analysis",
                                                      tabsetPanel(type = "tabs",
                                                                  tabPanel("Scree Plot", br(), span("Generated using R princomp and plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (Comp). The plot displays the eigenvalues against the number of dimensions. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaScreePlot')),
                                                                  tabPanel("Individuals Plot", br(), span("Generated using R princomp and R plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (Comp). The plot displays the eigenvalues for each individual (row) in the gene expression dataset for the top two principal components (Comp.1 and Comp.2). The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaIndividualsPlot')),
                                                                  tabPanel("Variables Plot", br(), span("Generated using R princomp and R plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (Comp). The plot displays the eigenvalues for each variable (column) in the gene expression dataset for the top two principal components (Comp.1 and Comp.2). The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaVariablesPlot'))                                                      )
                                             )
                                             
                                 )
                        ),
                        tabPanel("Differential Gene Expression Analysis",
                                 tabsetPanel(type = "tabs",
                                             tabPanel("Set Parameters", 
                                                      mainPanel(
                                                        fluidRow(
                                                        column(6,
                                                               br(), 
                                                               selectInput("columns1", "Group 1 Columns", choices=c(), multiple = TRUE),
                                                               selectInput("columns2", "Group 2 Columns", choices=c(), multiple = TRUE),
                                                               selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
                                                                           choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Holm", "None")) # "Hochberg" and "Hommel" were removed
                                                               ),
                                                        br(),
                                                        column(6,
                                                               radioButtons("limmaPrecisionWeights",
                                                                            label="Apply limma precision weights (vooma):",
                                                                            choices=list("Yes","No"),
                                                                            selected="No"),
                                                               radioButtons("forceNormalization",
                                                                            label="Force normalization:",
                                                                            choices=list("Yes","No"),
                                                                            selected="No"), 
                                                               sliderInput("significanceLevelCutOff", "Significance level cut-off:",
                                                                           min = 0, max = 1,
                                                                           value = 0.05),
                                                               actionButton("differentialExpressionButton", "Analyse")
                                                               )
                                                        )
                                                        )
                                                      ),
                                             tabPanel("Top Differentially Expressed Genes", br(), span("The table below displays the top differentially expressed genes between the groups selected."), br(), span("adj.P.Val is the P-value after adjustment for multiple testing. This column is generally recommended as the primary statistic by which to interpret results. Genes with the smallest P-values will be the most reliable."), br(), span("P.Value is the Raw P-value"), br(), span("t is the Moderated t-statistic (only available when two groups of Samples are defined)"), br(), span("B is the B-statistic or log-odds that the gene is differentially expressed (only available when two groups of Samples are defined)"), br(), span("logFC is the Log2-fold change between two experimental conditions (only available when two groups of Samples are defined)"), br(), span("F is the moderated F-statistic which combines the t-statistics for all the pair-wise comparisons into an overall test of significance for that gene (only available when more than two groups of Samples are defined)"), br(), dataTableOutput('dETable')),
                                             tabPanel("Histogram Plot", br(), span("Generated using hist. Use to view the distribution of the P-values in the analysis results. The P-value here is the same as in the Top differentially expressed genes table and computed using all selected contrasts. While the displayed table is limited by size this plot allows you to see the 'big picture' by showing the P-value distribution for all analyzed genes.") , plotOutput('dEHistogram')),
                                             tabPanel("Interactive Histogram Plot", br(), span("") , plotlyOutput('iDEHistogram')),
                                             tabPanel("Venn Diagram Plot", br(), span("Generated using limma (vennDiagram). Use to explore the overlap in significant genes between multiple contrasts.") , plotOutput('dEVennDiagram')),
                                             tabPanel("Q-Q Plot", br(), span("Generated using limma (qqt). Plots the quantiles of a data sample against the theoretical quantiles of a Student's t distribution. This plot helps to assess the quality of the limma test results. Ideally the points should lie along a straight line, meaning that the values for moderated t-statistic computed during the test follow their theoretically predicted distribution."), plotOutput('dEQQ')),
                                             tabPanel("Volcano Plot", br(), span("Generated using limma (volcanoplot). A volcano plot displays statistical significance (-log10 P value) versus magnitude of change (log2 fold change) and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff value") , plotOutput('dEVolcano')),
                                             tabPanel("Interactive Volcano Plot", br(), span("") , plotlyOutput('iDEVolcano')),
                                             tabPanel("Mean Difference Plot", br(), span("Generated using limma (plotMD). A mean difference (MD) plot displays log2 fold change versus average log2 expression values and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff (red = upregulated, blue = downregulated)."), plotOutput('dEMd')),
                                             tabPanel("Interactive Mean Difference Plot", br(), span(""), plotlyOutput('iDEMd'))
                                             
                                             )
                                 
                        )
  )
  )
)

server <- function(input, output, session){
  # Data Extraction Functions
  # Get the GEO2R data for all platforms
  allGset <- reactive({getGset(input$geoAccessionCode, input$platformAnnotation)})
  
  # Get a list of all the platforms
  platforms <- reactive({getPlatforms(allGset())})
  
  # Extract the GEO2R data from the specified platform
  gsetData <- reactive({getPlatformGset(allGset(), input$platform)})
  
  # Extract the experiment information 
  experimentInformation <- reactive({getExperimentInformation(gsetData())})
  
  # Extract expression data
  expressionData <- reactive({extractExpressionData(gsetData())})
  
  # Extract Column Information
  columnInfo <- reactive({getColumnDetails(gsetData())})
  
  # Is log transformation auto applied
  autoLogInformation <- reactive({isLogTransformAutoApplied(expressionData())})
  
  # Get a list of all the columns
  columns <- reactive({extractColumns(expressionData())})
  
  
  # Data Transformation Functions
  # Apply log transformation to expression data if necessary
  dataInput <- reactive({logTransformExpressionData(expressionData(), input$logTransformation)
  })
  
  # Perform KNN transformation on log expression data if necessary
  knnDataInput <- reactive({knnDataTransformation(dataInput(), input$knnTransformation)
  })
  
  # Remove all incomplete rows
  naOmitInput <-reactive({naOmitTransformation(knnDataInput())
  })
  
  # Perform PCA analysis on KNN transformation expression data using princomp
  pcaPrincompDataInput <- reactive({pcaPrincompAnalysis(naOmitInput())
  })
  
  
  # Data Visualisation Functions
  # Update Platform Options
  platformObserve <- observe({
    updateSelectInput(session, "platform",
                      choices = platforms(),
                      selected = tail(platforms(), 1))
  })
  
  # Update if log transformation took place
  output$logTransformationText <- renderUI({
    helpText(autoLogInformation())
  })
  
  # Experimental Information Display
  output$experimentInfo <- renderUI({
    extractExperimentInformation(experimentInformation())
  })
  
  # Column Set Plot
  output$columnTable <- renderDataTable({
    columnInfo()
    })
  
  # Data Set Plot
  output$table <- renderDataTable({
    knnDataInput()
  })
  
  # Interactive Box-and-Whisker Plot
  output$interactiveBoxAndWhiskerPlot <- renderPlotly({
    interactiveBoxAndWhiskerPlot(naOmitInput(), input$geoAccessionCode, input$platform)
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
    interactiveUmapPlot(naOmitInput(), input$knn, input$geoAccessionCode)
  })
  
  # Interactive Mean Variance Plot
  output$interactiveMeanVariancePlot <- renderPlotly({
    interactiveMeanVariancePlot(naOmitInput(),input$geoAccessionCode)
  })
  
  # Interactive PCA Scree Plot
  output$interactivePcaScreePlot <- renderPlotly({
    interactivePrincompPcaScreePlot(pcaPrincompDataInput(), input$geoAccessionCode)
  })
  
  # Interactive PCA Individual Plot
  output$interactivePcaIndividualsPlot <- renderPlotly({
    interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput(), input$geoAccessionCode)
  })
  
  # Interactive PCA Variables Plot
  output$interactivePcaVariablesPlot <- renderPlotly({
    interactivePrincompPcaVariablesPlot(pcaPrincompDataInput(), input$geoAccessionCode)
  })
  
  
  # Differential Gene Expression Functions
  # Update Column on UI
  columns1Observe <- observe({
    updateSelectInput(session, "columns1",
                      choices = columns())
  })
  
  columns2Observe <- observe({
    updateSelectInput(session, "columns2",
                      choices = columns())
  })

  observeEvent(input$differentialExpressionButton, {
    #dFAllGset <- getGset(input$geoAccessionCode, input$platformAnnotation, GSEMatrix=TRUE, getGPL=TRUE)
    #dFGsetData <- getPlatformGset(dFAllGset, input$platform)
    gsms <- calculateGsms(columns(),input$columns1, input$columns2)
    fit2 <- calculateFit2(input$geoAccessionCode, input$platform, input$platformAnnotation, gsms, input$logTransformation, input$limmaPrecisionWeights, input$forceNormalization, input$knnTransformation)
    #fit2 <- differentialGeneExpression(gsetData(), knnDataInput(), gsms, input$limmaPrecisionWeights, input$forceNormalization)
    adjustment <- adjustmentCalculation(input$pValueAdjustment)
    tT <- topDifferentiallyExpressedGenesTable(fit2, adjustment)
    dT <- dT(fit2, adjustment, input$significanceLevelCutOff)
    ct <- 1  
    
    output$dETable <- renderDataTable({
      as.data.frame(tT)
    })
    
    
    output$dEHistogram <- renderPlot({
      histogramPlot(fit2, adjustment)
    })
    
    output$iDEHistogram <- renderPlotly({
      interactiveHistogramPlot(fit2, adjustment)
    })
    
    output$dEVennDiagram <- renderPlot({
      vennDiagramPlot(dT)
    })
    
    output$dEQQ <- renderPlot({
      qqPlot(fit2)
    })
    
    output$dEVolcano <- renderPlot({
      volcanoPlot(fit2, dT, ct)
    })
    
    output$iDEVolcano <- renderPlotly({
      interactiveVolcanoPlot(fit2, adjustment, dT, ct)
    })
    
    output$dEMd <- renderPlot({
      mdPlot(fit2, dT, ct)
    })
    
    output$iDEMd <- renderPlotly({
      interactiveMeanDifferencePlot(fit2, adjustment, dT, ct)
    })
    
    updateActionButton(session, "differentialExpressionButton",
                       label = "Successful")
  })
  
}

shinyApp(ui, server)