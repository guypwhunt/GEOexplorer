#' A Function to Return the Side Bar Ui Component
#'
#' A Function to Return the Side Bar Ui Component
#' @import shiny
#' @examples sourceSideBarUi()
#' @importFrom shinyBS bsTooltip
#' @author Guy Hunt
sourceSideBarUi <- function() {
  sideBarUi <- sidebarPanel(
    helpText("Input a GEO series accession code (GSEXXXX format) to examine the gene expression data."),
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
    actionButton("exploratoryDataAnalysisButton", "Analyse")
  )
  return(sideBarUi)
}

#' A Function to Return the Dataset Information Ui Component
#'
#' A Function to Return the Dataset Information Ui Component
#' @import shiny
#' @examples sourceDatasetInformationUi()
#' @author Guy Hunt
sourceDatasetInformationUi <- function() {
  datasetInformationUi <- tabPanel("Dataset Information",
                                   tabsetPanel(type = "tabs",
                                               tabPanel("Experiment Information", br(), htmlOutput('experimentInfo')),
                                               tabPanel("Column Details", DT::dataTableOutput('columnTable')),
                                               tabPanel("Dataset", DT::dataTableOutput('table'))
                                   ))
  return(datasetInformationUi)
  }

#' A Function to Return the Differential Gene Expression Analysis Ui Component
#'
#' A Function to Return the Differential Gene Expression Analysis Ui Component
#' @import shiny plotly
#' @examples sourceDifferentialGeneExpressionAnalysisUi()
#' @importFrom DT dataTableOutput
#' @author Guy Hunt
sourceDifferentialGeneExpressionAnalysisUi <- function() {
  differentialGeneExpressionAnalysisUi <- tabPanel("Differential Gene Expression Analysis",
                                                   tabsetPanel(type = "tabs",
                                                               tabPanel("Set Parameters",
                                                                        mainPanel(
                                                                          DT::dataTableOutput('knnColumnTable'),
                                                                          fluidRow(
                                                                            column(6,
                                                                                   br(),
                                                                                   uiOutput("dyncolumns"),
                                                                                   selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
                                                                                               choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Holm", "None")), # "Hochberg" and "Hommel" were removed
                                                                                   radioButtons("limmaPrecisionWeights",
                                                                                                label="Apply limma precision weights (vooma):",
                                                                                                choices=list("Yes","No"),
                                                                                                selected="No")
                                                                            ),
                                                                            br(),
                                                                            column(6,
                                                                                   radioButtons("forceNormalization",
                                                                                                label="Force normalization:",
                                                                                                choices=list("Yes","No"),
                                                                                                selected="No"),
                                                                                   sliderInput("significanceLevelCutOff", "Significance level cut-off:",
                                                                                               min = 0, max = 1,
                                                                                               value = 0.05),
                                                                                   uiOutput("differentialExpressionButton")                                                                            )
                                                                          )
                                                                        )
                                                               ),
                                                               tabPanel("Top Differentially Expressed Genes", br(), span("The table below displays the top differentially expressed genes between the groups selected."), br(), br(), span("adj.P.Val is the P-value after adjustment for multiple testing. This column is generally recommended as the primary statistic by which to interpret results. Genes with the smallest P-values will be the most reliable."), br(), span("P.Value is the Raw P-value"), br(), span("t is the Moderated t-statistic (only available when two groups of Samples are defined)"), br(), span("B is the B-statistic or log-odds that the gene is differentially expressed (only available when two groups of Samples are defined)"), br(), span("logFC is the Log2-fold change between two experimental conditions (only available when two groups of Samples are defined)"), br(), span("F is the moderated F-statistic which combines the t-statistics for all the pair-wise comparisons into an overall test of significance for that gene (only available when more than two groups of Samples are defined)"),  br(), br(), downloadButton("downloadData", "Download"), br(), br(), dataTableOutput('dETable')),
                                                               tabPanel("Histogram Plot", br(), span("Generated using hist. Use to view the distribution of the P-values in the analysis results. The P-value here is the same as in the Top differentially expressed genes table and computed using all selected contrasts. While the displayed table is limited by size this plot allows you to see the 'big picture' by showing the P-value distribution for all analyzed genes."), br(), br(), plotlyOutput('iDEHistogram')),
                                                               tabPanel("Venn Diagram Plot", br(), span("Generated using limma (vennDiagram). Use to explore the overlap in significant genes between multiple contrasts.") , plotOutput('dEVennDiagram')),
                                                               tabPanel("Q-Q Plot", br(), span("Generated using limma (qqt) and R plotly. Plots the quantiles of a data sample against the theoretical quantiles of a Student's t distribution. This plot helps to assess the quality of the limma test results. Ideally the points should lie along a straight line, meaning that the values for moderated t-statistic computed during the test follow their theoretically predicted distribution."), br(), br(), plotlyOutput('iDEQQ')),
                                                               tabPanel("Volcano Plot", br(), span("Generated using R plotly. A volcano plot displays statistical significance (-log10 P value) versus magnitude of change (log2 fold change) and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff value") , br(), br(), plotlyOutput('iDEVolcano')),
                                                               tabPanel("Mean Difference Plot", br(), span("Generated using R plotly. A mean difference (MD) plot displays log2 fold change versus average log2 expression values and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff."), br(), br(), plotlyOutput('iDEMd'))

                                                   )

  )
  return(differentialGeneExpressionAnalysisUi)
}

#' A Function to Return the Exploratory Data Analysis Ui Component
#'
#' A Function to Return the Exploratory Data Analysis Ui Component
#' @import shiny plotly
#' @examples sourceExploratoryDataAnalysisUi()
#' @author Guy Hunt
sourceExploratoryDataAnalysisUi <- function() {
  exploratoryDataAnalysisUi <- tabPanel("Exploratory Data Analysis",
                                        tabsetPanel(type = "tabs",
                                                    tabPanel("Box-and-Whisper Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. The quartiles are calculated using the linear method. Viewing the distribution can be useful for determining if the data in the dataset is suitable for differential expression analysis. Generally, median-centred values are indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveBoxAndWhiskerPlot')),
                                                    tabPanel("Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveDensityPlot')),
                                                    tabPanel("3D Expression Density Plot", br(), span("Generated using R plotly. The plot below displays the distribution of the values of the genes in the dataset. This plot complements the boxplot in checking for data normalization before differential expression analysis. If density curves are similar from gene to gene, it is indicative that the data is normalized and cross-comparable. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveThreeDDensityPlot')),
                                                    tabPanel("Mean-Variance Plot", br(), span("Generated using R limma and plotly. The plot below is used to check the mean-variance relationship of the expression data, after fitting a linear model. It can help show if there is a lot of variation in the data. Each point represents a gene. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactiveMeanVariancePlot')),
                                                    tabPanel("UMAP Plot", br(), span("Generated using R umap and plotly. Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique useful for visualizing how genes are related to each other. The number of nearest neighbours used in the calculation is indicated in the graph. The plot shows data after log and KNN transformation if they were performed."), br(), br(), numericInput("knn", "Input the k-nearest neighbors value  to use:", 2, min = 2,step = 1), br(), plotlyOutput('interactiveUmapPlot')),
                                                    tabPanel("Heatmap Plot", br(), span("Generated using R cor and heatmaply. The plot below compares the correlation values of the samples in a heatmap."), br(), br(), plotlyOutput('interactiveHeatMapPlot')),
                                                    tabPanel("PCA Analysis",
                                                             tabsetPanel(type = "tabs",
                                                                         tabPanel("Scree Plot", br(), span("Generated using R prcomp and plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues against the number of dimensions. The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaScreePlot')),
                                                                         tabPanel("Individuals Plot", br(), span("Generated using R prcomp and R plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues for each individual (row) in the gene expression dataset for the top two principal components (PC1 and PC2). The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaIndividualsPlot')),
                                                                         tabPanel("Variables Plot", br(), span("Generated using R prcomp and R plotly. Principal component analysis (PCA) reduces the dimensionality of multivariate data to two dimensions that can be visualized graphically with minimal loss of information."), br(), span("Eigenvalues correspond to the amount of the variation explained by each principal component (PC). The plot displays the eigenvalues for each variable (column) in the gene expression dataset for the top two principal components (PC1 and PC2). The plot shows data after log and KNN transformation if they were performed."), br(), br(), plotlyOutput('interactivePcaVariablesPlot'))                                                      )
                                                    )

                                        )
  )
  return(exploratoryDataAnalysisUi)
}
