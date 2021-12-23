#' A Function to Return the Ui Components
#'
#' A Function to Return the Ui Components
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples navbarPage()
#' @author Guy Hunt
#' @import markdown
#' @noRd
sourceUi <- function() {
  uiComponents <- navbarPage(
    "GEOexplorer",
    id = "geoexplorerNavBar",
    tabPanel("Home", value = "Home",
             # Source the Side Bar UI Components
             sourceSideBarUi(),
             mainPanel(
               tabsetPanel(
                 type = "tabs",
                 # Source the Dataset Information
                 # UI Components
                 sourceDatasetInformationUi(),
                 # Source the Exploratory Data Analysis
                 # UI Components
                 sourceExploratoryDataAnalysisUi(),
                 # Source the Differential Gene Expression
                 # UI Components
                 sourceDifferentialGeneExpressionAnalysisUi(),
                 # Source the Enrichment
                 # UI Components
                 sourceEnrichmentnUi()
               )
             )),
    tabPanel("About",
             htmlOutput("about")
             ),
    tabPanel("GEO Search",
             sourceGeoSearchUi()),
    tabPanel("Tutorial",
             htmlOutput('tutorial')),
    tabPanel("Example Datasets",
             # Source example datasets
             sourceExampleUI())
  )
  return(uiComponents)
}

#' A Function to Return the Side Bar Ui Component
#'
#' A Function to Return the Side Bar Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceSideBarUi()
#' @author Guy Hunt
#' @noRd
sourceSideBarUi <- function() {
  sideBarUi <- sidebarPanel(
    radioButtons(
      "dataSetType",
      label = "Would you like analyse a single gene exression dataset or
      combine two gene exression datasets?",
      choices = list("Single", "Combine"),
      selected = "Single"
    ),
    uiOutput("output2"),
    radioButtons(
      "dataSource",
      label = "Would you like to upload the gene expression data
      or source the data from GEO?",
      choices = list("GEO", "Upload"),
      selected = "GEO"
    ),
    uiOutput("output4"),
    uiOutput("output5"),
    uiOutput("output6"),
    uiOutput("output7"),
    uiOutput("output8"),
    uiOutput("output9"),
    uiOutput("output10"),
    uiOutput("output11"),
    radioButtons(
      "logTransformation",
      label = "Apply log transformation to the data:",
      choices = list("Auto-Detect", "Yes", "No"),
      selected = "Auto-Detect"
    ),
    uiOutput("logTransformationText"),
    uiOutput("output13"),
    uiOutput("output14"),
    actionButton("exploratoryDataAnalysisButton", "Analyse")
  )
  return(sideBarUi)
}

#' A Function to Return the Dataset Information Ui Component
#'
#' A Function to Return the Dataset Information Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT dataTableOutput
#' @examples sourceDatasetInformationUi()
#' @author Guy Hunt
#' @noRd
sourceDatasetInformationUi <- function() {
  datasetInformationUi <- tabPanel("Dataset Information",
                                   tabsetPanel(
                                     type = "tabs",
                                     tabPanel(
                                       "Experiment Information",
                                       br(),
                                       span(
                                         "Summary information of the gene
                                         expression study is displayed below."
                                       ),
                                       br(),
                                       br(),
                                       htmlOutput('experimentInfo')
                                     ),
                                     tabPanel(
                                       "Experimental Conditions Information",
                                       br(),
                                       span(
                                         "A table containing information for
                                         each of the experimental conditions
                                         used in the gene expression study is
                                         displayed below.
                                         Each experimental condition relates to
                                         a column in the gene expression
                                         dataset in the 'Gene Expression
                                         Dataset' tab."
                                       ),
                                       br(),
                                       br(),
                                       dataTableOutput('columnTable')
                                     ),
                                     tabPanel(
                                       "Gene Expression Dataset",
                                       br(),
                                       span(
                                         "A table containing the gene
                                         expression data is displayed below.
                                         Each column
                                         relates to an experimental condition,
                                         each row relates to a gene, and each
                                         value relates to a gene expression
                                         value for that gene under that
                                         experimental condition. The values
                                         are displayed post KNN imputation,
                                         count per million transformation and
                                         log transformation if selected."
                                       ),
                                       br(),
                                       br(),
                                       downloadButton("downloadGeneExpression",
                                                      "Download"),
                                       br(),
                                       br(),
                                       dataTableOutput('table')
                                     )
                                   ))
  return(datasetInformationUi)
}

#' A Function to Return the Differential Gene Expression Analysis Ui Component
#'
#' A Function to Return the Differential Gene Expression Analysis Ui Component
#' @import plotly
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceDifferentialGeneExpressionAnalysisUi()
#' @importFrom DT dataTableOutput
#' @author Guy Hunt
#' @noRd
sourceDifferentialGeneExpressionAnalysisUi <- function() {
  differentialGeneExpressionAnalysisUi <-
    tabPanel(
      "Differential Gene Expression Analysis",
      tabsetPanel(
        type = "tabs",
        tabPanel(
          "Set Parameters",
          br(),
          span(
            "A table containing information for each of the experimental
              conditions used in the gene expression study is displayed below.
              In the group column, select the experimental conditions you want
              to include in group 1, group 2 or N/A if you want the
              experimental condition excluded from differential gene
              expression analysis. During differential gene expression
              analysis, group 1 is compared against group 2."
          ),
          br(),
          br(),
          tags$b("Select the experimental conditions to include in Group 1."),
          br(),
          br(),
          dataTableOutput('knnColumnTableOne'),
          br(),
          br(),
          tags$b("Select the experimental conditions to include in Group 2."),
          br(),
          br(),
          dataTableOutput('knnColumnTableTwo'),
          br(),
          br(),
          span(
            "The parameters for differential gene expression analysis are
              displayed below. Please select the appropriate parameters and
              click analyse to perform differential gene expression analysis."
          ),
          br(),
          br(),
          fluidRow(
            column(
              6,
              br(),
              uiOutput("dyncolumns"),
              selectInput(
                "pValueAdjustment",
                "Apply adjustment to the P-values:",
                choices = c(
                  "Benjamini & Hochberg (False discovery rate)",
                  "Benjamini & Yekutieli",
                  "Bonferroni",
                  "Holm",
                  "None"
                )
              ),
              # "Hochberg" and "Hommel" were removed
              radioButtons(
                "limmaPrecisionWeights",
                label =
                  "Apply limma precision weights (vooma):",
                choices =
                  list("Yes", "No"),
                selected =
                  "No"
              )
            ),
            br(),
            column(
              6,
              radioButtons(
                "forceNormalization",
                label =
                  "Force normalization:",
                choices =
                  list("Yes", "No"),
                selected =
                  "No"
              ),
              sliderInput(
                "significanceLevelCutOff",
                "Significance level cut-off:",
                min = 0,
                max = 1,
                value = 0.05
              ),
              uiOutput("output100"),
              br(),
              br()
            )
          )
        ),
        tabPanel(
          "Top Differentially Expressed Genes",
          br(),
          span(
            "The table below displays the top differentially expressed genes
            between the groups selected."
          ),
          br(),
          br(),
          span(
            "adj.P.Val is the P-value after adjustment for multiple testing.
            This column is generally recommended as the primary statistic by
            which to interpret results.
            Genes with the smallest P-values will be the most reliable."
          ),
          br(),
          br(),
          span("P.Value is the Raw P-value"),
          br(),
          br(),
          span("t is the Moderated t-statistic"),
          br(),
          br(),
          span(
            "B is the B-statistic or log-odds that the gene is
            differentially expressed"
          ),
          br(),
          br(),
          span(
            "logFC is the Log2-fold change between two experimental
            conditions"
          ),
          br(),
          br(),
          span(
            "F is the moderated F-statistic which combines
            the t-statistics for all the pair-wise comparisons into an
            overall test of significance for that gene"
          ),
          br(),
          br(),
          downloadButton("downloadData", "Download"),
          br(),
          br(),
          dataTableOutput('dETable')
        ),
        tabPanel(
          "Histogram Plot",
          br(),
          span(
            "Generated using hist. Use to view the distribution of
            the P-values in the analysis results.
            The P-value here is the same as in the Top differentially
            expressed genes table and computed
            using all selected contrasts. While the displayed table is
            limited by size this plot allows
            you to see the 'big picture' by showing the P-value
            distribution for all analyzed genes."
          ),
          br(),
          br(),
          plotlyOutput('iDEHistogram')
        ),
        tabPanel(
          "Venn Diagram Plot",
          br(),
          span(
            "Generated using limma (vennDiagram).
            Displays the number of differentially expressed genes versus the
            number of non-differentially expressed genes."
          ),
          plotOutput('dEVennDiagram')
        ),
        tabPanel(
          "Q-Q Plot",
          br(),
          span(
            "Generated using limma (qqt) and R plotly.
            Plots the quantiles of a data sample against the
            theoretical quantiles of a Student's t distribution.
            This plot helps to assess the quality of the limma
            test results. Ideally the points should lie along a straight line,
            meaning that the values for moderated t-statistic computed
            during the test follow their theoretically predicted distribution."
          ),
          br(),
          br(),
          plotlyOutput('iDEQQ')
        ),
        tabPanel(
          "Volcano Plot",
          br(),
          span(
            "Generated using R plotly.
            A volcano plot displays statistical significance
            (-log10 P value) versus magnitude of change (log2 fold change)
            and is useful for visualizing differentially expressed genes.
            Highlighted genes are significantly differentially expressed at
            the selected adjusted p-value cutoff value."
          ),
          br(),
          br(),
          plotlyOutput('iDEVolcano')
        ),
        tabPanel(
          "Mean Difference Plot",
          br(),
          span(
            "Generated using R plotly.
            A mean difference (MD) plot displays
            log2 fold change versus average
            log2 expression values and is useful for visualizing differentially
            expressed genes. Highlighted genes are significantly differentially
            expressed at the selected adjusted p-value cutoff."
          ),
          br(),
          br(),
          plotlyOutput('iDEMd')
        ),
        tabPanel(
          "Heatmap Plot",
          br(),
          span(
            "Generated using R heatmaply.
            A heatmap plot displaying the top differentially expressed genes
            expression values for each experimental condition. The expression
            values are displayed post KNN imputation, count per million
            transformation, log transformation, normalisation and limma
            precision weights if selected."
          ),
          br(),
          br(),
          numericInput(
            "numberOfGenes",
            "Input the number of genes to display:",
            2,
            min = 2,
            max = 250,
            step = 1
          ),
          br(),
          plotlyOutput('iHeatmap')
        )
      )
    )
  return(differentialGeneExpressionAnalysisUi)
}

#' A Function to Return the Exploratory Data Analysis Ui Component
#'
#' A Function to Return the Exploratory Data Analysis Ui Component
#' @import plotly
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceExploratoryDataAnalysisUi()
#' @author Guy Hunt
#' @noRd
sourceExploratoryDataAnalysisUi <- function() {
  exploratoryDataAnalysisUi <- tabPanel(
    "Exploratory Data Analysis",
    tabsetPanel(
      type = "tabs",
      tabPanel(
        "Expression Density Plot",
        br(),
        span(
          "Generated using R plotly.
          The plot below displays the distribution of the values
          of the genes in the dataset.
          This plot is useful for identifying if the data is normalised before
          performing differential expression analysis. If
          density curves are similar from gene
          to gene, it is indicative that the data
          is normalized and cross-comparable.
          The values are displayed post KNN imputation, count per million
          transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactiveDensityPlot')
      ),
      tabPanel(
        "3D Expression Density Plot",
        br(),
        span(
          "Generated using R plotly.
          The plot below displays the distribution of the values of the genes
          in the dataset. This plot is useful for identifying if the data is
          normalised before
          performing differential expression analysis.
          If density curves are similar from gene to gene, it is indicative
          that the data is normalized and cross-comparable. The values are
          displayed post KNN imputation, count per million
          transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactiveThreeDDensityPlot')
      ),
      tabPanel(
        "Box-and-Whisper Plot",
        br(),
        span(
          "Generated using R plotly.
          The plot below displays the distribution of the values
          of the genes in the dataset.
          The quartiles are calculated using the linear method.
          Viewing the distribution can be useful for determining if the
          data in the dataset is suitable for differential expression analysis.
          Generally, median-centred values are indicative that the data is
          normalized and cross-comparable.
          The values are displayed post KNN imputation, count per million
          transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactiveBoxAndWhiskerPlot'),
        plotOutput('nonInteractiveBoxAndWhiskerPlot')
      ),
      tabPanel(
        "PCA Scree Plot",
        br(),
        span(
          "Generated using R prcomp and plotly.
              Principal component analysis (PCA) reduces the
              dimensionality of multivariate data to two dimensions
              that can be visualized graphically with minimal loss
              of information."
        ),
        br(),
        span(
          "Eigenvalues correspond to the amount of the variation
              explained by each principal component (PC).
              The plot displays the eigenvalues against the number of
              dimensions.  The values are displayed post KNN imputation,
          count per million transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactivePcaScreePlot')
      ),
      tabPanel(
        "PCA Individuals Plot",
        br(),
        span(
          "Generated using R prcomp and R plotly.
              Principal component analysis (PCA)
              reduces the dimensionality of multivariate data to two dimensions
              that
              can be visualized graphically with minimal loss of information."
        ),
        br(),
        span(
          "Eigenvalues correspond to the amount of the variation explained
              by each principal component (PC). The plot displays the
              eigenvalues
              for each individual (row) in the gene expression dataset for the
              top two principal components (PC1 and PC2). The values are
              displayed post KNN imputation, count per million
              transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactivePcaIndividualsPlot')
      ),
      tabPanel(
        "Mean-Variance Plot",
        br(),
        span(
          "Generated using R limma and plotly.
          The plot below is used to check the mean-variance relationship
          of the expression data, after fitting a linear model.
          It can help show if there is a lot of variation in the data.
          Each point represents a gene.
          The values are displayed post KNN imputation, count per million
          transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactiveMeanVariancePlot')
      ),
      tabPanel(
        "Heatmap Plot",
        br(),
        span(
          "Generated using R cor and heatmaply.
          The plot below compares the correlation values of the samples in a
          heatmap. The values are displayed post KNN imputation, count per
          million transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactiveHeatMapPlot')
      ),
      tabPanel(
        "PCA Variables Plot",
        br(),
        span(
          "Generated using R prcomp and R plotly. Principal component
          analysis (PCA) reduces the dimensionality of multivariate data to two
          dimensions that can be visualized graphically with minimal loss
          of information."
        ),
        br(),
        span(
          "Eigenvalues correspond to the amount of the variation explained
              by each principal component (PC). The plot displays the
              eigenvalues for each variable (column) in the gene expression
              dataset for the top two principal components (PC1 and PC2).
              The values are displayed post KNN imputation, count per million
              transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactivePcaVariablesPlot')
      ),
      tabPanel(
        "3D PCA Variables Plot",
        br(),
        span(
          "Generated using R prcomp and R plotly. Principal component
          analysis
          (PCA) reduces the dimensionality of multivariate data to two
          dimensions that can be visualized graphically with minimal loss
          of information."
        ),
        br(),
        span(
          "Eigenvalues correspond to the amount of the variation explained
              by each principal component (PC). The plot displays the
              eigenvalues for each variable (column) in the gene expression
              dataset for the top three principal components
              (PC1, PC2 and PC3).
              The values are displayed post KNN imputation, count per million
              transformation and log transformation if selected."
        ),
        br(),
        br(),
        plotlyOutput('interactive3DPcaVariablesPlot')
      ),
      tabPanel(
        "UMAP Plot",
        br(),
        span(
          "Generated using R umap and plotly.
          Uniform Manifold Approximation and Projection (UMAP)
          is a dimension reduction technique useful for visualizing
          how genes are related to each other. The number of nearest
          neighbours used in the calculation is indicated in the graph.
          The values are displayed post KNN imputation, count per million
          transformation and log transformation if selected."
        ),
        br(),
        br(),
        numericInput(
          "knn",
          "Input the k-nearest neighbors value  to use:",
          2,
          min = 2,
          step = 1
        ),
        br(),
        plotlyOutput('interactiveUmapPlot')
      )
    )
  )
  return(exploratoryDataAnalysisUi)
}

#' A Function to Return the Side Bar Ui Component
#'
#' A Function to Return the Side Bar Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceSideBarUi()
#' @importFrom htmltools HTML
#' @author Guy Hunt
#' @noRd
sourceExampleUI <- function() {
  exampleUIComponents <- mainPanel(
    br(),
    br(),
    strong("This tab contains:"),
    helpText(
      "1. An example GEO accession code that can be
             loaded into GEOexplorer."
    ),
    helpText(
      "2. The required format of gene
             expression files to be processed by GEOexplorer."
    ),
    helpText(
      "2. An example of a microarray and an RNA seq gene expression
             file that can be processed by GEOexplorer."
    ),
    br(),
    br(),
    strong("Load an example GEO Accession Code"),
    br(),
    actionButton("loadExampleData", "Load"),
    br(),
    br(),
    strong("Download Gene Expression File Template"),
    br(),
    dataTableOutput('example'),
    downloadButton("downloadGeneExpressionFileTemplate", "Download"),
    br(),
    br(),
    strong("Download Example Microarray File"),
    br(),
    downloadButton("downloadMicroarrayExample", "Download"),
    br(),
    br(),
    strong("Download Example RNASeq File"),
    br(),
    downloadButton("downloadRnaSeqExample", "Download")
  )
  return(exampleUIComponents)
}

#' A Function to Return the Side Bar Ui Component
#'
#' A Function to Return the Side Bar Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceSideBarUi()
#' @importFrom DT dataTableOutput
#' @author Guy Hunt
#' @noRd
sourceGeoSearchUi <- function() {
  geoSearchUiComponents <- mainPanel(
    helpText("Please input a keyword or phrase (such as a paper
             title or author name) below, to search the GEO database for
             relevant datasets."),
    br(),
    br(),
    textInput("geoSearchTerm", "Keyword", value = ""),
    actionButton("searchGeo", "Search"),
    br(),
    br(),
    dataTableOutput('geoSearchResults')
  )
  return(geoSearchUiComponents)
}


#' A Function to Return the Enrichment Ui Component
#'
#' A Function to Return the Enrichment Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT dataTableOutput
#' @examples sourceEnrichmentnUi()
#' @author Guy Hunt
#' @noRd
sourceEnrichmentnUi <- function() {
  enrichmentUi <- tabPanel("Gene Enrichment Analysis",
                                   tabsetPanel(
                                     type = "tabs",
                                     tabPanel(
                                       "Set Parameters",
                                       br(),
                                       br(),
                                       selectInput("enrichDatabases",
                                                   "Select a database to use
                                                   for gene enrich analysis",
                                                   c("GO_Molecular_
                                                     Function_2015",
                                                     "GO_Cellular_
                                                     Component_2015",
                                                     "GO_Biological_
                                                     Process_2015"),
                                                   selected =
                                                     "GO_Biological
                                                   _Process_2015"),
                                       br(),
                                       br(),
                                       uiOutput("output101")
                                     ),
                                     tabPanel(
                                       "Differentially Expressed Genes Table",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       downloadButton(
                                         "downloadDifferentiallyExpressedGenesEnrichmentTable",
                                         "Download"),
                                       br(),
                                       br(),
                                       dataTableOutput(
                                         'differentiallyExpressedGenesEnrichmentTable')
                                     ),
                                     tabPanel(
                                       "Differentially Expressed Genes Plot",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       plotOutput('differentiallyExpressedGenesEnrichmentPlot')
                                     ),
                                     tabPanel(
                                       "Upregulated Genes Table",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       downloadButton(
                                         "downloadUpregulatedGenesEnrichmentTable",
                                         "Download"),
                                       br(),
                                       br(),
                                       dataTableOutput(
                                         'upregulatedGenesEnrichmentTable')
                                     ),
                                     tabPanel(
                                       "Upregulated Genes Plot",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       plotOutput(
                                         'upregulatedGenesEnrichmentPlot')
                                     ),
                                     tabPanel(
                                       "Downregulated Genes Table",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       downloadButton(
                                         "downloadDownregulatedGenesEnrichmentTable",
                                         "Download"),
                                       br(),
                                       br(),
                                       dataTableOutput(
                                         'downregulatedGenesEnrichmentTable')
                                     ),
                                     tabPanel(
                                       "Downregulated Genes Plot",
                                       br(),
                                       span(
                                         "Generated using R enrichR."
                                       ),
                                       br(),
                                       br(),
                                       plotOutput(
                                         'downregulatedGenesEnrichmentPlot')
                                     ),
                                   ))
  return(enrichmentUi)
}
