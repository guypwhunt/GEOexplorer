#' A Function to Return the Side Bar Ui Component
#'
#' A Function to Return the Side Bar Ui Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceSideBarUi()
#' @importFrom shinyBS bsTooltip
#' @importFrom shinybusy add_busy_spinner
#' @author Guy Hunt
#' @noRd
sourceSideBarUi <- function() {
  sideBarUi <- sidebarPanel(
    helpText(
      "Input a GEO series accession code (GSEXXXX format)
      to examine the gene expression data.
      This can be obtained from https://www.ncbi.nlm.nih.gov/gds."
    ),
    add_busy_spinner(spin = "fading-circle"),
    radioButtons(
      "inputDataType",
      label = "Do you want to upload a gene expression dataset or get one from
      GEO?",
      choices = list("GEO", "Upload"),
      selected = "GEO"
    ),
    uiOutput("dataInput"),
    uiOutput("dataInput2"),
    bsTooltip(
      id = "platform",
      title = "Each platform relates to a different microarray experiment
      performed in the study.",
      placement = "top",
      trigger = "hover"
    ),
    radioButtons(
      "logTransformation",
      label = "Apply log transformation to the data:",
      choices = list("Auto-Detect", "Yes", "No"),
      selected = "Auto-Detect"
    ),
    bsTooltip(
      id = "logTransformation",
      title = "The GEO database accepts a variety of data value types,
              including logged and unlogged data.
              Limma expects data values to be in log space.
              To address this, an auto-detect feature that checks
              the values of selected samples
              and automatically performs a log2 transformation on
              values determined not to be in log space.
              Alternatively, the user can select Yes
              to force log2 transformation,
              or No to override the auto-detect feature.
              The auto-detect feature only considers Sample values that
              have been assigned to a group, and applies the transformation in
              an all-or-none fashion",
      placement = "top",
      trigger = "hover"
    ),
    uiOutput("logTransformationText"),
    br(),
    radioButtons(
      "knnTransformation",
      label = "Apply k-nearest neighbors (KNN) algorithm to predict
      null data:",
      choices = list("Yes", "No"),
      selected = "No"
    ),
    bsTooltip(
      id = "knnTransformation",
      title = "Rows with over 50% missing values are imputed using the overall
              mean per sample. Columns with over
              80% will cause an error in the KNN
              computation.",
      placement = "top",
      trigger = "hover"
    ),
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
                                         experimental condition"
                                       ),
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
          dataTableOutput('knnColumnTable'),
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
              bsTooltip(
                id = "pValueAdjustment",
                title =
                  "P-value adjustments can be applied to reduce
                  the likelihood of
                  a false positive occurring. The P-value adjustment 'None'
                  indicates
                  no P-value adjustment will be applied and is the least
                  conservative
                  P-value adjustment. The Benjamini & Hochberg
                  (False discovery rate)
                  and Benjamini & Yekutieli methods are slightly more
                  conservative and aim to control the false discovery rate.
                  The Bonferroni
                  and Holm methods are the most conservative as they aim to
                  control the family-wise error rate.",
                placement = "top",
                trigger = "hover"
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
              ),
              bsTooltip(
                id = "limmaPrecisionWeights",
                title =
                  "Limma precision weights should be applied if there is
                  a strong
                  mean-variance trend as can be identified from the
                  'Mean-Variance Plot' tab. By applying limma precision
                  weights the
                  limma vooma function is used to estimate the mean-variance
                  relationship and uses this to compute appropriate
                  observational-level weights.",
                placement = "top",
                trigger = "hover"
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
              bsTooltip(
                id = "forceNormalization",
                title =
                  "Force normalisation should be selected if the gene
                  expression
                  dataset is not normally distributed, as can be identified
                  from the
                  'Box-and-Whisper Plot', the 'Expression Density Plot'
                  and the
                  '3D Expression Density Plot'. By selecting force
                  normalisation
                  quantile normalisation  is applied to the expression dataset
                  making all selected samples have identical value
                distributions.",
                placement = "top",
                trigger = "hover"
              ),
              sliderInput(
                "significanceLevelCutOff",
                "Significance level cut-off:",
                min = 0,
                max = 1,
                value = 0.05
              ),
              bsTooltip(
                id = "significanceLevelCutOff",
                title =
                  "The significance level cut-off is used to identify genes
                  that are
            differentially expressed between the two groups. Genes with
            adjusted P-values less than the significance level cut-off are
            determined to be differentially expressed.",
                placement = "top",
                trigger = "hover"
              ),
              uiOutput("differentialExpressionButton")
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
            the selected adjusted p-value cutoff value"
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
          The plot shows data after log and KNN transformation if they were
          performed."
        ),
        br(),
        br(),
        plotlyOutput('interactiveBoxAndWhiskerPlot')
      ),
      tabPanel(
        "Expression Density Plot",
        br(),
        span(
          "Generated using R plotly.
          The plot below displays the distribution of the values
          of the genes in the dataset.
          This plot complements the boxplot in
          checking for data normalization before
          differential expression analysis. If
          density curves are similar from gene
          to gene, it is indicative that the data
          is normalized and cross-comparable.
          The plot shows data after log and KNN
          transformation if they were performed."
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
          in the dataset. This plot complements the boxplot in checking
          for data normalization before differential expression analysis.
          If density curves are similar from gene to gene, it is indicative
          that the data is normalized and cross-comparable. The plot shows
          data after log and KNN transformation if they were performed."
        ),
        br(),
        br(),
        plotlyOutput('interactiveThreeDDensityPlot')
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
          The plot shows data after log and KNN transformation if they were
          performed."
        ),
        br(),
        br(),
        plotlyOutput('interactiveMeanVariancePlot')
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
          The plot shows data after log and KNN transformation
          if they were performed."
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
      ),
      tabPanel(
        "Heatmap Plot",
        br(),
        span(
          "Generated using R cor and heatmaply.
          The plot below compares the correlation values of the samples in a
          heatmap."
        ),
        br(),
        br(),
        plotlyOutput('interactiveHeatMapPlot')
      ),
      tabPanel(
        "PCA",
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Scree Plot",
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
              dimensions. The plot shows data after log and
              KNN transformation if they were performed."
            ),
            br(),
            br(),
            plotlyOutput('interactivePcaScreePlot')
          ),
          tabPanel(
            "Individuals Plot",
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
              top two principal components (PC1 and PC2). The plot shows data
              after log and KNN transformation if they were performed."
            ),
            br(),
            br(),
            plotlyOutput('interactivePcaIndividualsPlot')
          ),
          tabPanel(
            "Variables Plot",
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
              dataset for the top two principal components (PC1 and PC2).
              The plot shows data after log and KNN transformation if they were
              performed."
            ),
            br(),
            br(),
            plotlyOutput('interactivePcaVariablesPlot')
          )
        )
      )
    )
  )
  return(exploratoryDataAnalysisUi)
}
