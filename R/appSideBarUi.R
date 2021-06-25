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
