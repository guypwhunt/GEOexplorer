#library(GEOexplorer)

source("R/backendDifferentialGeneExpressionDataFunctions.R")
source("R/backendExploratoryDataAnalysisDataFunctions.R")
source("R/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R")
source("R/interactiveExploratoryDataAnalysisDataVisualizationFunctions.R")
source("R/nonInteractiveDifferentialGeneExpressionDataVisualizationFunctions.R")
source("R/nonInteractiveExploratoryDataAnalysisDataVisualizationFunctions.R")


# Input Values
logTransformation <- "Auto-Detect"  #
knnTransformation <- "Yes" # Values can also be "Yes"
knn <- 2
pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
limmaPrecisionWeights <- "Yes"
forceNormalization <- "Yes"
platformAnnotation <- "NCBI generated"
significanceLevelCutOff <- 0.05

# Get the GEO2R data for all platforms
geoAccessionCode <- "GSE18388"
allGset <- getGeoObject(geoAccessionCode)

# Extract platforms
platforms <- extractPlatforms(allGset)
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- extractPlatformGset(allGset, platform)

# Extract the experiment information
experimentInformation <- extractExperimentInformation(gsetData)

# Extract Sample Information
sampleInfo <- extractSampleInformation(gsetData)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Get column Details
columnInfo <- extractSampleDetails(gsetData)

# Is log transformation auto applied
autoLogInformation <- calculateAutoLogTransformApplication(expressionData)

# Get a list of all the columns
columns <- extractSampleNames(expressionData)

# Apply log transformation to expression data if necessary
dataInput <- calculateLogTransformation(expressionData, logTransformation)

# Perform KNN transformation on log expression data if necessary
knnDataInput <- calculateKnnImpute(dataInput, "Yes")

# Get a list of all the columns in the KNN output
knnColumns <- extractSampleNames(knnDataInput)

# Get knn output column Details
knnColumnInfo <- extractSampleDetails(gsetData)
knnColumnInfo <- knnColumnInfo[knnColumns,]

# Remove all incomplete rows
naOmitInput <- calculateNaOmit(knnDataInput)

# Perform Princomp PCA analysis on KNN transformation expression data
pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)


# Perform Prcomp PCA analysis on KNN transformation expression data
pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)

# Extract Experiment Information
extractedExperimentInformation <- convertExperimentInformation(experimentInformation)

# Non-Interactive Box-and-Whisker Plot
fig <- nonInteractiveBoxAndWhiskerPlot(ex = knnDataInput, geoAccessionCode = geoAccessionCode, platform = platform)

# Interactive Box-and-Whisker Plot
fig <- interactiveBoxAndWhiskerPlot(knnDataInput, geoAccessionCode, platform)
fig

# Non-Interactive Density Plot
fig <- nonInteractiveDensityPlot(ex = naOmitInput, geoAccessionCode = geoAccessionCode, platform = platform)

# Interactive Density Plot
fig <- interactiveDensityPlot(naOmitInput, geoAccessionCode, platform)
fig

# 3D Interactive Density Plot
fig <- interactiveThreeDDensityPlot(naOmitInput, geoAccessionCode, platform)
fig

# Interactive UMAP
fig <- interactiveUmapPlot(naOmitInput, knn, geoAccessionCode)
fig

# Interactive Mean Variance Plot
fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode, gsetData)
fig

# Interactive Princomp PCA Scree Plot
fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive Princomp PCA Individual Plot
fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput, geoAccessionCode, gsetData)
fig

# Interactive Princomp PCA Variables Plot
fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive Prcomp PCA Scree Plot
fig <- interactivePrcompPcaScreePlot(pcaPrcompDataInput, geoAccessionCode)
fig

# Interactive Prcomp PCA Individual Plot
fig <- interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput, geoAccessionCode, gsetData)
fig

# Interactive Prcomp PCA Variables Plot
fig <- interactivePrcompPcaVariablesPlot(pcaPrcompDataInput, geoAccessionCode)
fig

# Correlation Matrix of samples
fig <- interactiveHeatMapPlot(naOmitInput)
fig

# Non-Interactive UMAP
fig <- nonInteractiveUmapPlot(naOmitInput, knn, geoAccessionCode)

# Non-Interactive Mean Variance Plot
fig <- nonInteractiveMeanVariancePlot(naOmitInput, geoAccessionCode)

# Non-Interactive Princomp PCA Scree Plot
fig <- nonInteractivePcaScreePlot(pcaPrincompDataInput)
fig

# Non-Interactive Princomp PCA Individual Plot
fig <- nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
fig

# Non-Interactive Princomp PCA Variables Plot
fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
fig

# Non-Interactive Princomp PCA Individual and Variables Bilot
fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
fig

# Differential gene expression analysis functions
# Get column names
columnNames <- extractSampleNames(expressionData)

# Define Groups
numberOfColumns <- length(columnNames)
numberOfColumns <- numberOfColumns + 1
halfNumberOfColumns <- ceiling(numberOfColumns/2)
i <- 0

group1 <- c()
group2 <- c()

for (name in columnNames) {
  if (i < halfNumberOfColumns) {
    group1 <- c(group1, name)
    i <- i +1
  } else {
    group2 <- c(group2, name)
    i <- i +1
  }
}

# Select columns in group2
column2 <- calculateExclusiveColumns(columnNames, group1)

# Calculate gsms
gsms <- calculateEachGroupsSamples(columnNames,group1, group2)

# Convert P value adjustment
adjustment <- convertAdjustment(pValueAdjustment)

# Add error handling for subjects with 2 columns
# Get fit 2
fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, expressionData)

# Print Top deferentially expressed genes
tT <- calculateTopDifferentiallyExpressedGenes(fit2, adjustment)

# Non-Interactive Histogram
fig <- nonInteractiveHistogramPlot(fit2, adjustment)

# Interactive Histogram
fig <- interactiveHistogramPlot(fit2, adjustment)
fig

# summarize test results as "up", "down" or "not expressed"
dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)

# Non-Interactive Venn diagram
fig <- nonInteractiveVennDiagramPlot(dT)

# Non-Interactive Q-Q plot
fig <- nonInteractiveQQPlot(fit2)

# Interactive Q-Q plot
ct <- 1
fig <- interactiveQQPlot(fit2, dT, ct)
fig

# Non-Interactive volcano plot (log P-value vs log fold change)
# Unit testing not added
fig <- nonInteractiveVolcanoPlot(fit2, dT, ct)

# Interactive volcano plot (log P-value vs log fold change)
fig <- interactiveVolcanoPlot(fit2, dT, ct)
fig

# MD plot (log fold change vs mean log expression)
fig <- noninteractiveMeanDifferencePlot(fit2, dT, ct)

# Plot Interactive Mean Difference of fit 2 data
fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
fig
