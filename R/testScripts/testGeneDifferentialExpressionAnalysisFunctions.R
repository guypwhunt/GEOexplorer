# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/guypw/OneDrive/Documents/geo2rShinyApp')

# Import Functions
source("R/backendDifferentialGeneExpressionDataFunctions.R")
source("R/backendExploratoryDataAnalysisDataFunctions.R")
source("R/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R")
source("R/interactiveExploratoryDataAnalysisDataVisualizationFunctions.R")
source("R/nonInteractiveDifferentialGeneExpressionDataVisualizationFunctions.R")
source("R/nonInteractiveExploratoryDataAnalysisDataVisualizationFunctions.R")

# Import Libraries
library(plotly)
library(ggplot2)

# Input Values
logTransformations <- list("No", "Yes", "Auto-Detect")  # Values can also be "Yes" or "No"
knnTransformations <- list("Yes", "No") # Values can also be "Yes"
knn <- 2
pValueAdjustments <- list("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Holm", "None")
limmaPrecisionWeightss <- list("No", "Yes")
forceNormalizations <- list("No", "Yes")
platformAnnotations <- "NCBI generated"
significanceLevelCutOffs <- list(0.05, 0.2, 0.4, 0.6, 0.8, 1.0)



goodList <- list("GSE18385", "GSE18397", "GSE18400", "GSE18408", "GSE18423", "GSE18433", "GSE18434", "GSE18438","GSE18439","GSE18441", "GSE18442", "GSE18443","GSE18444","GSE18445","GSE18446","GSE18447","GSE18448", "GSE18449","GSE18450","GSE18451", "GSE25728", "GSE25734", "GSE25752", "GSE25721", "GSE25722", "GSE25724", "GSE25725", "GSE25727", "GSE25729", "GSE25731", "GSE25732", "GSE25736","GSE25737","GSE25741","GSE25742", "GSE25744","GSE25745","GSE25746", "GSE25755", "GSE25763","GSE25764","GSE25765","GSE25766","GSE25767","GSE25768","GSE25770","GSE25771","GSE25772","GSE25774","GSE25775","GSE25776","GSE25778","GSE25778", "GSE18380","GSE18382","GSE18383", "GSE18384", "GSE18386", "GSE18388","GSE18389","GSE18390","GSE18392","GSE18393","GSE18394","GSE18396", "GSE18399", "GSE18403","GSE18404","GSE18407", "GSE18409","GSE18411","GSE18412","GSE18413","GSE18414","GSE18415","GSE18416","GSE18417","GSE18419","GSE18420","GSE18421","GSE18422", "GSE18424","GSE18426","GSE18427","GSE18428","GSE18430","GSE18431","GSE18432","GSE18435","GSE18437", "GSE18452","GSE18453","GSE18454","GSE18456","GSE18457","GSE18458")
badList <- list("GSE25758", "GSE25762", "GSE25723", "GSE18459", "GSE25743", "GSE18387", "GSE18391", "GSE25733") # The first two have only 1 column, the third is just massive and the fourth errors on GEO2R

for(geoAccessionCode in goodList){
  for(logTransformation in logTransformations){
    for(knnTransformation in knnTransformations) {
      for(pValueAdjustment in pValueAdjustments) {
        for(limmaPrecisionWeights in limmaPrecisionWeightss){
          for(forceNormalization in forceNormalizations){
            for(significanceLevelCutOff in significanceLevelCutOffs) {
  tryCatch({
#geoAccessionCode <- "GSE25723"

# Get the GEO2R data for all platforms
allGset <- getGeoObject(geoAccessionCode)

# Get a list of all the platforms
platforms <- extractPlatforms(allGset)
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- extractPlatformGset(allGset, platform)

# Get column Details
getColumnDetails <- extractSampleDetails(gsetData)

# Extract the experiment information
experimentInformation <- extractExperimentInformation(gsetData)

# Get GEO2R data
gsetData <- extractGeoData(geoAccessionCode, platform)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Apply log transformation to expression data if necessary
expressionData <- calculateLogTransformation(expressionData, logTransformation)

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

ct <- 1

# Non-Interactive Venn diagram
fig <- nonInteractiveVennDiagramPlot(dT)

# Non-Interactive Q-Q plot
fig <- nonInteractiveQQPlot(fit2)

# Interactive Q-Q plot
fig <- interactiveQQPlot(fit2, dT, ct)
fig

# Non-Interactive volcano plot (log P-value vs log fold change)
fig <- nonInteractiveVolcanoPlot(fit2, dT, ct)

# Interactive volcano plot (log P-value vs log fold change)
fig <- interactiveVolcanoPlot(fit2, dT, ct)
fig

# Fix this
# MD plot (log fold change vs mean log expression)
fig <- noninteractiveMeanDifferencePlot(fit2, dT, ct)

# Plot Interactive Mean Difference of fit 2 data
fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
fig

  }, error = function(e) {
    outputFile <-file(paste0(geoAccessionCode,"DGEAOutput.txt"))
    write(as.character(geoAccessionCode),file = outputFile, append = TRUE, sep = "\n")
    close(outputFile)
  })
}
  }
}
  }
}
}
}
