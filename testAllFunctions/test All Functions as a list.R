# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/User/Documents/qmul_courses/ECS750PECS751PECS753PECS754PECS7500P - EECS MSC PROJECT - 202021/shiny_geo2r_visulisation')

# Import Functions
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")

# Import Libraries
library(plotly)
library(ggplot2)

# Input Values
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "Yes"
knn <- 2
outputFile <-file("output.txt")
geoAccessionCodes = list()

goodList <- list("GSE18385", "GSE18397", "GSE18400", "GSE18408", "GSE18423", "GSE18433", "GSE18434", "GSE18438","GSE18439","GSE18441", "GSE18442", "GSE18443","GSE18444","GSE18445","GSE18446","GSE18447","GSE18448", "GSE18449","GSE18450","GSE18451", "GSE25728", "GSE25734", "GSE25752", "GSE25721", "GSE25722", "GSE25724", "GSE25725", "GSE25727", "GSE25729", "GSE25731", "GSE25732", "GSE25733", "GSE25736","GSE25737","GSE25741","GSE25742","GSE25743","GSE25744","GSE25745","GSE25746", "GSE25755", "GSE25763","GSE25764","GSE25765","GSE25766","GSE25767","GSE25768","GSE25770","GSE25771","GSE25772","GSE25774","GSE25775","GSE25776","GSE25778","GSE25778", "GSE18380","GSE18382","GSE18383", "GSE18384", "GSE18386","GSE18387","GSE18388","GSE18389","GSE18390","GSE18391","GSE18392","GSE18393","GSE18394","GSE18396", "GSE18399", "GSE18403","GSE18404","GSE18407", "GSE18409","GSE18411","GSE18412","GSE18413","GSE18414","GSE18415","GSE18416","GSE18417","GSE18419","GSE18420","GSE18421","GSE18422", "GSE18424","GSE18426","GSE18427","GSE18428","GSE18430","GSE18431","GSE18432","GSE18435","GSE18437", "GSE18452","GSE18453","GSE18454","GSE18456","GSE18457","GSE18458")
badList <- list("GSE25758", "GSE25762", "GSE25723", "GSE18459") # The first two have only 1 column, the third is just massive and the fourth errors on GEO2R

for(geoAccessionCode in goodList)
{
  geoAccessionCode <- geoAccessionCode
  tryCatch({

# Get the GEO2R data for all platforms
allGset <- getGset(geoAccessionCode)

# Get a list of all the platforms
platforms <- getPlatforms(allGset) 
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- getPlatformGset(allGset, platform)

# Extract the experiment information 
experimentInformation <- getExperimentInformation(gsetData)

# This function was retired
# Get GEO2R data
#gsetData <- getGeoData(geoAccessionCode, platform)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Apply log transformation to expression data if necessary
dataInput <- logTransformExpressionData(expressionData, logTransformation)

# Is log transformation auto applied
autoLogInformation <- isLogTransformAutoApplied(expressionData)

# Perform KNN transformation on log expression data if necessary
knnDataInput <- knnDataTransformation(dataInput, knnTransformation)

# Remove all incomplete rows
naOmitInput <- naOmitTransformation(knnDataInput)

# Perform PCA analysis on KNN transformation expression data
pcaDataInput <- pcaAnalysis(naOmitInput)

# Perform PCA analysis on KNN transformation expression data
pcaPrincompDataInput <- pcaPrincompAnalysis(naOmitInput)

# Box-and-Whisker Plot
boxAndWhiskerPlot(geoAccessionCode, platform, knnDataInput)

# Expression Value Distribution Plot
expressionValueDistributionPlot(geoAccessionCode, platform, knnDataInput)

# Mean-Variance Plot
meanVariancePlot(geoAccessionCode, platform, naOmitInput)

# UMAP plot (multi-dimensional scaling)
umapPlot(geoAccessionCode, platform, naOmitInput, knn)

# Principal component analysis scree plot
pcaScreePlot(pcaPrincompDataInput)

# Principal component analysis individuals plot
pcaIndividualsPlot(pcaPrincompDataInput)

# Principal component analysis variables plot
pcaVariablesPlot(pcaPrincompDataInput)

# Principal component analysis biplot of individuals and variables
pcaBiplotPlot(pcaPrincompDataInput)

# Interactive Box-and-Whisker Plot
fig <- interactiveBoxAndWhiskerPlot(knnDataInput, geoAccessionCode, platform                                 )
fig

# Interactive Density Plot
fig <- interactiveDesnityPlot(naOmitInput, geoAccessionCode, platform)
fig

# 3D Interactive Density Plot
fig <- interactiveThreeDDesnityPlot(naOmitInput, geoAccessionCode, platform)
fig

# Interactive UMAP
fig <- interactiveUmapPlot(naOmitInput, knn, geoAccessionCode)
fig

# Interactive Mean Variance Plot
fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode)
fig

# Interactive PCA Scree Plot
fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive PCA Individual Plot
fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive PCA Variables Plot
fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput, geoAccessionCode)
fig

  }, error = function(e) {
    print(e);
    write(as.character(geoAccessionCode),file = outputFile, append = TRUE, sep = "\n")
    close(outputFile)
  })
}