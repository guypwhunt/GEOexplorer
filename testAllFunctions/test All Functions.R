# The purpose of this script is to test all the functions used in the shiny app
setwd('C:/Users/User/Documents/qmul_courses/ECS750PECS751PECS753PECS754PECS7500P - EECS MSC PROJECT - 202021/shiny_geo2r_visulisation')
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")

# Input Values
geoAccessionCode <- "GSE18380"
platform <- "GPL4694"
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "No"

# Get GEO2R data
gsetData <- getGeoData(geoAccessionCode, platform)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Apply log transformation to expression data if necessary
dataInput <- logTransformExpressionData(expressionData, logTransformation)

# Perform KNN transformation on log expression data if necessary
knnDataInput <- knnDataTransformation(dataInput, knnTransformation)

# Remove all incomplete rows
naOmitInput <- naOmitTransformation(knnDataInput)

# Perform PCA analysis on KNN transformation expression data
pcaDataInput <- pcaAnalysis(naOmitInput)

# Box-and-Whisker Plot
boxAndWhiskerPlot(geoAccessionCode, platform, knnDataInput)

# Expression Value Distribution Plot
expressionValueDistributionPlot(geoAccessionCode, platform, knnDataInput)

# Mean-Variance Plot
meanVariancePlot(geoAccessionCode, platform, naOmitInput)

# UMAP plot (multi-dimensional scaling)
umapPlot(geoAccessionCode, platform, naOmitInput)

# Principal component analysis scree plot
pcaScreePlot(pcaDataInput)

# Principal component analysis individuals plot
pcaIndividualsPlot(pcaDataInput)

# Principal component analysis variables plot
pcaVariablesPlot(pcaDataInput)

# Principal component analysis biplot of individuals and variables
pcaBiplotPlot(pcaDataInput)

# Interactive Box-and-Whisker Plot
library(plotly)
library(ggplot2)
fig <- interactiveBoxAndWhiskerPlot(knnDataInput, geoAccessionCode, platform                                 )
fig

# Interactive Density Plot
fig <- interactiveDesnityPlot(naOmitInput, geoAccessionCode, platform)
fig

# 3D Interactive Density Plot
fig <- interactiveThreeDDesnityPlot(naOmitInput, geoAccessionCode, platform)
fig

# Interactive UMAP
fig <- interactiveUmapPlot(naOmitInput)
fig