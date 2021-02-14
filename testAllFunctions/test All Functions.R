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
geoAccessionCode <- "GSE18380"
platform <- "GPL4694"
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "No"
clusters <- 5
knn <- 5

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
fig <- interactiveUmapPlot(naOmitInput, knn)
fig

# Interactive Mean Variance Plot
fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode)
fig

# Interactive PCA Scree Plot
fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive PCA Individual Plot
fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput, geoAccessionCode, clusters)
fig

# Interactive PCA Variables Plot
fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput, geoAccessionCode)
fig
