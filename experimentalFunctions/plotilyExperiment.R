  # The purpose of this script is to test all the functions used in the shiny app
  setwd('C:/Users/User/Documents/qmul_courses/ECS750PECS751PECS753PECS754PECS7500P - EECS MSC PROJECT - 202021/shiny_geo2r_visulisation')
  source("./geoIntegrationFunctions/geoIntegrationFunctions.R")
  source("./dataVisualizationFunctions/dataVisualizationFunctions.R")
  source("./dataTransformationFunctions/dataTransformationFunctions.R")
  source("./interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
  
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
  
  pcaDataInput <- pcaAnalysis(naOmitInput)
  
  
  library(plotly)
  library(dplyr)
  
  # Perform PCA analysis on KNN transformation expression data
  pcaDataInput <- pcaPrincompAnalysis(naOmitInput)
  
  
  carsHC <- hclust(dist(pcaDataInput$scores),method = "ward.D2")
  carsClusters <- cutree(carsHC,k=3)

  
  carsDf <- data.frame(pcaDataInput$scores,"cluster"=factor(carsClusters))
  carsDf <- transform(carsDf,cluster_name = paste("Cluster",carsClusters))

  library(plotly)
  p <- plot_ly(carsDf,x=~Comp.1,y=~Comp.2,text=rownames(carsDf), mode="markers", type = 'scatter'
               ,color = ~cluster_name,marker=list(size=3)
  )
  p <- layout(p,title="PCA Clusters from Hierachical Clustering of Cars Data",
              xaxis=list(title="PC1"),
              yaxis=list(title="PC2"))
  p
  
  
  