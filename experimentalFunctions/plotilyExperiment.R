  # The purpose of this script is to test all the functions used in the shiny app
  setwd('C:/Users/User/Documents/qmul_courses/ECS750PECS751PECS753PECS754PECS7500P - EECS MSC PROJECT - 202021/shiny_geo2r_visulisation')
  source("./geoIntegrationFunctions/geoIntegrationFunctions.R")
  source("./dataVisualizationFunctions/dataVisualizationFunctions.R")
  source("./dataTransformationFunctions/dataTransformationFunctions.R")
  source("./interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
  
  # Input Values
  geoAccessionCode <- "GSE18384"
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
  
  pcaDataInput <- pcaPrincompAnalysis(naOmitInput)
  
  
  library(plotly)
  library(dplyr)
  
  data <- pcaDataInput
  #pcaHC <- hclust(dist(data$scores),method = "ward.D2")
  #pcaClusters <- cutree(pcaHC,k=clusters)
  
  #pcaDf <- data.frame(data$scores,"cluster"=factor(pcaClusters))
  #pcaDf <- transform(pcaDf,cluster_name = paste("Cluster",pcaClusters))
  pcaDf <- data.frame(data$scores)
  pcaDf <- transform(pcaDf)
  
  fig <- plot_ly(pcaDf,x=~Comp.1,y=~Comp.2,text=rownames(pcaDf), mode="markers", type = 'scatter'
                 , marker = list(
                   color = 'rgb(17, 157, 255)',
                   size = 3,
                   line = list(
                     color = 'rgb(0, 0, 0)',
                     width = 1
                   ))
  )
  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Individuals Plot"),
                xaxis=list(title="PC1"),
                yaxis=list(title="PC2"))
  fig