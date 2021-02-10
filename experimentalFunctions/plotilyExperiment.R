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
  
  library(plotly)
  library(dplyr)
  
  # Perform PCA analysis on KNN transformation expression data
  pcaDataInput <- pcaAnalysis(naOmitInput)
  attributes(pcaDataInput)
  
  pcaDataInput$scale
  
  x <- pcaDataInput$x
  df <- data.frame(x)
  colnames(df)
  
  fig <- plot_ly(
    name = "Scree Plot",
    type = "bar"
  )
  
  i = 1 
  for(col in colnames(df)) {
    fig <- fig %>% add_trace(x = col ,y = df[,i], name = col)
    i <- i +1
    
  } 
  
  fig
  
