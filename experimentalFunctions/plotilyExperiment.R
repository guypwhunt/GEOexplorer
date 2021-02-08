  # The purpose of this script is to test all the functions used in the shiny app
  
  source(".\geoIntegrationFunctions/geoIntegrationFunctions.R")
  source(".\dataVisualizationFunctions/dataVisualizationFunctions.R")
  source(".\dataTransformationFunctions/dataTransformationFunctions.R")
  source(".\interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
  
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
  #library(umap)
  #library(maptools)
  #library(ggplot2)
  
  # Mean Varience Done
  data <- na.omit(knnDataInput)
  data <- lmFit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(data = data, x = ~Amean, y = ~sigma)
  fig
  
  # TBC
  # CREATE DATAFRAMES FROM THE DENSITY ATTRIBUTES
  # Scatter Plot Updates
  data <- na.omit(knnDataInput)
  data <- as.data.frame(data)
  density <- density(data)
  density <- as.data.frame(density)
  fig <- plot_ly(data = density, x = ~x, y = ~y, type = 'scatter', mode = 'lines', name = (paste(paste(geoAccessionCode,platform),'value distribution')))

  # Box plot updates
  data <- na.omit(knnDataInput)
  data <- as.data.frame(t(data))
  data <- as.data.frame(data)
  data$row.names <- row.names(data)
  attributes(data)
  
  fig <- plot_ly(data = data, x = ~row.names, y = data, type = "box", quartilemethod="linear")
  fig
  
  