  # The purpose of this script is to test all the functions used in the shiny app
  
  source(".\geoIntegrationFunctions/geoIntegrationFunctions.R")
  source(".\dataVisualizationFunctions/dataVisualizationFunctions.R")
  source(".\dataTransformationFunctions/dataTransformationFunctions.R")
  source(".\interactiveDataVisualizationFunctions\interactiveDataVisualizationFunctions.R")
  
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
  
  library(plotly)
  library(dplyr)
  
  knnDataInput <- as.data.frame(knnDataInput)
  
  knnDataInput <- na.omit(knnDataInput)
  ex <- data[!duplicated(knnDataInput), ]  # remove duplicates
  nNeighbors <- 5
  ump <- umap(t(ex), n_neighbors = nNeighbors, random_state = 123)
  
  attributes(ump)
  print(ump$layout)
  attributes(ump$layout)
  ump$layout
  ump$layout[1,]
  ump$layout[1,][1]
  row.names(ump$layout)
  
  library(plotly)
  i <- 1
  fig <- plot_ly(type = 'scatter', mode = 'markers')
  for(row in row.names(ump$layout)){
  fig <- fig %>% add_trace(x = ump$layout[i,][1], y = ump$layout[i,][2], name = row)
  i <- i+1
  }
  fig
  warnings()