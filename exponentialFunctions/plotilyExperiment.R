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

knnDataInput <- na.omit(knnDataInput)
data <- as.data.frame(knnDataInput)

fig <- plot_ly(data, type = "box", quartilemethod="linear")

i = 1
for(col in names(data)) {
  print(col)
  fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
  i <- i+1
}

fig

warnings()