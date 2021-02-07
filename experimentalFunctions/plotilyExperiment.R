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

fig <- plot_ly(type = 'scatter3d', mode = 'lines', name = 'Density Plot')

data <- as.data.frame(knnDataInput)
data <- na.omit(data)

density1 <- density(data[,1])
density2 <- density(data[,2])
density3 <- density(data[,3])

#fig <- fig %>% add_trace(x = ~density1$x, y = ~density1$y, z = 1, name = names(data)[1])
#fig <- fig %>% add_trace(x = ~density2$x, y = ~density2$y, z = 2, name = names(data)[2])
#fig <- fig %>% add_trace(x = ~density3$x, y = ~density3$y, z = 3, name = names(data)[3])
#fig

#fig <- plot_ly(type = 'scatter', mode = 'lines', name = 'Density Plot')
#fig <- fig %>% add_trace(x = ~density1$x, y = ~density1$y, name = names(data)[1])
#fig <- fig %>% add_trace(x = ~density2$x, y = ~density2$y, name = names(data)[2])
#fig <- fig %>% add_trace(x = ~density3$x, y = ~density3$y, name = names(data)[3])
#fig

fig <- plot_ly(data, type = "scatter3d", mode = "lines") 

i <- 1
for(col in names(data)) {
  density <- density(data[,i])
  print(i)
  print(density)
  x <- density$x
  print(x)
  fig <- fig %>% add_trace(x = density$x, y = i, z = density$y, name = col)
  i <- i+1
}
fig <- fig %>% layout(xaxis = list(title = 'Intensity'), yaxis = list(title = 'Density'))
fig