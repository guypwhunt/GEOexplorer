  # The purpose of this script is to test all the functions used in the shiny app
  setwd('C:/Users/User/Documents/qmul_courses/ECS750PECS751PECS753PECS754PECS7500P - EECS MSC PROJECT - 202021/shiny_geo2r_visulisation')
  source("./geoIntegrationFunctions/geoIntegrationFunctions.R")
  source("./dataVisualizationFunctions/dataVisualizationFunctions.R")
  source("./dataTransformationFunctions/dataTransformationFunctions.R")
  source("./interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
  
  library(scales)
  
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
  
  # First Graph
  data <- pcaDataInput
  
  eig.val <- get_eigenvalue(data)
  eig.val
  eig.val[1,2]
  eig.val[2,2]
  
  # Results for Variables
  res.var <- get_pca_var(data)
  res.var$coord          # Coordinates
  res.var$contrib[1]        # Contributions to the PCs
  res.var$cos2           # Quality of representation 
  # Results for individuals
  res.ind <- get_pca_ind(res.pca)
  res.ind$coord          # Coordinates
  res.ind$contrib        # Contributions to the PCs
  res.ind$cos2[,1]          # Quality of representation 
  
  pcaDf <- data.frame(data$scores)
  pcaDf <- transform(pcaDf)
  individualsStats <- get_pca_ind(res.pca)
  eigenValue <- get_eigenvalue(data)
  
  fig <- plot_ly(pcaDf,x=~Comp.1,y=~Comp.2,text=rownames(pcaDf), mode="markers", type = 'scatter'
                 , marker = list(
                   color = ~individualsStats$cos2[,1],
                   size = 3
                 )
  )
  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Individuals Plot"),
                xaxis=list(title=paste("PC1", label_percent()(eigenValue[1,2]/100))),
                yaxis=list(title=paste("PC2", label_percent()(eigenValue[2,2]/100))))
  fig