# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/guypw/OneDrive/Documents/geo2rShinyApp')

# Import Functions
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
source("differentialGeneExpressionAnalysis/differentialGeneExpressionAnalysis.R")


# Import Libraries
library(plotly)
library(ggplot2)
library(stringr)

# Input Values
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No"
knnTransformation <- "No" # Values can also be "Yes"
knn <- 2
geoAccessionCodes <- list("GSE18380")
pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
limmaPrecisionWeights <- "No"
forceNormalization <- "No"
platformAnnotation <- "NCBI generated"
significanceLevelCutOff <- 0.05

goodList <- list("GSE18385", "GSE18397", "GSE18400", "GSE18408", "GSE18423", "GSE18433", "GSE18434", "GSE18438","GSE18439","GSE18441", "GSE18442", "GSE18443","GSE18444","GSE18445","GSE18446","GSE18447","GSE18448", "GSE18449","GSE18450","GSE18451", "GSE25728", "GSE25734", "GSE25752", "GSE25721", "GSE25722", "GSE25724", "GSE25725", "GSE25727", "GSE25729", "GSE25731", "GSE25732", "GSE25733", "GSE25736","GSE25737","GSE25741","GSE25742","GSE25743","GSE25744","GSE25745","GSE25746", "GSE25755", "GSE25763","GSE25764","GSE25765","GSE25766","GSE25767","GSE25768","GSE25770","GSE25771","GSE25772","GSE25774","GSE25775","GSE25776","GSE25778","GSE25778", "GSE18380","GSE18382","GSE18383", "GSE18384", "GSE18386","GSE18387","GSE18388","GSE18389","GSE18390","GSE18391","GSE18392","GSE18393","GSE18394","GSE18396", "GSE18399", "GSE18403","GSE18404","GSE18407", "GSE18409","GSE18411","GSE18412","GSE18413","GSE18414","GSE18415","GSE18416","GSE18417","GSE18419","GSE18420","GSE18421","GSE18422", "GSE18424","GSE18426","GSE18427","GSE18428","GSE18430","GSE18431","GSE18432","GSE18435","GSE18437", "GSE18452","GSE18453","GSE18454","GSE18456","GSE18457","GSE18458")
badList <- list("GSE25758", "GSE25762", "GSE25723", "GSE18459") # The first two have only 1 column, the third is just massive and the fourth errors on GEO2R

#for(geoAccessionCode in geoAccessionCodes)
#{
geoAccessionCode <- "GSE50499"
#  tryCatch({

# Get the GEO2R data for all platforms
allGset <- getGset(geoAccessionCode)
attributes(allGset)


# Get a list of all the platforms
platforms <- getPlatforms(allGset)
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- getPlatformGset(allGset, platform)
gsetData

# Extract the experiment information
experimentInformation <- getExperimentInformation(gsetData)
experimentInformation

# Extract Sample Information
sampleInfo <- extractSampleInfo(gsetData)
sampleInfo

# Extract Sample Information
geneAnnotation <- extractGeneAnnotation(gsetData)
geneAnnotation

# Extract expression data
expressionData <- extractExpressionData(gsetData)
expressionData

# Get column Details
columnInfo <- getColumnDetails(gsetData)
columnInfo

# Is log transformation auto applied
autoLogInformation <- isLogTransformAutoApplied(expressionData)
autoLogInformation

# Get a list of all the columns
columns <- extractColumns(expressionData)
columns

# Apply log transformation to expression data if necessary
dataInput <- logTransformExpressionData(expressionData, logTransformation)
dataInput

# Perform KNN transformation on log expression data if necessary
knnDataInput <- knnDataTransformation(dataInput, knnTransformation)

# Remove all incomplete rows
naOmitInput <- naOmitTransformation(knnDataInput)

# Perform PCA analysis on KNN transformation expression data
pcaPrincompDataInput <- pcaPrincompAnalysis(naOmitInput)

# Extract Experiment Information
extractedExperimentInformation <- extractExperimentInformation(experimentInformation)

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
fig <- interactiveUmapPlot(naOmitInput, knn, geoAccessionCode)
fig

# Interactive Mean Variance Plot
fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode, gsetData)
fig

# Interactive PCA Scree Plot
fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Interactive PCA Individual Plot
fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput, geoAccessionCode, gsetData)
fig

# Interactive PCA Variables Plot
fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput, geoAccessionCode)
fig

# Correlation Matrix of samples
fig <- interactiveHeatMapPlot(naOmitInput)
fig

##########################
ex <- lmFit(naOmitInput)
ex <- as.data.frame(ex)
ex["ID"] <- rownames(ex)

colnames(gsetData@featureData)

gsetData@featureData@data["Gene symbol"]

geneData <- gsetData@featureData@data
geneData <- as.data.frame(geneData)
combineData <- merge(ex, geneData, by = "ID")
combineData %>% filter(ID %in% c(rownames(ex)))
colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")

if('ID' %in% colnames(combineData)){
  if('Gene symbol' %in% colnames(combineData)){
    if('Gene title' %in% colnames(combineData)){
      if('Gene ID' %in% colnames(combineData)){
        fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene symbol, '<br></br>', 'Gene Title: ', Gene title, '<br></br>', 'Gene ID: ', Gene ID, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                       hoverinfo = text,
                       mode = 'markers', marker = list(
                         color = 'rgb(17, 157, 255)',
                         size = 3,
                         line = list(
                           color = 'rgb(0, 0, 0)',
                           width = 1
                         )))
      } else {
        fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                       hoverinfo = text,
                       mode = 'markers', marker = list(
                         color = 'rgb(17, 157, 255)',
                         size = 3,
                         line = list(
                           color = 'rgb(0, 0, 0)',
                           width = 1
                         )))
      }
    } else {
      fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                     text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                     hoverinfo = text,
                     mode = 'markers', marker = list(
                       color = 'rgb(17, 157, 255)',
                       size = 3,
                       line = list(
                         color = 'rgb(0, 0, 0)',
                         width = 1
                       )))
    }
  } else{
    fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                   text = ~paste('ID: ', ID, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                   hoverinfo = text,
                   mode = 'markers', marker = list(
                     color = 'rgb(17, 157, 255)',
                     size = 3,
                     line = list(
                       color = 'rgb(0, 0, 0)',
                       width = 1
                     )))
  }
} else{
  fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                 text = ~paste('Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                 hoverinfo = text,
                 mode = 'markers', marker = list(
                   color = 'rgb(17, 157, 255)',
                   size = 3,
                   line = list(
                     color = 'rgb(0, 0, 0)',
                     width = 1
                   )))
}
fig <- fig %>% layout(
  title = (paste('Mean variance trend, ',geoAccessionCode)))
fig
##################

#  }, error = function(e) {
#    outputFile <-file("output.txt")
#    write(as.character(geoAccessionCode),file = outputFile, append = TRUE, sep = "\n")
#    close(outputFile)
#  })
#}
