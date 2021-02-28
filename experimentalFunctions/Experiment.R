# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/guypw/Documents/geo2rShinyApp')

# Import Functions
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
source("differentialGeneExpressionAnalysis/differentialGeneExpressionAnalysis.R")

library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE18388", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

attributes(gset)
attributes(phenoData(gset))

getColumnDetails <- function(gset){
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data
  df <- data.frame(
    column=row.names(phenoData),
    title=phenoData["title"],
                   source=phenoData["source_name_ch1"], 
                   characteristic1=phenoData["characteristics_ch1"], 
                   characteristic2=phenoData["characteristics_ch1.1"]) 
  return(df)
}
x <- getColumnDetails(gset)
x
