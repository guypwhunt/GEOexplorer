# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/guypw/OneDrive/Documents/geo2rShinyApp')

# Import Functions
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
source("differentialGeneExpressionAnalysis/differentialGeneExpressionAnalysis.R")
source("interactiveDataVisualizationFunctions/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R")


# Import Libraries
library(plotly)
library(ggplot2)

# Input Values
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "Yes"
knn <- 2
geoAccessionCodes <- list('GSE25770', 'GSE25771','GSE25774',"GSE18385","GSE18397", "GSE18400", "GSE18408", "GSE18423", "GSE18433", 'GSE25778', 'GSE25776', 'GSE25775', "GSE18388",'GSE25768','GSE25767','GSE25766','GSE25765','GSE25764','GSE25763','GSE25755','GSE25752','GSE25746','GSE25745','GSE25744','GSE25743','GSE25742','GSE25741','GSE25737','GSE25736','GSE25734','GSE25733','GSE25732','GSE25731','GSE25729','GSE25728','GSE25727','GSE25725','GSE25724','GSE25722','GSE25721','GSE18458','GSE18457','GSE18456','GSE18454','GSE18453','GSE18452','GSE18451','GSE18450','GSE18449','GSE18448','GSE18447','GSE18446','GSE18445','GSE18444','GSE18443','GSE18442','GSE18441','GSE18439','GSE18438','GSE18437','GSE18435','GSE18434','GSE18432','GSE18431','GSE18430','GSE18428','GSE18427','GSE18426','GSE18424','GSE18422','GSE18421','GSE18420','GSE18419','GSE18417','GSE18416','GSE18415','GSE18414','GSE18413','GSE18412','GSE18411','GSE18409','GSE18407','GSE18404','GSE18403','GSE18399','GSE18396','GSE18394','GSE18393','GSE18392','GSE18391','GSE18390','GSE18389','GSE18388','GSE18387','GSE18386','GSE18384','GSE18383','GSE18382','GSE18380')
pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
limmaPrecisionWeights <- "No"
forceNormalization <- "No"
platformAnnotation <- "NCBI generated"
significanceLevelCutOff <- 0.50
outputFile <-file("Output.txt") 

goodList <- list('GSE25772', 'GSE18388')
badList <- list("GSE25758", "GSE25762", "GSE25723", "GSE18459") # The first two have only 1 column, the third is just massive and the fourth errors on GEO2R
investigateList <- list()

for(geoAccessionCode in geoAccessionCodes)
{
#geoAccessionCode <- "GSE18388"
  tryCatch({

# Get the GEO2R data for all platforms
allGset <- getGset(geoAccessionCode)

# Get a list of all the platforms
platforms <- getPlatforms(allGset) 
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- getPlatformGset(allGset, platform)

# Get column Details
getColumnDetails <- getColumnDetails(gsetData)

# Extract the experiment information 
experimentInformation <- getExperimentInformation(gsetData)

# This function was retired
# Get GEO2R data
gsetData <- getGeoData(geoAccessionCode, platform)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Apply log transformation to expression data if necessary
#dataInput <- logTransformExpressionData(expressionData, logTransformation)

# Is log transformation auto applied
#autoLogInformation <- isLogTransformAutoApplied(expressionData)

# Perform KNN transformation on log expression data if necessary
#knnDataInput <- knnDataTransformation(dataInput, knnTransformation)

# Remove all incomplete rows
#naOmitInput <- naOmitTransformation(knnDataInput)

# This function was retired
# Perform PCA analysis on KNN transformation expression data
# pcaDataInput <- pcaAnalysis(naOmitInput)

# Perform PCA analysis on KNN transformation expression data
#pcaPrincompDataInput <- pcaPrincompAnalysis(naOmitInput)


# Differential gene expression analysis functions
# Get column names
columnNames <- extractColumns(expressionData)

# Define Groups
numberOfColumns <- length(columnNames)
numberOfColumns <- numberOfColumns + 1
halfNumberOfColumns <- ceiling(numberOfColumns/2)
i <- 0

group1 <- c() 
group2 <- c() 

for (name in columnNames) {
  if (i < halfNumberOfColumns) {
    group1 <- c(group1, name)
    i <- i +1
  } else {
    group2 <- c(group2, name)
    i <- i +1
  }
}

# Calculate gsms
gsms <- calculateGsms(columnNames,group1, group2)

# Convert P value adjustment
adjustment <- adjustmentCalculation(pValueAdjustment)

# Get fit2 data
#fit2 <- differentialGeneExpression(gsetData, knnDataInput, gsms, limmaPrecisionWeights, forceNormalization)

# Get fit 2 2
fit2 <- calculateFit2(geoAccessionCode, platform, gsms, logTransformation, limmaPrecisionWeights, forceNormalization, knnTransformation)

# Print Top deferentially expressed genes
tT <- topDifferentiallyExpressedGenesTable(fit2, adjustment)

# Plot Histogram of fit 2 data
fig <- histogramPlot(fit2, adjustment)

# Plot Interactive Histogram of fit 2 data
fig <- interactiveHistogramPlot(fit2, adjustment)
fig

# summarize test results as "up", "down" or "not expressed"
dT <- calculateDT(fit2, adjustment, significanceLevelCutOff)
ct <- 1 

# Venn diagram of results
fig <- vennDiagramPlot(dT)

# create Q-Q plot for t-statistic
fig <- qqPlot(fit2)

# Plot Interactive Q-Q plot for t-statistic
fig <- interactiveQQPlot(fit2, dT, ct)
fig

# volcano plot (log P-value vs log fold change)
fig <- volcanoPlot(fit2, dT, ct) 

# Interactive volcano plot (log P-value vs log fold change)
fig <- interactiveVolcanoPlot(fit2, dT, ct)
fig

# MD plot (log fold change vs mean log expression)
fig <- mdPlot(fit2, dT, ct)

# Plot Interactive Mean Difference of fit 2 data
fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
fig

t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqData <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic", plot.it = FALSE)

attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
final_attributes_list <- c()

for (attribute in attributes_list) {
  if (attribute %in% colnames(fit2$genes))
    final_attributes_list <- c(final_attributes_list, attribute)
}

qqData2 <- data.frame(qqData, dT[t.good,ct], fit2$genes[final_attributes_list][t.good,])

colnames(qqData2) <- c("x", "y", "regulation", final_attributes_list)
qqData2$regulation <- as.character(qqData2$regulation)
qqData2$regulation[qqData2$regulation == "1"] <- "Upregulated"
qqData2$regulation[qqData2$regulation == "0"] <- "Similar Expression"
qqData2$regulation[qqData2$regulation == "-1"] <- "Downregulation"

fig <- plot_ly()
fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"), 
                          hovertext = qqData2[final_attributes_list],
                          marker = list(size = 3))
fig <- fig %>% layout(
  title = ('Moderated t statistic'),
  xaxis = list(
    title = "Theoretical Quantiles"
  ),
  yaxis = list(
    title = "Sample Quantiles"
  ))
fig

  }, error = function(e) {
    write(as.character(geoAccessionCode),file = outputFile, append = TRUE, sep = "\n")
    close(outputFile)
  })
}