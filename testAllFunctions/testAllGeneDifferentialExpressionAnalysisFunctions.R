# The purpose of this script is to test all the functions used in the shiny app

# Change the below file path to the file path you save the repo to
setwd('C:/Users/guypw/Documents/geo2rShinyApp')

# Import Functions
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataVisualizationFunctions/dataVisualizationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")
source("interactiveDataVisualizationFunctions/interactiveDataVisualizationFunctions.R")
source("differentialGeneExpressionAnalysis/differentialGeneExpressionAnalysis.R")


# Import Libraries
library(plotly)
library(ggplot2)

# Input Values
logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "Yes"
knn <- 2
outputFile <-file("output.txt")
geoAccessionCodes <- list("GSE18388")
pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
limmaPrecisionWeights <- "No"
forceNormalization <- "No"
platformAnnotation <- "NCBI generated"
significanceLevelCutOff <- 0.05

goodList <- list("GSE18385", "GSE18397", "GSE18400", "GSE18408", "GSE18423", "GSE18433", "GSE18434", "GSE18438","GSE18439","GSE18441", "GSE18442", "GSE18443","GSE18444","GSE18445","GSE18446","GSE18447","GSE18448", "GSE18449","GSE18450","GSE18451", "GSE25728", "GSE25734", "GSE25752", "GSE25721", "GSE25722", "GSE25724", "GSE25725", "GSE25727", "GSE25729", "GSE25731", "GSE25732", "GSE25733", "GSE25736","GSE25737","GSE25741","GSE25742","GSE25743","GSE25744","GSE25745","GSE25746", "GSE25755", "GSE25763","GSE25764","GSE25765","GSE25766","GSE25767","GSE25768","GSE25770","GSE25771","GSE25772","GSE25774","GSE25775","GSE25776","GSE25778","GSE25778", "GSE18380","GSE18382","GSE18383", "GSE18384", "GSE18386","GSE18387","GSE18388","GSE18389","GSE18390","GSE18391","GSE18392","GSE18393","GSE18394","GSE18396", "GSE18399", "GSE18403","GSE18404","GSE18407", "GSE18409","GSE18411","GSE18412","GSE18413","GSE18414","GSE18415","GSE18416","GSE18417","GSE18419","GSE18420","GSE18421","GSE18422", "GSE18424","GSE18426","GSE18427","GSE18428","GSE18430","GSE18431","GSE18432","GSE18435","GSE18437", "GSE18452","GSE18453","GSE18454","GSE18456","GSE18457","GSE18458")
badList <- list("GSE25758", "GSE25762", "GSE25723", "GSE18459") # The first two have only 1 column, the third is just massive and the fourth errors on GEO2R

#for(geoAccessionCode in geoAccessionCodes)
#{
geoAccessionCode <- "GSE18397"#"GSE18388"
#  tryCatch({

# Get the GEO2R data for all platforms
allGset <- getGset(geoAccessionCode, platformAnnotation)

# Get a list of all the platforms
platforms <- getPlatforms(allGset) 
platform <- platforms[1]

# Extract the GEO2R data from the specified platform
gsetData <- getPlatformGset(allGset, platform)

# Get column Details
getColumnDetails <- getColumnDetails(gsetData)
getColumnDetails

# Extract the experiment information 
experimentInformation <- getExperimentInformation(gsetData)

# This function was retired
# Get GEO2R data
#gsetData <- getGeoData(geoAccessionCode, platform)

# Extract expression data
expressionData <- extractExpressionData(gsetData)

# Apply log transformation to expression data if necessary
dataInput <- logTransformExpressionData(expressionData, logTransformation)

# Is log transformation auto applied
autoLogInformation <- isLogTransformAutoApplied(expressionData)

# Perform KNN transformation on log expression data if necessary
knnDataInput <- knnDataTransformation(dataInput, knnTransformation)

# Remove all incomplete rows
naOmitInput <- naOmitTransformation(knnDataInput)

# This function was retired
# Perform PCA analysis on KNN transformation expression data
# pcaDataInput <- pcaAnalysis(naOmitInput)

# Perform PCA analysis on KNN transformation expression data
pcaPrincompDataInput <- pcaPrincompAnalysis(naOmitInput)


# Differential gene expression analysis functions
# Get column names
columnNames <- extractColumns(expressionData)

# Define Groups
group1 <- list ("GSM458594", "GSM458595", "GSM458596") #, "GSM458597") 
group2 <- list ("GSM458598", "GSM458599", "GSM458600", "GSM458601") 

# Calculate gsms
gsms <- calculateGsms(columnNames,group1, group2)

# Convert P value adjustment
adjustment <- adjustmentCalculation(pValueAdjustment)

# Get annotated GSET
dEAllGset <- getGset(geoAccessionCode, platformAnnotation, GSEMatrix=TRUE, getGPL=TRUE)
dEGsetData <- getPlatformGset(dEAllGset, platform)

# Extract expression data
dEExpressionData <- extractExpressionData(dEGsetData)

# Apply log transformation to expression data if necessary
dEDataInput <- logTransformExpressionData(dEExpressionData, logTransformation)

# Perform KNN transformation on log expression data if necessary
dEKnnDataInput <- knnDataTransformation(dEDataInput, knnTransformation)



  # make proper column names to match toptable 
  fvarLabels(dEGsetData) <- make.names(fvarLabels(dEGsetData))
  
  # group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]
  
  # filter out excluded samples (marked as "X")
  sel <- which(sml != "X")
  sml <- sml[sel]
  dEGsetData <- dEGsetData[ ,sel]
  dEKnnDataInput <- dEKnnDataInput[ ,sel]
  exprs(dEGsetData) <- dEKnnDataInput
  
  if(forceNormalization == "Yes"){
    exprs(dEGsetData) <- normalizeBetweenArrays(exprs(dEGsetData)) # normalize data
  }
  
  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("Group1","Group2"))
  levels(gs) <- groups
  dEGsetData$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)
  
  if (limmaPrecisionWeights == "Yes"){
    nall <- nrow(dEGsetData)
    dEGsetData <- dEGsetData[complete.cases(exprs(dEGsetData)), ]
    
    # calculate precision weights and show plot of mean-variance trend
    v <- vooma(dEGsetData, design, plot=T)
    # OR weights by group
    # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
    v$genes <- fData(dEGsetData) # attach gene annotations
    
    # fit linear model
    fit  <- lmFit(v)
  } else if (limmaPrecisionWeights == "No"){
    fit <- lmFit(dEGsetData, design)  # fit linear model 
  }
  
  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  return(fit2)
}


# Get fit2 data
fit2 <- differentialGeneExpression1(dEGsetData, dEKnnDataInput, gsms, limmaPrecisionWeights, forceNormalization)

# Print Top deferentially expressed genes
tT <- topDifferentiallyExpressedGenesTable(fit2, adjustment)

# Plot Histogram of fit 2 data
fig <- histogramPlot(fit2, adjustment)
fig

# summarize test results as "up", "down" or "not expressed"
dT <- dT(fit2, adjustment, significanceLevelCutOff)

# Venn diagram of results
fig <- vennDiagramPlot(dT)
fig

# create Q-Q plot for t-statistic
fig <- qqPlot(fit2)
fig

# volcano plot (log P-value vs log fold change)
ct <- 1 
fig <- volcanoPlot(fit2, dT, ct) 
fig

# MD plot (log fold change vs mean log expression)
fig <- mdPlot(fit2, dT, ct)
fig

#  }, error = function(e) {
#    write(as.character(geoAccessionCode),file = outputFile, append = TRUE, sep = "\n")
#    close(outputFile)
#  })
#}