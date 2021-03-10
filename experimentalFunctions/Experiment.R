# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

setwd('C:/Users/guypw/Documents/geo2rShinyApp')
source("geoIntegrationFunctions/geoIntegrationFunctions.R")
source("dataTransformationFunctions/dataTransformationFunctions.R")

logTransformation <- "Auto-Detect"  # Values can also be "Yes" or "No" 
knnTransformation <- "No" # Values can also be "Yes"
knn <- 2
outputFile <-file("output.txt")
geoAccessionCode <- "GSE18388"
pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
limmaPrecisionWeights <- "No"
forceNormalization <- "No"
platformAnnotation <- "NCBI generated"
significanceLevelCutOff <- 0.05
platform <- "GPL6246"
gsms <- "00001111"

fit2 <- calculateFit2(geoAccessionCode, platform, platformAnnotation, gsms, logTransformation, limmaPrecisionWeights, forceNormalization, knnTransformation)


#calculateFit2 <- function(geoAccessionCode, platform, gsms, logTransformation){
if (platformAnnotation == "Submitter supplied") {
  platformAnnotation <- FALSE
} else if (platformAnnotation == "NCBI generated") {
  platformAnnotation <- TRUE
} else {
  platformAnnotation <- TRUE
}
gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, AnnotGPL=platformAnnotation)
gset <- getPlatformGset(gset, platform)

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
sml <- strsplit(gsms, split="")[[1]]

ex <- extractExpressionData(gset)

ex <- logTransformExpressionData(ex, logTransformation)

ex <- knnDataTransformation(ex, knnTransformation)

sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
ex <- ex[ ,sel]
exprs(gset) <- ex

if(forceNormalization == "Yes"){
  exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
}



# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Group1","Group2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

if (limmaPrecisionWeights == "Yes"){
  nall <- nrow(gset)
  gset <- gset[complete.cases(exprs(gset)), ]
  
  # calculate precision weights and show plot of mean-variance trend
  v <- vooma(gset, design, plot=T)
  # OR weights by group
  # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
  v$genes <- fData(gset) # attach gene annotations
  
  # fit linear model
  fit  <- lmFit(v)
} else if (limmaPrecisionWeights == "No"){
  fit <- lmFit(gset, design)  # fit linear model 
}

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

fit2
#return(fit2)
#}