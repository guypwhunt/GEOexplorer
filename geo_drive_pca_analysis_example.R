# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)
source("analyticsFunctions.R")

# load series and platform data from GEO

gset <- getGEO("GSE18384", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

pca_analysis2 <- function(ex){
  pca <- prcomp(ex, scale = TRUE)
  return(pca)
}

pca_data <- pca_analysis2(ex)
print(pca_data)

pca_data <- pca_analysis(ex)
print(pca_data)