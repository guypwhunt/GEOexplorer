library(GEOquery)
library(limma)
library(umap)
library("maptools")
library(ggplot2)

boxAndWhiskerPlot <- function(geoAccessionCode, platform, data) {
  par(mar=c(7,4,2,1))
  title <- paste (geoAccessionCode, "/", platform, sep ="")
  boxplot(data, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
}

expressionValueDistributionPlot <- function(geoAccessionCode, platform, data) {
  par(mar=c(4,4,2,1))
  title <- paste (geoAccessionCode, "/", platform, " value distribution", sep ="")
  plotDensities(data, main=title, legend=F)
}
  
meanVariancePlot <- function(geoAccessionCode, platform, data) {
  ex <- na.omit(data) # eliminate rows with NAs
  plotSA(lmFit(ex), main= paste("Mean variance trend,", geoAccessionCode))
}

umapPlot <- function(geoAccessionCode, platform, data) {
  ex <- data[!duplicated(data), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
  plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", pch=20, cex=1.5)
  pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
  }