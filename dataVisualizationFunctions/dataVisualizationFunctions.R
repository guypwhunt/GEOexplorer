library(limma)
library(umap)
library(maptools)
library(ggplot2)
library(factoextra)

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
  plotSA(lmFit(ex), main= paste("Mean variance trend,", geoAccessionCode))
}

umapPlot <- function(geoAccessionCode, platform, data) {
  ex <- data[!duplicated(data), ]  # remove duplicates
  nNeighbors <- 5
  ump <- umap(t(ex), n_neighbors = nNeighbors, random_state = 123)
  plot(ump$layout, main=paste("UMAP plot, number of nearest neighbors used =", nNeighbors), xlab="", ylab="", pch=20, cex=1.5)
  pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
}

pcaScreePlot <- function(data) {
  fviz_eig(data)
}

pcaIndividualsPlot <- function(data) {
  fviz_pca_ind(data,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               geom = "point",
               repel = TRUE     # Avoid text overlapping
  )
}

pcaVariablesPlot <- function(data) {
  fviz_pca_var(data,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
}

pcaBiplotPlot <- function(data) {
  fviz_pca_biplot(data, repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  geom = "point",
                  col.ind = "#696969",  # Individuals color
  )
}