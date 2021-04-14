library(limma)
library(umap)
library(maptools)
library(ggplot2)
library(factoextra)
library(pheatmap)


boxAndWhiskerPlot <- function(geoAccessionCode, platform, ex) {
  par(mar=c(7,4,2,1))
  title <- paste (geoAccessionCode, "/", platform, sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
}

expressionValueDistributionPlot <- function(geoAccessionCode, platform, ex) {
  par(mar=c(4,4,2,1))
  title <- paste(geoAccessionCode, "/", platform, " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)
}
  
meanVariancePlot <- function(geoAccessionCode, platform, ex) {
  plotSA(lmFit(ex), main= paste("Mean variance trend,", geoAccessionCode))
}

umapPlot <- function(geoAccessionCode, platform, ex, knn) {
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
  plot(ump$layout, main=paste("UMAP plot, number of nearest neighbors used =", knn), xlab="", ylab="", pch=20, cex=1.5)
  pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
}

pcaScreePlot <- function(ex) {
  fviz_eig(ex)
}

pcaIndividualsPlot <- function(ex) {
  fviz_pca_ind(ex,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               geom = "point",
               repel = TRUE     # Avoid text overlapping
  )
}

pcaVariablesPlot <- function(ex) {
  fviz_pca_var(ex,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
}

pcaBiplotPlot <- function(ex) {
  fviz_pca_biplot(ex, repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  geom = "point",
                  col.ind = "#696969",  # Individuals color
  )
}

correlationMatrixPlot <- function(ex){
  corMatrix <- cor(ex,use="c")
  fig <- pheatmap(corMatrix)  
  return(fig)
}