#' A Function to Create a Box and Whisker Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Box and Whisker Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- boxAndWhiskerPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
boxAndWhiskerPlot <- function(geoAccessionCode = "", platform ="", ex) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  par(mar=c(7,4,2,1))
  title <- paste (geoAccessionCode, "/", platform, sep ="")
  boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
}

#' A Function to Create a Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Density Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- expressionValueDistributionPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
expressionValueDistributionPlot <- function(geoAccessionCode = "", platform = "", ex) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  par(mar=c(4,4,2,1))
  title <- paste(geoAccessionCode, "/", platform, " value distribution", sep ="")
  plotDensities(ex, main=title, legend=F)
}

#' A Function to Create a Mean Variance Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Mean Variance Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- meanVariancePlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
meanVariancePlot <- function(geoAccessionCode = "", platform = "", ex) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  plotSA(lmFit(ex), main= paste("Mean variance trend,", geoAccessionCode))
}

#' A Function to Create a UMAP Plot from an Expression Object
#'
#' This function allows you to plot expression data into a UMAP Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- umapPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
umapPlot <- function(geoAccessionCode = "", platform = "", ex, knn) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
  plot(ump$layout, main=paste("UMAP plot, number of nearest neighbors used =", knn), xlab="", ylab="", pch=20, cex=1.5)
  pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
}

#' A Function to Create a Histogram of the Principle Components from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into a Histogram of the Principle Components
#' @param pcaEx A PCA object which can be obtained from the pcaPrincompAnalysis() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- pcaScreePlot(pcaPrincompDataInput)
#' @author Guy Hunt
pcaScreePlot <- function(pcaEx) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  library(factoextra)
  library(pheatmap)
  fviz_eig(pcaEx)
}

#' A Function to Create an Individuals Scatter Plot from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Individuals Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the pcaPrincompAnalysis() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- pcaIndividualsPlot(pcaPrincompDataInput)
#' @author Guy Hunt
pcaIndividualsPlot <- function(pcaEx) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  library(factoextra)
  library(pheatmap)
  fviz_pca_ind(pcaEx,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               geom = "point",
               repel = TRUE     # Avoid text overlapping
  )
}

#' A Function to Create an Variables Scatter Plot from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Variables Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the pcaPrincompAnalysis() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- pcaVariablesPlot(pcaPrincompDataInput)
#' @author Guy Hunt
pcaVariablesPlot <- function(pcaEx) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  library(factoextra)
  library(pheatmap)
  fviz_pca_var(pcaEx,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
}

#' A Function to Create a Scatter Plot that contains both the Individuals and Variables from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Scatter Plot contains both the Individuals and Variables
#' @param pcaEx A PCA object which can be obtained from the pcaPrincompAnalysis() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- pcaBiplotPlot(pcaPrincompDataInput)
#' @author Guy Hunt
pcaBiplotPlot <- function(pcaEx) {
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  library(factoextra)
  library(pheatmap)
  fviz_pca_biplot(pcaEx, repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  geom = "point",
                  col.ind = "#696969",  # Individuals color
  )
}

#' A Function to Create a Correlation Matrix that contains both the Correlations between Samples
#'
#' This function allows you to plot PCA expression results into an Scatter Plot contains both the Individuals and Variables
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- correlationMatrixPlot(expressionData)
#' @author Guy Hunt
correlationMatrixPlot <- function(ex){
  library(limma)
  library(umap)
  library(maptools)
  library(ggplot2)
  library(factoextra)
  library(pheatmap)
  corMatrix <- cor(ex,use="c")
  fig <- pheatmap(corMatrix)
  return(fig)
}
