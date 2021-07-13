#' A Function to Create a Box and Whisker Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Box and Whisker Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import limma
#' @examples fig <- nonInteractiveBoxAndWhiskerPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveBoxAndWhiskerPlot <- function(ex, geoAccessionCode = "", platform ="") {
  library(limma)
  par(mar=c(7,4,2,1))
  title <- paste (geoAccessionCode, "/", platform, sep ="")
  fig <- boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
  return(fig)
}

#' A Function to Create a Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Density Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import limma
#' @examples fig <- nonInteractiveDesnityPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveDesnityPlot <- function(ex, geoAccessionCode = "", platform = "") {
  library(limma)
  par(mar=c(4,4,2,1))
  title <- paste(geoAccessionCode, "/", platform, " value distribution", sep ="")
  fig <- plotDensities(ex, main=title, legend=F)
  return(fig)
}

#' A Function to Create a Mean Variance Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Mean Variance Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import limma
#' @examples fig <- nonInteractiveMeanVariancePlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveMeanVariancePlot <- function(ex, geoAccessionCode = "", platform = "") {
  library(limma)
  fig <- plotSA(lmFit(ex), main= paste("Mean variance trend,", geoAccessionCode))
  return(fig)
}

#' A Function to Create a UMAP Plot from an Expression Object
#'
#' This function allows you to plot expression data into a UMAP Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import umap limma
#' @importFrom maptools pointLabel
#' @examples fig <- nonInteractiveUmapPlot(expressionData, 3, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveUmapPlot <- function(ex, knn, geoAccessionCode = "", platform = "") {
  library(limma)
  library(umap)
  library(maptools)
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
  plot(ump$layout, main=paste("UMAP plot, number of nearest neighbors used =", knn), xlab="", ylab="", pch=20, cex=1.5)
  fig <- pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
  return(fig)
}

#' A Function to Create a Histogram of the Principle Components from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into a Histogram of the Principle Components
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @importFrom factoextra fviz_eig
#' @examples fig <- nonInteractivePcaScreePlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaScreePlot <- function(pcaEx) {
  library(factoextra)
  fig <- fviz_eig(pcaEx)
  return(fig)
}

#' A Function to Create an Individuals Scatter Plot from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Individuals Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @importFrom factoextra fviz_pca_ind
#' @examples fig <- nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaIndividualsPlot <- function(pcaEx) {
  library(factoextra)
  fig <- fviz_pca_ind(pcaEx,
               col.ind = "cos2", # Color by the quality of representation
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               geom = "point",
               repel = TRUE     # Avoid text overlapping
  )
  return(fig)
}

#' A Function to Create an Variables Scatter Plot from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Variables Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @importFrom factoextra fviz_pca_var
#' @examples fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaVariablesPlot <- function(pcaEx) {
  library(factoextra)
  fig <- fviz_pca_var(pcaEx,
               col.var = "contrib", # Color by contributions to the PC
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE     # Avoid text overlapping
  )
  return(fig)
}

#' A Function to Create a Scatter Plot that contains both the Individuals and Variables from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results into an Scatter Plot contains both the Individuals and Variables
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @importFrom factoextra fviz_pca_biplot
#' @examples fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaBiplotPlot <- function(pcaEx) {
  library(factoextra)
  fig <- fviz_pca_biplot(pcaEx, repel = TRUE,
                  col.var = "#2E9FDF", # Variables color
                  geom = "point",
                  col.ind = "#696969",  # Individuals color
  )
  return(fig)
}

#' A Function to Create a Correlation Matrix that contains both the Correlations between Samples
#'
#' This function allows you to plot a heatmap of the correlations between experimental conditions
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import pheatmap
#' @examples fig <- nonInteractiveCorrelationMatrixPlot(expressionData)
#' @author Guy Hunt
nonInteractiveCorrelationMatrixPlot <- function(ex){
  library(pheatmap)
  corMatrix <- cor(ex,use="c")
  fig <- pheatmap(corMatrix)
  return(fig)
}
