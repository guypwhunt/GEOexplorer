#' A Function to Create a Box and Whisker Plot from an Expression Object
#'
#' This function allows you to plot expression data into a Box and Whisker Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- nonInteractiveBoxAndWhiskerPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveBoxAndWhiskerPlot <- function(geoAccessionCode = "", platform ="", ex) {
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
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- nonInteractiveDesnityPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveDesnityPlot <- function(geoAccessionCode = "", platform = "", ex) {
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
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- nonInteractiveMeanVariancePlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveMeanVariancePlot <- function(geoAccessionCode = "", platform = "", ex) {
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
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2
#' @examples fig <- nonInteractiveUmapPlot(expressionData, "GSE18380", "GPL4694")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
nonInteractiveUmapPlot <- function(geoAccessionCode = "", platform = "", ex, knn) {
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
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- nonInteractivePcaScreePlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaScreePlot <- function(pcaEx) {
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
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaIndividualsPlot <- function(pcaEx) {
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
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaVariablesPlot <- function(pcaEx) {
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
#' @param pcaEx A PCA object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaBiplotPlot <- function(pcaEx) {
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
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute umap maptools ggplot2 factoextra pheatmap
#' @examples fig <- nonInteractiveCorrelationMatrixPlot(expressionData)
#' @author Guy Hunt
nonInteractiveCorrelationMatrixPlot <- function(ex){
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

#' A Function to Plot a Venn Diagram with the Number of Genes that were and were not Differentially Expressed
#'
#' This function creates a venndigram containing the number of genes that were and were not differentially expressed
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples fig <- nonInteractiveVennDiagramPlot(dT)
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object
nonInteractiveVennDiagramPlot <- function(dT) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  vennDiagram(dT, circle.col=palette())
}

#' A Function to Create a QQ Plot of the Quantiles of a Data Sample Against the Theoretical Quantiles of a Student's T Distribution from Differential Gene Expression Analysis
#'
#' This function allows you to plot a QQ plot of the quantiles of a data sample against the theoretical quantiles of a Student's t distribution from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @keywords GEO
#' @export
#' @examples fig <- nonInteractiveQQPlot(fit2)
#' @import GEOquery limma umap data.table
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveQQPlot <- function(fit2) {
  # create Q-Q plot for t-statistic
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
}

#' A Function to Create a Volcano Plot of the Statistical Significance (-log10 P Value) Versus Magnitude of Change (log2 Fold Change) from Differential Gene Expression Analysis
#'
#' This function allows you to plot a volcano plot of the statistical significance (-log10 P value) versus magnitude of change (log2 fold change) from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @export
#' @examples fig <- nonInteractiveVolcanoPlot(fit2, dT, ct)
#' @import GEOquery limma umap data.table
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveVolcanoPlot <- function(fit2, dT, ct) {
  # volcano plot (log P-value vs log fold change)
  colnames(fit2) # list contrast names
  volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
              highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
}

#' A Function to Create a Mean Difference Plot of the log2 Fold Change Versus Average log2 Expression Values from Differential Gene Expression Analysis
#'
#' This function allows you to plot a mean difference plot of the log2 fold change versus average log2 expression values from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @export
#' @examples fig <- noninteractiveMeanDifferencePlot(fit2, dT, ct)
#' @import plotly ggplot2 limma scales
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
noninteractiveMeanDifferencePlot <- function(fit2, dT, ct) {
  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  abline(h=0)
}

#' A Function to Create a Histogram of the P values from Differential Gene Expression Analysis
#'
#' This function allows you to plot a histogram of the P values from differential gene expression analysis
#' @param fit2 An object containing the results of differntial gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A character string containing the adjustment to P-values
#' @keywords GEO
#' @export
#' @examples fig <- nonInteractiveHistogramPlot(fit2, adjustment)
#' @import plotly ggplot2 limma scales
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveHistogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
}
