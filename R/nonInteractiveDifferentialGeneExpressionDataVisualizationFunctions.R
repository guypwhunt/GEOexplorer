
#' A Function to Plot a Venn Diagram with the Number of Genes that were and were not Differentially Expressed
#'
#' This function creates a venndigram containing the number of genes that were and were not differentially expressed
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @keywords GEO
#' @export
#' @examples fig <- nonInteractiveVennDiagramPlot(dT)
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object
nonInteractiveVennDiagramPlot <- function(dT) {
  library(limma)
  vennDiagram(dT, circle.col=palette())
}

#' A Function to Create a QQ Plot of the Quantiles of a Data Sample Against the Theoretical Quantiles of a Student's T Distribution from Differential Gene Expression Analysis
#'
#' This function allows you to plot a QQ plot of the quantiles of a data sample against the theoretical quantiles of a Student's t distribution from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @keywords GEO
#' @export
#' @examples fig <- nonInteractiveQQPlot(fit2)
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveQQPlot <- function(fit2) {
  library(limma)
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
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveVolcanoPlot <- function(fit2, dT, ct) {
  library(limma)
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
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
noninteractiveMeanDifferencePlot <- function(fit2, dT, ct) {
  library(limma)
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
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveHistogramPlot <- function(fit2, adjustment) {
  library(limma)
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
}
