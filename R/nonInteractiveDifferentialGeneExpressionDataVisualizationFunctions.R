
#' A Function to Plot a Venn Diagram with the Number of Genes that were and were not Differentially Expressed
#'
#' This function creates a venndigram containing the number of genes that were and were not differentially expressed
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @keywords GEO
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#'
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#' group1 <- c(group1, name)
#' i <- i +1
#' } else {
#' group2 <- c(group2, name)
#' i <- i +1
#' }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
#'
#' # Non-Interactive Venn diagram
#' fig <- nonInteractiveVennDiagramPlot(dT)
#'
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object
nonInteractiveVennDiagramPlot <- function(dT) {
  fig <- vennDiagram(dT, circle.col=palette())
  return(fig)
}

#' A Function to Create a QQ Plot of the Quantiles of a Data Sample Against the Theoretical Quantiles of a Student's T Distribution from Differential Gene Expression Analysis
#'
#' This function allows you to plot a QQ plot of the quantiles of a data sample against the theoretical quantiles of a Student's t distribution from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @keywords GEO
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#'
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#'   if (i < halfNumberOfColumns) {
#'     group1 <- c(group1, name)
#'     i <- i +1
#' } else {
#'     group2 <- c(group2, name)
#'     i <- i +1
#'   }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Non-Interactive Q-Q plot
#' fig <- nonInteractiveQQPlot(fit2)
#'
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveQQPlot <- function(fit2) {
  # create Q-Q plot for t-statistic
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  fig <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
  return(fig)
}

#' A Function to Create a Volcano Plot of the Statistical Significance (-log10 P Value) Versus Magnitude of Change (log2 Fold Change) from Differential Gene Expression Analysis
#'
#' This function allows you to plot a volcano plot of the statistical significance (-log10 P value) versus magnitude of change (log2 fold change) from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#'
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#'   if (i < halfNumberOfColumns) {
#'     group1 <- c(group1, name)
#'     i <- i +1
#'   } else {
#'     group2 <- c(group2, name)
#'     i <- i +1
#'   }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
#' ct <- 1
#'
#' # Non-Interactive volcano plot (log P-value vs log fold change)
#'fig <- nonInteractiveVolcanoPlot(fit2, dT, ct)
#'
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveVolcanoPlot <- function(fit2, dT, ct) {
  # volcano plot (log P-value vs log fold change)
  colnames(fit2) # list contrast names
  fig <- volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
              highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
  return(fig)
}

#' A Function to Create a Mean Difference Plot of the log2 Fold Change Versus Average log2 Expression Values from Differential Gene Expression Analysis
#'
#' This function allows you to plot a mean difference plot of the log2 fold change versus average log2 expression values from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#'
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#'   if (i < halfNumberOfColumns) {
#'     group1 <- c(group1, name)
#'     i <- i +1
#'   } else {
#'     group2 <- c(group2, name)
#'     i <- i +1
#'   }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
#' ct <- 1
#'
#' # MD plot (log fold change vs mean log expression)
#' fig <- noninteractiveMeanDifferencePlot(fit2, dT, ct)
#'
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpressionSummary()] for differential gene expression summary object, [calculateDifferentialGeneExpression()] for differential gene expression object
noninteractiveMeanDifferencePlot <- function(fit2, dT, ct) {
  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  fig <- plotMD(fit2, column=ct, status=dT[,ct], legend=FALSE, pch=20, cex=1)
  abline(h=0)
  return(fig)
}

#' A Function to Create a Histogram of the P values from Differential Gene Expression Analysis
#'
#' This function allows you to plot a histogram of the P values from differential gene expression analysis
#' @param fit2 An object containing the results of differntial gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A character string containing the adjustment to P-values
#' @keywords GEO
#' @import limma
#' @examples
#' #' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#'
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#'   if (i < halfNumberOfColumns) {
#'     group1 <- c(group1, name)
#'     i <- i +1
#'   } else {
#'     group2 <- c(group2, name)
#'     i <- i +1
#'   }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
#' ct <- 1
#'
#' # Non-Interactive Histogram
#' fig <- nonInteractiveHistogramPlot(fit2, adjustment)
#'
#' fig <- nonInteractiveHistogramPlot(fit2, adjustment)
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
nonInteractiveHistogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  fig <- hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
  return(fig)
}
