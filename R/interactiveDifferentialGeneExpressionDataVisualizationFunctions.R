#' A Function to Create an Interactive Histogram of the
#' P values from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive
#' histogram of the P values from differential gene expression
#' analysis
#' @param fit2 An object containing the results of
#' differential gene expression analysis which can be
#' obtained from the calculateDifferentialGeneExpression()
#' function
#' @param adjustment A character string containing the
#' adjustment to P-values
#' @keywords GEO
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if
#' # necessary
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Interactive Histogram
#' fig <- interactiveHistogramPlot(fit2, adjustment)
#' fig
#'
#' @import plotly limma
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object
interactiveHistogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2,
                  adjust.method = adjustment,
                  sort.by = "B",
                  number = Inf)
  fig15 <- plot_ly(x = tT2$adj.P.Val,
                 type = "histogram",
                 nbinsx = 30)
  fig15 <- fig15 %>% layout(
    title = 'Adjusted P-value distribution',
    xaxis = list(title = 'Adjusted P-value'),
    yaxis = list(title = 'Number of genes')
  )
  
  try(fig15 <- toWebGL(fig15))
  
  return(fig15)
}

#' A Function to Create an Interactive Mean Difference
#' Plot of the log2 Fold Change Versus Average log2 Expression
#' Values from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive mean
#' difference plot of the log2 fold change versus average
#' log2 expression values from differential gene expression
#' analysis
#' @param fit2 An object containing the results of differential
#' gene expression analysis which can be obtained from the
#' calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is
#' unregulated, down regulated or has a similar level of
#' expression which can be obtained from the
#' calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select
#' from the dT object
#' @keywords GEO
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if
#' # necessary
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2,
#' adjustment, significanceLevelCutOff)
#'
#' # Plot Interactive Mean Difference of fit 2 data
#' ct <- 1
#' fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
#' fig
#'
#' fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object,
#' [calculateDifferentialGeneExpressionSummary()]
#' for summary differential gene expression object
interactiveMeanDifferencePlot <- function(fit2, dT, ct) {
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()

  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list,
                                 attribute)
  }

  if (is.null(fit2$genes[final_attributes_list]) == TRUE) {
    fit2Df <- data.frame(fit2$Amean, fit2$coefficients,
                         dT[, ct])
  } else {
    fit2Df <-
      data.frame(fit2$Amean, fit2$coefficients, dT[, ct],
                 fit2$genes[final_attributes_list])
  }

  colnames(fit2Df) <-
    c("aMean",
      "coefficients",
      "regulation",
      final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <-
    "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"

  fig16 <-
    plot_ly(
      data = fit2Df,
      x = ~ aMean,
      y = ~ coefficients,
      color = ~ regulation,
      colors = c("blue", "black", "red"),
      type = 'scatter',
      mode = 'markers',
      text = if ('ID' %in% final_attributes_list) {
        if ('Gene.symbol' %in% final_attributes_list) {
          if ('Gene.title' %in% final_attributes_list) {
            if ('Gene.ID' %in% final_attributes_list) {
              ~ paste(
                'ID: ',
                ID,
                '<br>',
                'Gene Symbol: ',
                Gene.symbol,
                '<br>',
                'Gene Title: ',
                Gene.title,
                '<br>',
                'Gene ID: ',
                Gene.ID,
                '<br>',
                'Average Log-Expression: ',
                aMean,
                '<br>',
                'Log-Fold-Change: ',
                coefficients
              )
            } else {
              ~ paste(
                'ID: ',
                ID,
                '<br>',
                'Gene Symbol: ',
                Gene.symbol,
                '<br>',
                'Gene Title: ',
                Gene.title,
                '<br>',
                'Average Log-Expression: ',
                aMean,
                '<br>',
                'Log-Fold-Change: ',
                coefficients
              )
            }
          } else {
            ~ paste(
              'ID: ',
              ID,
              '<br>',
              'Gene Symbol: ',
              Gene.symbol,
              '<br>',
              '<br>',
              'Average Log-Expression: ',
              aMean,
              '<br>',
              'Log-Fold-Change: ',
              coefficients
            )
          }
        } else {
          ~ paste(
            'ID: ',
            ID,
            '<br>',
            'Average Log-Expression: ',
            aMean,
            '<br>',
            'Log-Fold-Change: ',
            coefficients
          )
        }
      } else {
        ~ paste('Average Log-Expression: ',
                aMean,
                '<br>',
                'Log-Fold-Change: ',
                coefficients)
      }
      ,
      hoverinfo = text,
      marker = list(size = 3)
    )

  fig16 <- fig16 %>% layout(
    title = ('Group1-Group2'),
    xaxis = list(title = "Average log-expression"),
    yaxis = list(title = "log-fold-change")
  )
  
  return(fig16)
}

#' A Function to Create an Interactive Volcano Plot of
#' the Statistical Significance (-log10 P Value) Versus
#' Magnitude of Change (log2 Fold Change) from Differential
#' Gene Expression Analysis
#'
#' This function allows you to plot an interactive volcano plot
#' of the statistical significance (-log10 P value) versus
#' magnitude of change (log2 fold change) from differential
#' gene expression analysis
#' @param fit2 An object containing the results of differential
#' gene expression analysis which can be obtained from the
#' calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is
#' unregulated, down regulated or has a similar level of
#' expression which can be obtained from the
#' calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to
#' select from the dT object
#' @keywords GEO
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2,
#' adjustment, significanceLevelCutOff)
#' ct <- 1
#'
#' # Interactive volcano plot (log P-value vs log fold change)
#' fig <- interactiveVolcanoPlot(fit2, dT, ct)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object,
#' [calculateDifferentialGeneExpressionSummary()]
#' for summary differential gene expression object
interactiveVolcanoPlot <- function(fit2, dT, ct) {
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title',
                       'Gene.ID')
  final_attributes_list <- c()

  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list,
                                 attribute)
  }

  if (is.null(fit2$genes[final_attributes_list]) == TRUE) {
    fit2Df <-
      data.frame((0 - log10(fit2$p.value)), fit2$coefficients,
                 dT[, ct])
  } else {
    fit2Df <-
      data.frame((0 - log10(fit2$p.value)), fit2$coefficients,
                 dT[, ct], fit2$genes[final_attributes_list])
  }

  colnames(fit2Df) <-
    c("pValues",
      "coefficients",
      "regulation",
      final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <-
    "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"

  fig17 <-
    plot_ly(
      data = fit2Df,
      x = ~ coefficients,
      y = ~ pValues,
      color = ~ regulation,
      colors = c("blue", "black", "red"),
      type = 'scatter',
      mode = 'markers',
      text =
        if ('ID' %in% final_attributes_list) {
          if ('Gene.symbol' %in% final_attributes_list) {
            if ('Gene.title' %in% final_attributes_list) {
              if ('Gene.ID' %in% final_attributes_list) {
                ~ paste(
                  'ID: ',
                  ID,
                  '<br>',
                  'Gene Symbol: ',
                  Gene.symbol,
                  '<br>',
                  'Gene Title: ',
                  Gene.title,
                  '<br>',
                  'Gene ID: ',
                  Gene.ID,
                  '<br>',
                  'Log2 Fold Change: ',
                  coefficients,
                  '<br>',
                  '-Log10(P-Value): ',
                  pValues
                )
              } else {
                ~ paste(
                  'ID: ',
                  ID,
                  '<br>',
                  'Gene Symbol: ',
                  Gene.symbol,
                  '<br>',
                  'Gene Title: ',
                  Gene.title,
                  '<br>',
                  'Log2 Fold Change: ',
                  coefficients,
                  '<br>',
                  '-Log10(P-Value): ',
                  pValues
                )
              }
            } else {
              ~ paste(
                'ID: ',
                ID,
                '<br>',
                'Gene Symbol: ',
                Gene.symbol,
                '<br>',
                '<br>',
                'Log2 Fold Change: ',
                coefficients,
                '<br>',
                '-Log10(P-Value): ',
                pValues
              )
            }
          } else {
            ~ paste(
              'ID: ',
              ID,
              '<br>',
              'Log2 Fold Change: ',
              coefficients,
              '<br>',
              '-Log10(P-Value): ',
              pValues
            )
          }
        } else {
          ~ paste('Log2 Fold Change: ',
                  coefficients,
                  '<br>',
                  '-Log10(P-Value): ',
                  pValues)
        }
      ,
      hoverinfo = text,
      marker = list(size = 3)
    )

  fig17 <- fig17 %>% layout(
    title = ('Group1-Group2'),
    xaxis = list(title = "Log2 Fold Change"),
    yaxis = list(title = "-log10(P-value)")
  )

  return(fig17)
}

#' A Function to Create an Interactive QQ Plot of the
#' Quantiles of a Data Sample Against the Theoretical
#' Quantiles of a Student's T Distribution from Differential
#' Gene Expression Analysis
#'
#' This function allows you to plot an interactive QQ plot
#' of the quantiles of a data sample against the theoretical
#' quantiles of a Student's t distribution from differential
#' gene expression analysis
#' @param fit2 An object containing the results of differential
#' gene expression analysis which can be obtained from the
#' calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is
#' unregulated, down regulated or has a similar level of
#' expression which can be obtained from the
#' calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from
#' the dT object
#' @keywords GEO
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2,
#' adjustment, significanceLevelCutOff)
#'
#' # Interactive Q-Q plot
#' ct <- 1
#' fig <- interactiveQQPlot(fit2, dT, ct)
#' fig
#'
#' @import plotly limma
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object,
#' [calculateDifferentialGeneExpressionSummary()]
#' for summary differential gene expression object
interactiveQQPlot <- function(fit2, dT, ct) {
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqData <-
    qqt(fit2$t[t.good],
        fit2$df.total[t.good],
        main = "Moderated t statistic",
        plot.it = FALSE)
  
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title',
                       'Gene.ID')
  final_attributes_list <- c()
  
  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list,
                                 attribute)
  }
  
  genes <- tryCatch({
    fit2$genes[,final_attributes_list][t.good,]
  }, error=function(cond) {
    fit2$genes[,final_attributes_list][t.good]
  }
  )
  
  if (is.null(genes)) {
    qqData2 <- data.frame(qqData, dT[t.good, ct])
  } else {
    qqData2 <-
      data.frame(qqData, dT[t.good, ct],
                 genes)
  }
  
  colnames(qqData2) <-
    c("x", "y", "regulation", final_attributes_list)
  qqData2$regulation <- as.character(qqData2$regulation)
  qqData2$regulation[qqData2$regulation == "1"] <- "Upregulated"
  qqData2$regulation[qqData2$regulation == "0"] <-
    "Similar Expression"
  qqData2$regulation[qqData2$regulation == "-1"] <- "Downregulation"
  
  fig18 <- plot_ly()
  fig18 <-
    fig18 %>% add_trace(
      data = qqData2,
      x = ~ x,
      y = ~ y,
      type = 'scatter',
      mode = 'markers',
      color = ~ regulation,
      colors = c("blue", "black", "red"),
      text =
        if ('ID' %in% final_attributes_list) {
          if ('Gene.symbol' %in% final_attributes_list) {
            if ('Gene.title' %in% final_attributes_list) {
              if ('Gene.ID' %in% final_attributes_list) {
                ~ paste(
                  'ID: ',
                  ID,
                  '<br>',
                  'Gene Symbol: ',
                  Gene.symbol,
                  '<br>',
                  'Gene Title: ',
                  Gene.title,
                  '<br>',
                  'Gene ID: ',
                  Gene.ID,
                  '<br>',
                  'Theoretical Quantiles: ',
                  x,
                  '<br>',
                  'Sample Quantiles: ',
                  y
                )
              } else {
                ~ paste(
                  'ID: ',
                  ID,
                  '<br>',
                  'Gene Symbol: ',
                  Gene.symbol,
                  '<br>',
                  'Gene Title: ',
                  Gene.title,
                  '<br>',
                  'Theoretical Quantiles: ',
                  x,
                  '<br>',
                  'Sample Quantiles: ',
                  y
                )
              }
            } else {
              ~ paste(
                'ID: ',
                ID,
                '<br>',
                'Gene Symbol: ',
                Gene.symbol,
                '<br>',
                'Theoretical Quantiles: ',
                x,
                '<br>',
                'Sample Quantiles: ',
                y
              )
            }
          } else {
            ~ paste(
              'ID: ',
              ID,
              '<br>',
              'Theoretical Quantiles: ',
              x,
              '<br>',
              'Sample Quantiles: ',
              y
            )
          }
        } else {
          ~ paste('Theoretical Quantiles: ',
                  x,
                  '<br>',
                  'Sample Quantiles: ',
                  y)
        }
      ,
      hoverinfo = text,
      marker = list(size = 3)
    )
  fig18 <- fig18 %>% layout(
    title = ('Moderated t statistic'),
    xaxis = list(title = "Theoretical Quantiles"),
    yaxis = list(title = "Sample Quantiles")
  )

  return(fig18)
}


#' A Function to Create an Heatmap Plot of the
#' Top Differentially expressed genes for each experimental condition
#'
#' This function allows you to Create an Heatmap Plot of the
#' Top Differentially expressed genes for each experimental condition
#' @import heatmaply
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object,
#' [calculateDifferentialGeneExpressionSummary()]
#' for summary differential gene expression object
interactiveDGEHeatMapPlot <- function(ex, limmaPrecisionWeights,
                                      numberOfGenes, tT) {
  # Select the genes of interest
  i <- row.names(tT[seq_len(numberOfGenes),])

  # Create the heatmap
  if (limmaPrecisionWeights == "Yes") {
    heatmapFig <- heatmaply(ex$E[i, ])
  } else if (limmaPrecisionWeights == "No") {
    heatmapFig <- heatmaply(ex[i, ])
  }
  
  return(heatmapFig)

}
