#' A Function to Create an Interactive Histogram of the P values from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive histogram of the P values from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A character string containing the adjustment to P-values
#' @keywords GEO
#' @export
#' @examples fig <- interactiveHistogramPlot(fit2, adjustment)
#' @import plotly limma
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
interactiveHistogramPlot <- function(fit2, adjustment) {
  library(plotly)
  library(limma)
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  fig <- plot_ly(x = tT2$adj.P.Val, type = "histogram", nbinsx = 30)
  fig <- fig %>% layout(title = 'Adjusted P-value distribution',
                        xaxis = list(title = 'Adjusted P-value'),
                        yaxis = list(title = 'Number of genes'))
  fig
}

#' A Function to Create an Interactive Mean Difference Plot of the log2 Fold Change Versus Average log2 Expression Values from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive mean difference plot of the log2 fold change versus average log2 expression values from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @export
#' @examples fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
#' @import plotly
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object, [calculateDifferentialGeneExpressionSummary()] for summary differential gene expression object
interactiveMeanDifferencePlot <- function(fit2, dT, ct) {
  library(plotly)
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()

  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }

  if (is.null(fit2$genes[final_attributes_list]) == TRUE) {
    fit2Df <- data.frame(fit2$Amean, fit2$coefficients, dT[,ct])
  } else {
    fit2Df <- data.frame(fit2$Amean, fit2$coefficients, dT[,ct], fit2$genes[final_attributes_list])
  }

  colnames(fit2Df) <- c("aMean", "coefficients", "regulation", final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <- "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"

  if('ID' %in% final_attributes_list){
    if('Gene.symbol' %in% final_attributes_list){
      if('Gene.title' %in% final_attributes_list){
        if('Gene.ID' %in% final_attributes_list){
          fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                         hoverinfo = text,
                         marker = list(size = 3))
        } else {
          fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                         hoverinfo = text,
                         marker = list(size = 3))
        }
      } else {
        fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', '<br></br>', 'Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                       hoverinfo = text,
                       marker = list(size = 3))
      }
    } else{
      fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                     text = ~paste('ID: ', ID, '<br></br>', 'Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                     hoverinfo = text,
                     marker = list(size = 3))
    }
  } else{
    fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                   text = ~paste('Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                   hoverinfo = text,
                   marker = list(size = 3))
  }

  fig <- fig %>% layout(
    title = ('Group1-Group2'),
    xaxis = list(
      title = "Average log-expression"
    ),
    yaxis = list(
      title = "log-fold-change"
    ))
  fig
  }

#' A Function to Create an Interactive Volcano Plot of the Statistical Significance (-log10 P Value) Versus Magnitude of Change (log2 Fold Change) from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive volcano plot of the statistical significance (-log10 P value) versus magnitude of change (log2 fold change) from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @export
#' @examples fig <- interactiveVolcanoPlot(fit2, dT, ct)
#' @import plotly
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object, [calculateDifferentialGeneExpressionSummary()] for summary differential gene expression object
interactiveVolcanoPlot <- function(fit2, dT, ct) {
  library(plotly)
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()

  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }

  if (is.null(fit2$genes[final_attributes_list]) == TRUE) {
    fit2Df <- data.frame((0-log10(fit2$p.value)), fit2$coefficients, dT[,ct])
  } else {
    fit2Df <- data.frame((0-log10(fit2$p.value)), fit2$coefficients, dT[,ct], fit2$genes[final_attributes_list])
  }

  colnames(fit2Df) <- c("pValues", "coefficients", "regulation", final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <- "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"

  if('ID' %in% final_attributes_list){
    if('Gene.symbol' %in% final_attributes_list){
      if('Gene.title' %in% final_attributes_list){
        if('Gene.ID' %in% final_attributes_list){
          fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                         hoverinfo = text,
                         marker = list(size = 3))
        } else {
          fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                         hoverinfo = text,
                         marker = list(size = 3))
        }
      } else {
        fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', '<br></br>', 'Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                       hoverinfo = text,
                       marker = list(size = 3))
      }
    } else{
      fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                     text = ~paste('ID: ', ID, '<br></br>', 'Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                     hoverinfo = text,
                     marker = list(size = 3))
    }
  } else{
    fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers',
                   text = ~paste('Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                   hoverinfo = text,
                   marker = list(size = 3))
  }
  fig <- fig %>% layout(
    title = ('Group1-Group2'),
    xaxis = list(
      title = "Log2 Fold Change"
    ),
    yaxis = list(
      title = "-log10(P-value)"
    ))
  fig
}

#' A Function to Create an Interactive QQ Plot of the Quantiles of a Data Sample Against the Theoretical Quantiles of a Student's T Distribution from Differential Gene Expression Analysis
#'
#' This function allows you to plot an interactive QQ plot of the quantiles of a data sample against the theoretical quantiles of a Student's t distribution from differential gene expression analysis
#' @param fit2 An object containing the results of differential gene expression analysis which can be obtained from the calculateDifferentialGeneExpression() function
#' @param dT An object that summarises if each gene is unregulated, down regulated or has a similar level of expression which can be obtained from the calculateDifferentialGeneExpressionSummary() function
#' @param ct A integer indicating the column to select from the dT object
#' @keywords GEO
#' @export
#' @examples fig <- interactiveQQPlot(fit2, dT, ct)
#' @import plotly limma
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object, [calculateDifferentialGeneExpressionSummary()] for summary differential gene expression object
interactiveQQPlot <- function(fit2, dT, ct) {
  library(plotly)
  library(limma)
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqData <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic", plot.it = FALSE)

  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()

  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }

  if (is.null(fit2$genes[final_attributes_list][t.good,]) == TRUE) {
    qqData2 <- data.frame(qqData, dT[t.good,ct])
  } else {
    qqData2 <- data.frame(qqData, dT[t.good,ct], fit2$genes[final_attributes_list][t.good,])
  }


  colnames(qqData2) <- c("x", "y", "regulation", final_attributes_list)
  qqData2$regulation <- as.character(qqData2$regulation)
  qqData2$regulation[qqData2$regulation == "1"] <- "Upregulated"
  qqData2$regulation[qqData2$regulation == "0"] <- "Similar Expression"
  qqData2$regulation[qqData2$regulation == "-1"] <- "Downregulation"

  fig <- plot_ly()
  if('ID' %in% final_attributes_list){
    if('Gene.symbol' %in% final_attributes_list){
      if('Gene.title' %in% final_attributes_list){
        if('Gene.ID' %in% final_attributes_list){
          fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"),
                                    text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                                    hoverinfo = text,
                                    marker = list(size = 3))
        } else {
          fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"),
                                    text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                                    hoverinfo = text,
                                    marker = list(size = 3))
        }
      } else {
        fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"),
                                  text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                                  hoverinfo = text,
                                  marker = list(size = 3))
      }
    } else{
      fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"),
                                text = ~paste('ID: ', ID, '<br></br>', 'Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                                hoverinfo = text,
                                marker = list(size = 3))
    }
  } else{
    fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"),
                              text = ~paste('Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                              hoverinfo = text,
                              marker = list(size = 3))
  }
  fig <- fig %>% layout(
    title = ('Moderated t statistic'),
    xaxis = list(
      title = "Theoretical Quantiles"
    ),
    yaxis = list(
      title = "Sample Quantiles"
    ))
  fig
}

