library(plotly)
library(ggplot2)
library(limma)
library(scales)

interactiveHistogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  fig <- plot_ly(x = tT2$adj.P.Val, type = "histogram", nbinsx = 30)
  fig <- fig %>% layout(title = 'P-adj value distribution',
                        xaxis = list(title = 'P-adj'),
                        yaxis = list(title = 'Number of genes'))
  fig
}


interactiveMeanDifferencePlot <- function(fit2, dT, ct) {
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()
  
  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }

  fit2Df <- data.frame(fit2$Amean, fit2$coefficients, dT[,ct], fit2$genes[final_attributes_list])
  colnames(fit2Df) <- c("aMean", "coefficients", "regulation", final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <- "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"
  
  fig <- plot_ly(data = fit2Df, x = ~aMean, y = ~coefficients, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers', 
                 text = ~paste('ID: ', ID, '<br></br>', 'Average Log-Expression: ', aMean, '<br></br>', 'Log-Fold-Change: ', coefficients),
                 hoverinfo = text,
                 marker = list(size = 3))
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
  
interactiveVolcanoPlot <- function(fit2, dT, ct) {
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()
  
  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }
  fit2Df <- data.frame((0-log10(fit2$p.value)), fit2$coefficients, dT[,ct], fit2$genes[final_attributes_list])
  colnames(fit2Df) <- c("pValues", "coefficients", "regulation", final_attributes_list)
  fit2Df$regulation[fit2Df$regulation == "1"] <- "Upregulated"
  fit2Df$regulation[fit2Df$regulation == "0"] <- "Similar Expression"
  fit2Df$regulation[fit2Df$regulation == "-1"] <- "Downregulation"
  
  fig <- plot_ly(data = fit2Df, x = ~coefficients, y = ~pValues, color = ~regulation, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers', 
                 text = ~paste('ID: ', ID, '<br></br>', 'Log2 Fold Change: ', coefficients, '<br></br>', '-Log10(P-Value): ', pValues),
                 hoverinfo = text,
                 marker = list(size = 3))
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

interactiveQQPlot <- function(fit2, dT, ct) {
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqData <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic", plot.it = FALSE)
  
  attributes_list <- c('ID', 'Gene.symbol', 'Gene.title', 'Gene.ID')
  final_attributes_list <- c()
  
  for (attribute in attributes_list) {
    if (attribute %in% colnames(fit2$genes))
      final_attributes_list <- c(final_attributes_list, attribute)
  }

  qqData2 <- data.frame(qqData, dT[t.good,ct], fit2$genes[final_attributes_list][t.good,])
  
  colnames(qqData2) <- c("x", "y", "regulation", final_attributes_list)
  qqData2$regulation <- as.character(qqData2$regulation)
  qqData2$regulation[qqData2$regulation == "1"] <- "Upregulated"
  qqData2$regulation[qqData2$regulation == "0"] <- "Similar Expression"
  qqData2$regulation[qqData2$regulation == "-1"] <- "Downregulation"
  
  fig <- plot_ly()
  fig <- fig %>% add_trace( data = qqData2, x = ~x, y = ~y, type = 'scatter', mode = 'markers', color = ~regulation, colors = c("blue", "black", "red"), 
                            text = ~paste('ID: ', ID,'<br></br>', 'Theoretical Quantiles: ', x, '<br></br>', 'Sample Quantiles: ', y),
                            hoverinfo = text,
                            marker = list(size = 3))
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
  
  