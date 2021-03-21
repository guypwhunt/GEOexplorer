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


interactiveMeanDifferencePlot <- function(fit2, adjustment, dT, ct) {
  fit2Df <- data.frame(Amean = fit2$Amean, sigma = fit2$sigma, coefficients = fit2$coefficients, expression = dT[,ct])
  fig <- plot_ly(data = fit2Df, x = ~Amean, y = ~Group1.Group2, color = ~Group1.Group2.1, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers', 
                 marker = list(size = 3))
  fig <- fig %>% layout(
    title = ('Group1-Group2'),
    xaxis = list(
      title = "Average log-expression"
    ),
    yaxis = list(
      title = "log-fold-change"
    ))
  fig}
  
interactiveVolcanoPlot <- function(fit2, adjustment, dT, ct) {
  fit2Df <- data.frame((0-log10(fit2$p.value)), fit2$coefficients, expression = dT[,ct])
  
  fig <- plot_ly(data = fit2Df, x = ~Group1.Group2.1, y = ~Group1.Group2, color = ~Group1.Group2.2, colors = c("blue", "black", "red"), type = 'scatter', mode = 'markers', 
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

interactiveQQPlot <- function(fit2) {
  t.good <- which(!is.na(fit2$F)) # filter out bad probes
  qqData <- qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic", plot.it = FALSE)
  
  fig <- plot_ly(x = qqData$x, y = qqData$y, type = 'scatter', mode = 'markers', 
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
  
  