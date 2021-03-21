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
  