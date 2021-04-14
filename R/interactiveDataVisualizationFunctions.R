#' @export
#' @import plotly ggplot2 limma
interactiveBoxAndWhiskerPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
  library(ggplot2)
  library(limma)
  ex <- as.data.frame(ex)
  fig <- plot_ly(type = "box", quartilemethod="linear")
  i = 1
  for(col in names(ex)) {
    fig <- fig %>% add_trace(x = names(ex)[i], y = ex[,i], quartilemethod="linear", name=names(ex)[i])
    i <- i+1
  }
  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode, "/"), platform)))
  fig
}

#' @export
#' @import plotly ggplot2 limma
interactiveDesnityPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
  library(ggplot2)
  library(limma)
  ex <- as.data.frame(ex)
  fig <- plot_ly(type = 'scatter', mode = 'lines', name = (paste(paste(geoAccessionCode,platform),'value distribution')))
  i <- 1
  for(col in names(ex)) {
    density <- density(ex[,i])
    fig <- fig %>% add_trace(x = density$x, y = density$y, name = col)
    i <- i+1
  }

  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode,platform),'value distribution')),
                        xaxis = list(title = 'Intensity'),
                        yaxis = list(title = 'Density'))
  fig
}

#' @export
#' @import plotly ggplot2 limma
interactiveThreeDDesnityPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
  library(ggplot2)
  library(limma)
  ex <- as.data.frame(ex)
  fig <- plot_ly(type = 'scatter3d', mode = 'lines', name = (paste(paste(geoAccessionCode,platform),'value distribution')))
  i <- 1
  for(col in names(ex)) {
    density <- density(ex[,i])
    fig <- fig %>% add_trace(x = density$x, y = i, z = density$y, name = col)
    i <- i+1
  }

  fig <- fig %>% layout(
    title = (paste(paste(geoAccessionCode,platform),'value distribution')),
    scene = list(
      xaxis = list(title = "Intensity"),
      yaxis = list(title = ""),
      zaxis = list(title = "Density")
    ))
  fig
}

#' @export
#' @import plotly ggplot2 limma
interactiveUmapPlot <- function(ex, knn, geoAccessionCode) {
  library(plotly)
  library(ggplot2)
  library(limma)
  ex <- ex[!duplicated(ex), ]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
  i <- 1
  fig <- plot_ly(type = 'scatter', mode = 'markers')
  for(row in row.names(ump$layout)){
    fig <- fig %>% add_trace(x = ump$layout[i,][1], y = ump$layout[i,][2], name = row)
    i <- i+1
  }
  fig <- fig %>% layout(
    title = (paste(geoAccessionCode, paste('UMAP plot, number of nearest neighbors used =',knn))))
  fig
}

#' @export
#' @import plotly ggplot2 limma stringr
interactiveMeanVariancePlot <- function(ex, geoAccessionCode, gset) {
  library(plotly)
  library(ggplot2)
  library(limma)
  library(stringr)
  ex <- lmFit(ex)
  ex <- as.data.frame(ex)
  ex["ID"] <- rownames(ex)
  geneData <- gset@featureData@data
  geneData <- as.data.frame(geneData)
  combineData <- merge(ex, geneData, by = "ID")
  colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")
  combineData %>% filter(ID %in% c(rownames(ex)))

  if('ID' %in% colnames(combineData)){
    if('Gene.symbol' %in% colnames(combineData)){
      if('Gene.title' %in% colnames(combineData)){
        if('Gene.ID' %in% colnames(combineData)){
          fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                         hoverinfo = text,
                         mode = 'markers', marker = list(
                           color = 'rgb(17, 157, 255)',
                           size = 3,
                           line = list(
                             color = 'rgb(0, 0, 0)',
                             width = 1
                           )))
        } else {
          fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                         hoverinfo = text,
                         mode = 'markers', marker = list(
                           color = 'rgb(17, 157, 255)',
                           size = 3,
                           line = list(
                             color = 'rgb(0, 0, 0)',
                             width = 1
                           )))
        }
      } else {
        fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                       hoverinfo = text,
                       mode = 'markers', marker = list(
                         color = 'rgb(17, 157, 255)',
                         size = 3,
                         line = list(
                           color = 'rgb(0, 0, 0)',
                           width = 1
                         )))
      }
    } else{
      fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                     text = ~paste('ID: ', ID, '<br></br>', 'Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                     hoverinfo = text,
                     mode = 'markers', marker = list(
                       color = 'rgb(17, 157, 255)',
                       size = 3,
                       line = list(
                         color = 'rgb(0, 0, 0)',
                         width = 1
                       )))
    }
  } else{
    fig <- plot_ly(data = combineData, x = ~Amean, y = ~sigma, type = 'scatter',
                   text = ~paste('Amean: ', Amean, '<br></br>', 'Sigma: ', sigma),
                   hoverinfo = text,
                   mode = 'markers', marker = list(
                     color = 'rgb(17, 157, 255)',
                     size = 3,
                     line = list(
                       color = 'rgb(0, 0, 0)',
                       width = 1
                     )))
  }
  fig <- fig %>% layout(
    title = (paste('Mean variance trend, ',geoAccessionCode)))
  fig
}

#' @export
#' @import plotly ggplot2 limma
interactivePcaScreePlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
  library(ggplot2)
  library(limma)
  columnNames <- colnames(pcaData$x)
  proportionOfVariance <- pcaData$sdev^2/sum(pcaData$sdev^2)
  pcaDataFrame <- data.frame(columnNames, proportionOfVariance)

  fig <- plot_ly(data = pcaDataFrame, x = ~columnNames, y = ~proportionOfVariance, type = "bar") %>%
    layout(
      title = paste(geoAccessionCode, "Scree Plot"),
      xaxis = list(title = "Principal Components/Dimensions",
                   categoryorder = "array",
                   categoryarray = ~columnNames),
      yaxis = list(title = "Percentage of Explained Variances",
                   tickformat = "%")
    )

  fig
}

#' @export
#' @import plotly ggplot2 limma
interactivePrincompPcaScreePlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
  library(ggplot2)
  library(limma)
  columnNames <-   colnames(pcaData$loadings)
  proportionOfVariance <- pcaData$sdev^2/sum(pcaData$sdev^2)
  pcaDataFrame <- data.frame(columnNames, proportionOfVariance)
  fig <- plot_ly(data = pcaDataFrame, x = ~columnNames, y = ~proportionOfVariance, type = "bar") %>%
    layout(
      title = paste(geoAccessionCode, "Scree Plot"),
      xaxis = list(title = "Principal Components/Dimensions",
                   categoryorder = "array",
                   categoryarray = ~columnNames),
      yaxis = list(title = "Percentage of Explained Variances",
                   tickformat = "%")
    )

  fig
}

#' @export
#' @import plotly ggplot2 limma stringr
interactivePrincompPcaIndividualsPlot <- function(pcaData, geoAccessionCode, gset) {
  library(plotly)
  library(ggplot2)
  library(limma)
  library(stringr)
  pcaDf <- data.frame(pcaData$scores)
  pcaDf <- transform(pcaDf)
  pcaDf["ID"] <- rownames(pcaDf)
  geneData <- gset@featureData@data
  geneData <- as.data.frame(geneData)
  combineData <- merge(pcaDf, geneData, by = "ID")
  combineData %>% filter(ID %in% c(rownames(pcaDf)))
  colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")

  individualsStats <- get_pca_ind(pcaData)
  eigenValue <- get_eigenvalue(pcaData)

  if('ID' %in% colnames(combineData)){
    if('Gene.symbol' %in% colnames(combineData)){
      if('Gene.title' %in% colnames(combineData)){
        if('Gene.ID' %in% colnames(combineData)){
          fig <- plot_ly(combineData,x=~Comp.1,y=~Comp.2, mode="markers", type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Dimension 1: ', Comp.1, '<br></br>', 'Dimension 2: ', Comp.2),
                         hoverinfo = text,
                         marker = list(
                           color = ~individualsStats$cos2[,1],
                           size = 3
                         ))
        } else {
          fig <- plot_ly(combineData,x=~Comp.1,y=~Comp.2, mode="markers", type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Dimension 1: ', Comp.1, '<br></br>', 'Dimension 2: ', Comp.2),
                         hoverinfo = text,
                         marker = list(
                           color = ~individualsStats$cos2[,1],
                           size = 3
                         ))
        }
      } else {
        fig <- plot_ly(combineData,x=~Comp.1,y=~Comp.2, mode="markers", type = 'scatter',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Dimension 1: ', Comp.1, '<br></br>', 'Dimension 2: ', Comp.2),
                       hoverinfo = text,
                       marker = list(
                         color = ~individualsStats$cos2[,1],
                         size = 3
                       ))
      }
    } else{
      fig <- plot_ly(combineData,x=~Comp.1,y=~Comp.2, mode="markers", type = 'scatter',
                     text = ~paste('ID: ', ID, '<br></br>', 'Dimension 1: ', Comp.1, '<br></br>', 'Dimension 2: ', Comp.2),
                     hoverinfo = text,
                     marker = list(
                       color = ~individualsStats$cos2[,1],
                       size = 3
                     ))
    }
  } else{
    fig <- plot_ly(combineData,x=~Comp.1,y=~Comp.2, mode="markers", type = 'scatter',
                   text = ~paste('Dimension 1: ', Comp.1, '<br></br>', 'Dimension 2: ', Comp.2),
                   hoverinfo = text,
                   marker = list(
                     color = ~individualsStats$cos2[,1],
                     size = 3
                   ))
  }

  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Individuals Plot"),
                xaxis=list(title=paste("Comp.1", label_percent(accuracy=0.1)(eigenValue[1,2]/100))),
                yaxis=list(title=paste("Comp.2", label_percent(accuracy=0.1)(eigenValue[2,2]/100))))
  fig

}

#' @export
#' @import plotly ggplot2 limma
interactivePrincompPcaVariablesPlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
  library(ggplot2)
  library(limma)
  variableStats <- get_pca_var(pcaData)
  eigenValue <- get_eigenvalue(pcaData)
  pcaData <- as.data.frame(unclass(pcaData$loadings))

  fig <- plot_ly(pcaData,x=~Comp.1,y=~Comp.2,text=rownames(pcaData), mode="markers", type = 'scatter'
                 ,marker=list(size=10, color = ~variableStats$contrib[,1]), name = rownames(pcaData))

  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Variables Plot"),
                xaxis=list(title=paste("Comp.1", label_percent(accuracy=0.1)(eigenValue[1,2]/100))),
                yaxis=list(title=paste("Comp.2", label_percent(accuracy=0.1)(eigenValue[2,2]/100))))
  fig
}

#' @export
#' @import plotly ggplot2 limma scales pheatmap heatmaply stringr
interactiveHeatMapPlot <- function(ex) {
  library(plotly)
  library(ggplot2)
  library(limma)
  library(scales)
  library(pheatmap)
  library(heatmaply)
  library(stringr)
  corMatrix <- cor(ex,use="c")
  df <- data.frame(corMatrix[1,])
  df <- df[-c(1)]
  i = 1
  while(i <= length(colnames(corMatrix))){
    df[i] <- data.frame(corMatrix[,i])
    i <- i + 1
  }
  colnames(df) <- colnames(corMatrix)
  fig <- heatmaply(df)
  fig
}
