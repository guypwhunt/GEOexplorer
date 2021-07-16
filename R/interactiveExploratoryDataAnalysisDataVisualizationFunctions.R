#' A Function to Create an Interactive Box and Whisker Plot from an Expression Object
#'
#' This function allows you to plot expression data into an interactive Box and Whisker Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveBoxAndWhiskerPlot(expressionData, "GSE18380", "GPL4694")
#' @import plotly
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
interactiveBoxAndWhiskerPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
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

#' A Function to Create an Interactive Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into an interactive Density Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveDensityPlot(expressionData, "GSE18380", "GPL4694")
#' @import plotly
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
interactiveDensityPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
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

#' A Function to Create an Interactive Three Dimensional Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into an interactive Three Dimensional Density Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform A character string representing the study's platform
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveThreeDDensityPlot(expressionData, "GSE18380", "GPL4694")
#' @import plotly
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
interactiveThreeDDensityPlot <- function(ex, geoAccessionCode, platform) {
  library(plotly)
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

#' A Function to Create an Interactive UMAP Plot from an Expression Object
#'
#' This function allows you to plot expression data into an interactive UMAP Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param knn A integer representing the number of k nearest neighbor's to use
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveUmapPlot(expressionData, 2, "GSE18380")
#' @import plotly umap
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
interactiveUmapPlot <- function(ex, knn, geoAccessionCode) {
  library(plotly)
  library(umap)
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

#' A Function to Create an Interactive Mean Variance Plot from an Expression Object
#'
#' This function allows you to plot expression data into an interactive Mean Variance Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveMeanVariancePlot(expressionData, "GSE18380", gset)
#' @import plotly limma
#' @importFrom stringr str_replace_all
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object, [extractPlatformGset()] for GEO object
interactiveMeanVariancePlot <- function(ex, geoAccessionCode, gset) {
  library(plotly)
  library(limma)
  library(stringr)
  exData <- lmFit(ex)
  exData <- as.data.frame(exData)
  exData["ID"] <- rownames(ex)
  geneData <- gsetData@featureData@data
  # Error handling to catch gset without featureData
  if(ncol(geneData) > 0){
    geneData <- as.data.frame(geneData)
    combineData <- merge(exData, geneData, by = "ID")
    colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")
    combineData %>% filter(ID %in% c(rownames(exData)))
  } else{
    combineData <- exData
  }

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

#' A Function to Create an Interactive Histogram of the Prcomp Principle Components from the PCA outputs of an Expression Object
#'
#' This function allows you to plot Prcomp PCA expression results into an interactive Histogram of the Principle Components
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePca() function
#' @keywords GEO
#' @export
#' @examples fig <- interactivePrcompPcaScreePlot(pcaData, "GSE18380")
#' @import plotly
#' @author Guy Hunt
#' @seealso [calculatePca()] for PCA expression object
interactivePrcompPcaScreePlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
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

#' A Function to Create an Interactive Histogram of the princomp Principle Components from the PCA outputs of an Expression Object
#'
#' This function allows you to plot princomp PCA expression results into an interactive Histogram of the Principle Components
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @examples fig <- interactivePrincompPcaScreePlot(pcaData, "GSE18380")
#' @import plotly
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp PCA expression object
interactivePrincompPcaScreePlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
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

#' A Function to Create an Interactive Scatter Plot of the princomp Principle Components Analysis of each of the Genes in an Expression Object
#'
#' This function allows you to plot princomp PCA expression results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePrincompPca() function
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @examples fig <- interactivePrincompPcaIndividualsPlot(pcaData, "GSE18380", gset)
#' @import plotly
#' @importFrom stringr str_replace_all
#' @importFrom factoextra get_pca_ind get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp PCA expression object, [extractPlatformGset()] for GEO object
interactivePrincompPcaIndividualsPlot <- function(pcaData, geoAccessionCode, gset) {
  library(plotly)
  library(stringr)
  library(factoextra)
  library(scales)

  pcaDf <- data.frame(pcaData$scores)
  pcaDf <- transform(pcaDf)
  pcaDf["ID"] <- rownames(pcaDf)
  geneData <- gset@featureData@data
  # Error handling for gset without featureData@data
  if(ncol(geneData) > 0){
  geneData <- as.data.frame(geneData)
  combineData <- merge(pcaDf, geneData, by = "ID")
  combineData %>% filter(ID %in% c(rownames(pcaDf)))
  colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")
  } else {
    combineData <- pcaDf
  }


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

#' A Function to Create an Interactive Scatter Plot of the princomp Principle Components Analysis of each of the Samples in an Expression Object
#'
#' This function allows you to plot princomp PCA expression results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @examples fig <- interactivePrincompPcaVariablesPlot(pcaData, "GSE18380")
#' @import plotly
#' @importFrom factoextra get_pca_var get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @seealso [calculatePrincompPca()] for Princomp PCA expression object, [extractPlatformGset()] for GEO object
interactivePrincompPcaVariablesPlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
  library(factoextra)
  library(scales)

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

#' A Function to Create an Interactive Heat Map of the Correlations between Samples
#'
#' This function allows you to plot an interactive heat map of the correlations between samples in an expression object
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveHeatMapPlot(expressionData)
#' @import plotly
#' @importFrom heatmaply heatmaply
#' @author Guy Hunt
interactiveHeatMapPlot <- function(ex) {
  library(plotly)
  library(heatmaply)
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

#' A Function to Create an Interactive Scatter Plot of the prcomp Principle Components Analysis of each of the Genes in an Expression Object
#'
#' This function allows you to plot prcomp PCA expression results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePrincompPca() function
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @examples fig <- interactiveprcompPcaIndividualsPlot(pcaData, "GSE18380", gset)
#' @import plotly
#' @importFrom stringr str_replace_all
#' @importFrom factoextra get_pca_ind get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @seealso [calculatePrcompPca()] for Princomp PCA expression object, [extractPlatformGset()] for GEO object
interactivePrcompPcaIndividualsPlot <- function(pcaData, geoAccessionCode, gset) {
  library(plotly)
  library(factoextra)
  library(stringr)
  library(scales)
  pcaDf <- data.frame(pcaData$x)
  pcaDf <- transform(pcaDf)
  pcaDf["ID"] <- rownames(pcaDf)
  geneData <- gset@featureData@data
  # Error handling for gset without featureData@data
  if(ncol(geneData) > 0){
    geneData <- as.data.frame(geneData)
    combineData <- merge(pcaDf, geneData, by = "ID")
    combineData %>% filter(ID %in% c(rownames(pcaDf)))
    colnames(combineData) <- str_replace_all(colnames(combineData), " ", ".")
  } else {
    combineData <- pcaDf
  }


  individualsStats <- get_pca_ind(pcaData)
  eigenValue <- get_eigenvalue(pcaData)


  if('ID' %in% colnames(combineData)){
    if('Gene.symbol' %in% colnames(combineData)){
      if('Gene.title' %in% colnames(combineData)){
        if('Gene.ID' %in% colnames(combineData)){
          fig <- plot_ly(combineData,x=~PC1,y=~PC2, mode="markers", type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Gene ID: ', Gene.ID, '<br></br>', 'Dimension 1: ', PC1, '<br></br>', 'Dimension 2: ', PC2),
                         hoverinfo = text,
                         marker = list(
                           color = ~individualsStats$cos2[,1],
                           size = 3
                         ))
        } else {
          fig <- plot_ly(combineData,x=~PC1,y=~PC2, mode="markers", type = 'scatter',
                         text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Gene Title: ', Gene.title, '<br></br>', 'Dimension 1: ', PC1, '<br></br>', 'Dimension 2: ', PC2),
                         hoverinfo = text,
                         marker = list(
                           color = ~individualsStats$cos2[,1],
                           size = 3
                         ))
        }
      } else {
        fig <- plot_ly(combineData,x=~PC1,y=~PC2, mode="markers", type = 'scatter',
                       text = ~paste('ID: ', ID, '<br></br>', 'Gene Symbol: ', Gene.symbol, '<br></br>', 'Dimension 1: ', PC1, '<br></br>', 'Dimension 2: ', PC2),
                       hoverinfo = text,
                       marker = list(
                         color = ~individualsStats$cos2[,1],
                         size = 3
                       ))
      }
    } else{
      fig <- plot_ly(combineData,x=~PC1,y=~PC2, mode="markers", type = 'scatter',
                     text = ~paste('ID: ', ID, '<br></br>', 'Dimension 1: ', PC1, '<br></br>', 'Dimension 2: ', PC2),
                     hoverinfo = text,
                     marker = list(
                       color = ~individualsStats$cos2[,1],
                       size = 3
                     ))
    }
  } else{
    fig <- plot_ly(combineData,x=~PC1,y=~PC2, mode="markers", type = 'scatter',
                   text = ~paste('Dimension 1: ', PC1, '<br></br>', 'Dimension 2: ', PC2),
                   hoverinfo = text,
                   marker = list(
                     color = ~individualsStats$cos2[,1],
                     size = 3
                   ))
  }


  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Individuals Plot"),
                xaxis=list(title=paste("PC1", label_percent(accuracy=0.1)(eigenValue[1,2]/100))),
                yaxis=list(title=paste("PC2", label_percent(accuracy=0.1)(eigenValue[2,2]/100))))
  fig
}

#' A Function to Create an Interactive Scatter Plot of the Principle Components Analysis of each of the Samples in an Expression Object
#'
#' This function allows you to plot PCA expression results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a Geo Expression object which can be obtained from the calculatePrincompPca() function
#' @keywords GEO
#' @export
#' @examples fig <- interactivePrcompPcaVariablesPlot(pcaData, "GSE18380")
#' @import plotly
#' @importFrom factoextra get_pca_var get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @seealso [calculatePrcompPca()] for Princomp PCA expression object, [extractPlatformGset()] for GEO object
interactivePrcompPcaVariablesPlot <- function(pcaData, geoAccessionCode) {
  library(plotly)
  library(factoextra)
  library(scales)
  variableStats <- get_pca_var(pcaData)
  eigenValue <- get_eigenvalue(pcaData)
  pcaData <- as.data.frame(unclass(pcaData$rotation))

  fig <- plot_ly(pcaData,x=~PC1,y=~PC2,text=rownames(pcaData), mode="markers", type = 'scatter'
                 ,marker=list(size=10, color = ~variableStats$contrib[,1]), name = rownames(pcaData))

  fig <- layout(fig,title= paste(geoAccessionCode, "PCA Variables Plot"),
                xaxis=list(title=paste("PC1", label_percent(accuracy=0.1)(eigenValue[1,2]/100))),
                yaxis=list(title=paste("PC2", label_percent(accuracy=0.1)(eigenValue[2,2]/100))))
  fig
}
