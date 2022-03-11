#' A Function to Create an Interactive Box and Whisker Plot
#' from an Expression Object
#'
#' This function allows you to plot expression data into an
#' interactive Box and Whisker Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Interactive Box-and-Whisker Plot
#' fig <- interactiveBoxAndWhiskerPlot(knnDataInput,
#' geoAccessionCode, platform)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
interactiveBoxAndWhiskerPlot <-
  function(ex) {
    ex <- as.data.frame(ex)

    fig1 <- plot_ly(type = "box", quartilemethod = "linear",
                   boxpoints = FALSE)
    i = 1
    for (col in names(ex)) {
      fig1 <-
        fig1 %>% add_trace(
          x = names(ex)[i],
          y = ex[, i],
          quartilemethod = "linear",
          name = names(ex)[i]
        )
      i <- i + 1
    }
    fig1 <-
      fig1 %>% layout(title = "Box And Whisker Plot")
    
    try(fig1 <- toWebGL(fig1))
    
    return(fig1)
  }

#' A Function to Create an Interactive Density Plot from
#' an Expression Object
#'
#' This function allows you to plot expression data into an
#' interactive Density Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Interactive Density Plot
#' fig <- interactiveDensityPlot(naOmitInput, geoAccessionCode,
#' platform)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
interactiveDensityPlot <- function(ex) {
  ex <- as.data.frame(ex)
  fig2 <-
    plot_ly(type = 'scatter',
            mode = 'lines')
  i <- 1
  for (col in names(ex)) {
    density <- density(ex[, i])
    fig2 <-
      fig2 %>% add_trace(x = density$x,
                        y = density$y,
                        name = col)
    i <- i + 1
  }

  fig2 <-
    fig2 %>% layout(
      title = "Density Plot",
      xaxis = list(title = 'Intensity'),
      yaxis = list(title = 'Density')
    )
  
  try(fig2 <- toWebGL(fig2))
  
  return(fig2)
}

#' A Function to Create an Interactive Three Dimensional
#' Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into an
#' interactive Three Dimensional Density Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # 3D Interactive Density Plot
#' fig <- interactiveThreeDDensityPlot(naOmitInput,
#' geoAccessionCode, platform)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
interactiveThreeDDensityPlot <-
  function(ex) {
    ex <- as.data.frame(ex)
    fig3 <-
      plot_ly(type = 'scatter3d',
              mode = 'lines')
    i <- 1
    for (col in names(ex)) {
      density <- density(ex[, i])
      fig3 <-
        fig3 %>% add_trace(
          x = density$x,
          y = i,
          z = density$y,
          name = col
        )
      i <- i + 1
    }

    fig3 <- fig3 %>% layout(
      title = "Density Plot",
      scene = list(
        xaxis = list(title = "Intensity"),
        yaxis = list(title = ""),
        zaxis = list(title = "Density")
      )
    )
    
    try(fig3 <- toWebGL(fig3))
    
    return(fig3)
  }

#' A Function to Create an Interactive UMAP Plot from
#' an Expression Object
#'
#' This function allows you to plot expression data into an
#' interactive UMAP Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param knn A integer representing the number of k nearest
#' neighbor's to use
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Interactive UMAP
#' knn <- 2
#' fig <- interactiveUmapPlot(naOmitInput, knn,
#' geoAccessionCode)
#' fig
#'
#' @import plotly umap
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
interactiveUmapPlot <- function(ex, knn) {
  ex <- ex[!duplicated(ex),]  # remove duplicates
  ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
  i <- 1
  fig5 <- plot_ly(type = 'scatter', mode = 'markers')
  for (row in row.names(ump$layout)) {
    fig5 <-
      fig5 %>% add_trace(x = ump$layout[i, ][1],
                        y = ump$layout[i, ][2],
                        name = row)
    i <- i + 1
  }
  fig5 <- fig5 %>% layout(title = (
    paste('UMAP plot, number of nearest neighbors used =', knn)
  ))
  
  try(fig5 <- toWebGL(fig5))
  
  return(fig5)
}

#' A Function to Create an Interactive Mean Variance Plot
#' from an Expression Object
#'
#' This function allows you to plot expression data into an
#' interactive Mean Variance Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Interactive Mean Variance Plot
#' fig <- interactiveMeanVariancePlot(naOmitInput,
#' geoAccessionCode, gsetData)
#' fig
#'
#' @import plotly limma
#' @importFrom stringr str_replace_all
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()]
#' for expression object, [extractPlatformGset()]
#' for GEO object
interactiveMeanVariancePlot <-
  function(ex, gset) {
    # Create linear model
    exData <- lmFit(ex)

    # Convert to dataframe
    exData <- as.data.frame(exData)
    exData["ID"] <- rownames(ex)

    if (is.null(gset)) {
      combineData <- exData
    } else if (length(gset@featureData@data)==0)
      {
      combineData <- exData
    } else {
      # Extract gene data
      geneData <- gset@featureData@data

      # Error handling to catch gset without featureData
      if (ncol(geneData) > 0) {
        geneData <- as.data.frame(geneData)
        combineData <- merge(exData, geneData, by = "ID")
        colnames(combineData) <-
          str_replace_all(colnames(combineData), " ", ".")
        combineData %>% filter("ID" %in% c(rownames(exData)))
      } else{
        combineData <- exData
      }
    }
    # Plot mean variance
    fig6 <-
      plot_ly(
        data = combineData,
        x = ~ Amean,
        y = ~ sigma,
        type = 'scatter',
        text =
          # Error handling to add gene data
          if ('ID' %in% colnames(combineData)) {
            if ('Gene.symbol' %in% colnames(combineData)) {
              if ('Gene.title' %in% colnames(combineData)) {
                if ('Gene.ID' %in% colnames(combineData)) {
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
                    'Amean: ',
                    Amean,
                    '<br>',
                    'Sigma: ',
                    sigma
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
                    'Amean: ',
                    Amean,
                    '<br>',
                    'Sigma: ',
                    sigma
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
                  'Amean: ',
                  Amean,
                  '<br>',
                  'Sigma: ',
                  sigma
                )
              }
            } else {
              ~ paste('ID: ',
                      ID,
                      '<br>',
                      'Amean: ',
                      Amean,
                      '<br>',
                      'Sigma: ',
                      sigma)
            }
          } else{
            ~ paste('Amean: ', Amean, '<br>', 'Sigma: ', sigma)
          }
        ,
        hoverinfo = text,
        mode = 'markers',
        marker = list(
          color = 'rgb(17, 157, 255)',
          size = 3,
          line = list(color = 'rgb(0, 0, 0)',
                      width = 1)
        )
      )
    fig6 <- fig6 %>% layout(title =
      'Mean Variance Plot')
    
    try(fig6 <- toWebGL(fig6))
      
    return(fig6)
  }

#' A Function to Create an Interactive Histogram of the
#' Prcomp Principle Components from the PCA outputs of an
#' Expression Object
#'
#' This function allows you to plot Prcomp PCA expression
#' results into an interactive Histogram of the Principle
#' Components
#' @param geoAccessionCode A character string representing
#' a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA
#' on a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Prcomp PCA analysis on KNN transformation
#' # expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' # Interactive Prcomp PCA Scree Plot
#' fig <- interactivePrcompPcaScreePlot(pcaPrcompDataInput,
#' geoAccessionCode)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()] for PCA expression object
interactivePrcompPcaScreePlot <-
  function(pcaData) {
    columnNames <- colnames(pcaData$x)
    proportionOfVariance <- pcaData$sdev ^ 2 / sum(pcaData$sdev ^ 2)
    pcaDataFrame <- data.frame(columnNames, proportionOfVariance)

    fig7 <-
      plot_ly(
        data = pcaDataFrame,
        x = ~ columnNames,
        y = ~ proportionOfVariance,
        type = "bar"
      ) %>%
      layout(
        title = "Scree Plot",
        xaxis = list(
          title = "Principal Components/Dimensions",
          categoryorder = "array",
          categoryarray = ~ columnNames
        ),
        yaxis = list(title = "Percentage of Explained Variances",
                     tickformat = ".0%")
      )

    try(fig7 <- toWebGL(fig7))
    
    return(fig7)
  }

#' A Function to Create an Interactive Histogram of the
#' princomp Principle Components from the PCA outputs of an
#' Expression Object
#'
#' This function allows you to plot princomp PCA expression
#' results into an interactive Histogram of the Principle
#' Components
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a
#' Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Princomp PCA analysis on KNN transformation expression data#
#' pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
#'
#' # Interactive Princomp PCA Scree Plot
#' fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput,
#' geoAccessionCode)
#' fig
#'
#' @import plotly
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()] for Princomp PCA expression object
interactivePrincompPcaScreePlot <-
  function(pcaData) {
    columnNames <-   colnames(pcaData$loadings)
    proportionOfVariance <- pcaData$sdev ^ 2 / sum(pcaData$sdev ^ 2)
    pcaDataFrame <- data.frame(columnNames, proportionOfVariance)
    fig8 <-
      plot_ly(
        data = pcaDataFrame,
        x = ~ columnNames,
        y = ~ proportionOfVariance,
        type = "bar"
      ) %>%
      layout(
        title = "Scree Plot",
        xaxis = list(
          title = "Principal Components/Dimensions",
          categoryorder = "array",
          categoryarray = ~ columnNames
        ),
        yaxis = list(title = "Percentage of Explained Variances",
                     tickformat = ".0%")
      )
    
    try(fig8 <- toWebGL(fig8))

    return(fig8)
  }

#' A Function to Create an Interactive Scatter Plot of the
#' princomp Principle Components Analysis of each of the
#' Genes in an Expression Object
#'
#' This function allows you to plot princomp PCA expression
#' results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on a
#' Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Princomp PCA analysis on KNN transformation expression data#
#' pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
#'
#' # Interactive Princomp PCA Individual Plot
#' fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput,
#' geoAccessionCode, gsetData)
#' fig
#'
#' @import plotly
#' @importFrom stringr str_replace_all
#' @importFrom factoextra get_pca_ind get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()]
#' for Princomp PCA expression object, [extractPlatformGset()]
#' for GEO object
interactivePrincompPcaIndividualsPlot <-
  function(pcaData, gset) {
    pcaDf <- data.frame(pcaData$scores)
    pcaDf <- transform(pcaDf)
    pcaDf["ID"] <- rownames(pcaDf)
    if (is.null(gset)) {
      combineData <- pcaDf
    } else if (length(gset@featureData@data)==0)
    {
      combineData <- pcaDf
    } else {
      geneData <- gset@featureData@data
      # Error handling for gset without featureData@data
      if (ncol(geneData) > 0) {
        geneData <- as.data.frame(geneData)
        combineData <- merge(pcaDf, geneData, by = "ID")
        combineData %>% filter("ID" %in% c(rownames(pcaDf)))
        colnames(combineData) <-
          str_replace_all(colnames(combineData), " ", ".")
      } else {
        combineData <- pcaDf
      }
    }


    individualsStats <- get_pca_ind(pcaData)
    eigenValue <- get_eigenvalue(pcaData)

    fig9 <-
      plot_ly(
        combineData,
        x =  ~ Comp.1,
        y =  ~ Comp.2,
        mode = "markers",
        type = 'scatter',
        text =
          if ('ID' %in% colnames(combineData)) {
            if ('Gene.symbol' %in% colnames(combineData)) {
              if ('Gene.title' %in% colnames(combineData)) {
                if ('Gene.ID' %in% colnames(combineData)) {
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
                    'Dimension 1: ',
                    Comp.1,
                    '<br>',
                    'Dimension 2: ',
                    Comp.2
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
                    'Dimension 1: ',
                    Comp.1,
                    '<br>',
                    'Dimension 2: ',
                    Comp.2
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
                  'Dimension 1: ',
                  Comp.1,
                  '<br>',
                  'Dimension 2: ',
                  Comp.2
                )
              }
            } else {
              ~ paste(
                'ID: ',
                ID,
                '<br>',
                'Dimension 1: ',
                Comp.1,
                '<br>',
                'Dimension 2: ',
                Comp.2
              )
            }
          } else {
            ~ paste('Dimension 1: ', Comp.1, '<br>',
                    'Dimension 2: ', Comp.2)
          }
        ,
        hoverinfo = text,
        marker = list(color = ~ individualsStats$cos2[, 1],
                      size = 3)
      )

    fig9 <-
      layout(
        fig9,
        title = "PCA Individuals Plot",
        xaxis = list(title = paste(
          "Comp.1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
        )),
        yaxis = list(title = paste(
          "Comp.2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
        ))
      )
    
    try(fig9 <- toWebGL(fig9))
    
    return(fig9)
  }

#' A Function to Create an Interactive Scatter Plot of the
#' princomp Principle Components Analysis of each of the
#' Samples in an Expression Object
#'
#' This function allows you to plot princomp PCA expression
#' results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param pcaData An object containing the results of PCA
#' on a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Princomp PCA analysis on KNN transformation expression data
#' pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
#'
#' # Interactive Princomp PCA Variables Plot
#' fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput,
#' geoAccessionCode)
#' fig
#'
#' @import plotly
#' @importFrom factoextra get_pca_var get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()]
#' for Princomp PCA expression object,
#' [extractPlatformGset()] for GEO object
interactivePrincompPcaVariablesPlot <-
  function(pcaData) {
    variableStats <- get_pca_var(pcaData)
    eigenValue <- get_eigenvalue(pcaData)
    pcaData <- as.data.frame(unclass(pcaData$loadings))

    fig10 <-
      plot_ly(
        pcaData,
        x =  ~ Comp.1,
        y =  ~ Comp.2,
        text = rownames(pcaData),
        mode = "markers",
        type = 'scatter'
        ,
        marker = list(size = 10, color = ~ variableStats$contrib[, 1]),
        name = rownames(pcaData)
      )

    fig10 <-
      layout(
        fig10,
        title ="PCA Variables Plot",
        xaxis = list(title = paste(
          "Comp.1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
        )),
        yaxis = list(title = paste(
          "Comp.2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
        ))
      )
    
    try(fig10 <- toWebGL(fig10))
    
    return(fig10)
  }

#' A Function to Create an Interactive Heat Map of the
#' Correlations between Samples
#'
#' This function allows you to plot an interactive heat map
#' of the correlations between samples in an expression object
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Correlation Matrix of samples
#' fig <- interactiveHeatMapPlot(naOmitInput)
#' fig
#'
#' @import plotly
#' @importFrom heatmaply heatmaply
#' @author Guy Hunt
#' @noRd
interactiveHeatMapPlot <- function(ex) {
  corMatrix <- cor(ex, use = "c")
  df <- data.frame(corMatrix[1, ])
  df <- df[-c(1)]
  i = 1
  while (i <= length(colnames(corMatrix))) {
    df[i] <- data.frame(corMatrix[, i])
    i <- i + 1
  }
  colnames(df) <- colnames(corMatrix)
  heatmapFig <- heatmaply(df)
  
  return(heatmapFig)
}

#' A Function to Create an Interactive Scatter Plot of the
#' prcomp Principle Components Analysis of each of the Genes
#' in an Expression Object
#'
#' This function allows you to plot prcomp PCA expression
#' results into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing
#' a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA
#' on a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Prcomp PCA analysis on KNN transformation expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' # Interactive Prcomp PCA Individual Plot
#' fig <- interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
#' geoAccessionCode, gsetData)
#' fig
#'
#' @import plotly
#' @importFrom stringr str_replace_all
#' @importFrom factoextra get_pca_ind get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrcompPca()]
#' for Princomp PCA expression object,
#' [extractPlatformGset()] for GEO object
interactivePrcompPcaIndividualsPlot <-
  function(pcaData, gset) {
    pcaDf <- data.frame(pcaData$x)
    pcaDf <- transform(pcaDf)
    pcaDf["ID"] <- rownames(pcaDf)
    if (is.null(gset)) {
      combineData <- pcaDf
    } else if (length(gset@featureData@data)==0)
    {
      combineData <- pcaDf
    } else {
      geneData <- gset@featureData@data
      # Error handling for gset without featureData@data
      if (ncol(geneData) > 0) {
        geneData <- as.data.frame(geneData)
        combineData <- merge(pcaDf, geneData, by = "ID")
        combineData %>% filter("ID" %in% c(rownames(pcaDf)))
        colnames(combineData) <-
          str_replace_all(colnames(combineData), " ", ".")
      } else {
        combineData <- pcaDf
      }
    }


    individualsStats <- get_pca_ind(pcaData)
    eigenValue <- get_eigenvalue(pcaData)


    fig11 <-
      plot_ly(
        combineData,
        x =  ~ PC1,
        y =  ~ PC2,
        mode = "markers",
        type = 'scatter',
        text =
          if ('ID' %in% colnames(combineData)) {
            if ('Gene.symbol' %in% colnames(combineData)) {
              if ('Gene.title' %in% colnames(combineData)) {
                if ('Gene.ID' %in% colnames(combineData)) {
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
                    'Dimension 1: ',
                    PC1,
                    '<br>',
                    'Dimension 2: ',
                    PC2
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
                    'Dimension 1: ',
                    PC1,
                    '<br>',
                    'Dimension 2: ',
                    PC2
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
                  'Dimension 1: ',
                  PC1,
                  '<br>',
                  'Dimension 2: ',
                  PC2
                )
              }
            } else {
              ~ paste('ID: ',
                      ID,
                      '<br>',
                      'Dimension 1: ',
                      PC1,
                      '<br>',
                      'Dimension 2: ',
                      PC2)
            }
          } else {
            ~ paste('Dimension 1: ', PC1, '<br>',
                    'Dimension 2: ', PC2)
          }
        ,
        hoverinfo = text,
        marker = list(color = ~ individualsStats$cos2[, 1],
                      size = 3)
      )

    fig11 <-
      layout(
        fig11,
        title = "PCA Individuals Plot",
        xaxis = list(title = paste(
          "PC1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
        )),
        yaxis = list(title = paste(
          "PC2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
        ))
      )
    
    try(fig11 <- toWebGL(fig11))
    
    return(fig11)
  }

#' A Function to Create an Interactive 3D Scatter Plot of the
#' prcomp Principle Components Analysis of each of the Genes
#' in an Expression Object
#'
#' This function allows you to plot prcomp PCA expression
#' results into an interactive 3D Scatter Plot
#' @param geoAccessionCode A character string representing
#' a GEO object for download and parsing
#' @param pcaData An object containing the results of PCA
#' on a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Prcomp PCA analysis on KNN transformation expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' # Interactive Prcomp PCA Individual Plot
#' fig <- interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
#' geoAccessionCode, gsetData)
#' fig
#'
#' @import plotly
#' @importFrom stringr str_replace_all
#' @importFrom factoextra get_pca_ind get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrcompPca()]
#' for Princomp PCA expression object,
#' [extractPlatformGset()] for GEO object
interactive3DPrcompPcaIndividualsPlot <-
  function(pcaData, gset) {
    pcaDf <- data.frame(pcaData$x)
    pcaDf <- transform(pcaDf)
    pcaDf["ID"] <- rownames(pcaDf)
    gset <- NULL
    if (is.null(gset)) {
      combineData <- pcaDf
    } else if (length(gset@featureData@data)==0)
    {
      combineData <- pcaDf
    } else {
      geneData <- gset@featureData@data
      # Error handling for gset without featureData@data
      if (ncol(geneData) > 0) {
        geneData <- as.data.frame(geneData)
        combineData <- merge(pcaDf, geneData, by = "ID")
        combineData %>% filter("ID" %in% c(rownames(pcaDf)))
        colnames(combineData) <-
          str_replace_all(colnames(combineData), " ", ".")
      } else {
        combineData <- pcaDf
      }
    }

    individualsStats <- get_pca_ind(pcaData)
    eigenValue <- get_eigenvalue(pcaData)

    fig12 <-
      plot_ly(
        combineData,
        x =  ~ PC1,
        y =  ~ PC2,
        z =  ~ PC3,
        mode = "markers",
        type = 'scatter3d',
        text =
          if ('ID' %in% colnames(combineData)) {
            if ('Gene.symbol' %in% colnames(combineData)) {
              if ('Gene.title' %in% colnames(combineData)) {
                if ('Gene.ID' %in% colnames(combineData)) {
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
                    'Dimension 1: ',
                    PC1,
                    '<br>',
                    'Dimension 2: ',
                    PC2
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
                    'Dimension 1: ',
                    PC1,
                    '<br>',
                    'Dimension 2: ',
                    PC2
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
                  'Dimension 1: ',
                  PC1,
                  '<br>',
                  'Dimension 2: ',
                  PC2
                )
              }
            } else {
              ~ paste('ID: ',
                      ID,
                      '<br>',
                      'Dimension 1: ',
                      PC1,
                      '<br>',
                      'Dimension 2: ',
                      PC2)
            }
          } else {
            ~ paste('Dimension 1: ', PC1, '<br>',
                    'Dimension 2: ', PC2)
          }
        ,
        hoverinfo = text,
        marker = list(color = ~ individualsStats$cos2[, 1],
                      size = 3)
      )

    fig12 <-
      layout(
        fig12,
        title = "PCA Individuals Plot",
        scene = list(
          xaxis = list(title = paste(
            "PC1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
          )),
          yaxis = list(title = paste(
            "PC2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
          )),
          zaxis = list(title = paste(
            "PC3", label_percent(accuracy = 0.1)(eigenValue[3, 2] / 100)
          ))
        )
      )
    
    try(fig12 <- toWebGL(fig12))
    
    return(fig12)
  }

#' A Function to Create an Interactive Scatter Plot
#' of the Principle Components Analysis of each of the
#' Samples in an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into an interactive Scatter Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on
#' a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#' # Perform Prcomp PCA analysis on KNN transformation expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' # Interactive Prcomp PCA Variables Plot
#' fig <- interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
#' geoAccessionCode)
#' fig
#'
#' @import plotly
#' @importFrom factoextra get_pca_var get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrcompPca()]
#' for Princomp PCA expression object,
#' [extractPlatformGset()] for GEO object
interactivePrcompPcaVariablesPlot <-
  function(pcaData) {
    variableStats <- get_pca_var(pcaData)
    eigenValue <- get_eigenvalue(pcaData)
    pcaData <- as.data.frame(unclass(pcaData$rotation))

    fig13 <-
      plot_ly(
        pcaData,
        x =  ~ PC1,
        y =  ~ PC2,
        text = rownames(pcaData),
        mode = "markers",
        type = 'scatter'
        ,
        marker = list(size = 10, color = ~ variableStats$contrib[, 1]),
        name = rownames(pcaData)
      )

    fig13 <-
      layout(
        fig13,
        title = "PCA Variables Plot",
        xaxis = list(title = paste(
          "PC1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
        )),
        yaxis = list(title = paste(
          "PC2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
        ))
      )
    
    try(fig13 <- toWebGL(fig13))
    
    return(fig13)
  }

#' A Function to Create an Interactive 3D Scatter Plot
#' of the Principle Components Analysis of each of the
#' Samples in an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into an interactive 3D Scatter Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param pcaData An object containing the results of PCA on
#' a Geo Expression object which can be obtained from the
#' calculatePrincompPca() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#' # Perform Prcomp PCA analysis on KNN transformation expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' # Interactive Prcomp PCA Variables Plot
#' fig <- interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
#' geoAccessionCode)
#' fig
#'
#' @import plotly
#' @importFrom factoextra get_pca_var get_eigenvalue
#' @importFrom scales label_percent
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrcompPca()]
#' for Princomp PCA expression object,
#' [extractPlatformGset()] for GEO object
interactive3DPrcompPcaVariablesPlot <-
  function(pcaData) {
    variableStats <- get_pca_var(pcaData)
    eigenValue <- get_eigenvalue(pcaData)
    pcaData <- as.data.frame(unclass(pcaData$rotation))

    fig14 <-
      plot_ly(
        pcaData,
        x =  ~ PC1,
        y =  ~ PC2,
        z =  ~ PC2,
        text = rownames(pcaData),
        mode = "markers",
        type = 'scatter3d'
        ,
        marker = list(size = 10, color = ~ variableStats$contrib[, 1]),
        name = rownames(pcaData)
      )

    fig14 <-
      layout(
        fig14,
        title = "PCA Variables Plot",
        scene = list(
        xaxis = list(title = paste(
          "PC1", label_percent(accuracy = 0.1)(eigenValue[1, 2] / 100)
        )),
        yaxis = list(title = paste(
          "PC2", label_percent(accuracy = 0.1)(eigenValue[2, 2] / 100)
        )),
        zaxis = list(title = paste(
          "PC3", label_percent(accuracy = 0.1)(eigenValue[3, 2] / 100)
        ))
      )
      )
    
    try(fig14 <- toWebGL(fig14))
    
    return(fig14)
  }
