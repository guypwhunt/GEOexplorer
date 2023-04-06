#' A Function to Create a Box and Whisker Plot from an
#' Expression Object
#'
#' This function allows you to plot expression data into a
#' Box and Whisker Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Non-Interactive Box-and-Whisker Plot
#' fig <- nonInteractiveBoxAndWhiskerPlot(ex = knnDataInput,
#' geoAccessionCode = geoAccessionCode, platform = platform)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
nonInteractiveBoxAndWhiskerPlot <-
  function(ex) {
    par(mar = c(7, 4, 2, 1))
    fig <-
      boxplot(
        ex,
        boxwex = 0.7,
        #xlab = "Experimental Conditions",
        ylab = "Gene Expression",
        notch = TRUE,
        main = "",
        outline = FALSE,
        las = 2
      )
    return(fig)
  }

#' A Function to Create a Density Plot from an Expression Object
#'
#' This function allows you to plot expression data into a
#' Density Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Non-Interactive Density Plot
#' fig <- nonInteractiveDensityPlot(ex = naOmitInput,
#' geoAccessionCode = geoAccessionCode, platform = platform)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
nonInteractiveDensityPlot <-
  function(ex,
           geoAccessionCode = "",
           platform = "") {
    par(mar = c(4, 4, 2, 1))
    title <-
      paste(geoAccessionCode, "/", platform, " value distribution", sep = "")
    fig <- plotDensities(ex, main = title, legend = FALSE)
    return(fig)
  }

#' A Function to Create a Mean Variance Plot from an
#' Expression Object
#'
#' This function allows you to plot expression data into a
#' Mean Variance Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Non-Interactive Mean Variance Plot
#' fig <- nonInteractiveMeanVariancePlot(naOmitInput,
#' geoAccessionCode)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
nonInteractiveMeanVariancePlot <-
  function(ex,
           geoAccessionCode = "",
           platform = "") {
    fig <-
      plotSA(lmFit(ex), main = paste("Mean variance trend,", geoAccessionCode))
    return(fig)
  }

#' A Function to Create a UMAP Plot from an Expression Object
#'
#' This function allows you to plot expression data into a
#' UMAP Plot
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param platform A character string representing the study's
#' platform
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
#' @keywords GEO
#' @import umap limma
#' @importFrom maptools pointLabel
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
#' # Non-Interactive UMAP
#' knn <- 2
#' fig <- nonInteractiveUmapPlot(naOmitInput, knn,
#' geoAccessionCode)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
nonInteractiveUmapPlot <-
  function(ex,
           knn,
           geoAccessionCode = "",
           platform = "") {
    ex <- ex[!duplicated(ex),]  # remove duplicates
    ump <- umap(t(ex), n_neighbors = knn, random_state = 123)
    plot(
      ump$layout,
      main = paste("UMAP plot, number of nearest neighbors used =",
                   knn),
      xlab = "",
      ylab = "",
      pch = 20,
      cex = 1.5
    )
    fig <- pointLabel(
      ump$layout,
      labels = rownames(ump$layout),
      method = "SANN",
      cex = 0.6
    )
    return(fig)
  }

#' A Function to Create a Histogram of the Principle Components
#' from the PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into a Histogram of the Principle Components
#' @param pcaEx A PCA object which can be obtained from the
#' calculatePrincompPca() function
#' @keywords GEO
#' @importFrom factoextra fviz_eig
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
#' # Non-Interactive Princomp PCA Scree Plot
#' fig <- nonInteractivePcaScreePlot(pcaPrincompDataInput)
#' fig
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()]
#' for Princomp Pca expression object
nonInteractivePcaScreePlot <- function(pcaEx) {
  fig <- fviz_eig(pcaEx)
  return(fig)
}

#' A Function to Create an Individuals Scatter Plot from the
#' PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into an Individuals Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the
#' calculatePrincompPca() function
#' @keywords GEO
#' @importFrom factoextra fviz_pca_ind
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
#' # Non-Interactive Princomp PCA Individual Plot
#' fig <- nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
#' fig
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaIndividualsPlot <- function(pcaEx) {
  fig <- fviz_pca_ind(
    pcaEx,
    # Color by the quality of representation
    col.ind = "cos2",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    geom = "point",
    repel = TRUE     # Avoid text overlapping
  )
  return(fig)
}

#' A Function to Create an Variables Scatter Plot from the
#' PCA outputs of an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into an Variables Scatter Plot
#' @param pcaEx A PCA object which can be obtained from the
#' calculatePrincompPca() function
#' @keywords GEO
#' @importFrom factoextra fviz_pca_var
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
#' # Non-Interactive Princomp PCA Variables Plot
#' fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
#' fig
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaVariablesPlot <- function(pcaEx) {
  fig <- fviz_pca_var(
    pcaEx,
    # Color by contributions to the PC
    col.var = "contrib",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE     # Avoid text overlapping
  )
  return(fig)
}

#' A Function to Create a Scatter Plot that contains
#' both the Individuals and Variables from the PCA outputs
#' of an Expression Object
#'
#' This function allows you to plot PCA expression results
#' into an Scatter Plot contains both the Individuals and
#' Variables
#' @param pcaEx A PCA object which can be obtained from the
#' calculatePrincompPca() function
#' @keywords GEO
#' @importFrom factoextra fviz_pca_biplot
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
#' # Non-Interactive Princomp PCA Variables Plot
#' fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
#' fig
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [calculatePrincompPca()] for Princomp Pca expression object
nonInteractivePcaBiplotPlot <- function(pcaEx) {
  fig <- fviz_pca_biplot(
    pcaEx,
    repel = TRUE,
    # Variables color
    col.var = "#2E9FDF",
    geom = "point",
    # Individuals color
    col.ind = "#696969",
  )
  return(fig)
}

#' A Function to Create a Correlation Matrix that contains
#' both the Correlations between Samples
#'
#' This function allows you to plot a heatmap of the
#' correlations between experimental conditions
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
#' @keywords GEO
#' @import pheatmap
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Non-Interactive Heatmap
#' fig <- nonInteractiveCorrelationMatrixPlot(naOmitInput)
#'
#' @author Guy Hunt
#' @noRd
nonInteractiveCorrelationMatrixPlot <- function(ex) {
  corMatrix <- cor(ex, use = "c")
  fig <- pheatmap(corMatrix)
  return(fig)
}
