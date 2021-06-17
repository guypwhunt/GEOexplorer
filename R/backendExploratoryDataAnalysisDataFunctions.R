#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object from GEO when you know the specific Geo Accession code and platform
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform The study's platform
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples extractGeoData("GSE18380", "GPL4694")
#' @author Guy Hunt
extractGeoData <- function(geoAccessionCode, platform) {
  library(GEOquery)
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object from GEO when you only know the Geo Accession code
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param GSEMatrix A boolean telling GEOquery whether or not to use GSE Series Matrix files from GEO. The parsing of these files can be many orders-of-magnitude faster than parsing the GSE SOFT format files. Defaults to TRUE, meaning that the SOFT format parsing will not occur; set to FALSE if you for some reason need other columns from the GSE records.
#' @param getGPL A boolean defaulting to TRUE as to whether or not to download and include GPL information when getting a GSEMatrix file. You may want to set this to FALSE if you know that you are going to annotate your featureData using Bioconductor tools rather than relying on information provided through NCBI GEO. Download times can also be greatly reduced by specifying FALSE.
#' @param platformAnnotation A string defaulting to "NCBI generated" meaning true as to whether or not to use the Annotation GPL information. These files are nice to use because they contain up-to-date information remapped from Entrez Gene on a regular basis. However, they do not exist for all GPLs; in general, they are only available for GPLs referenced by a GDS. Input "Submitter supplied" for FALSE
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples allGset <- getGeoObject("GSE18380", GSEMatrix=TRUE, getGPL=TRUE, platformAnnotation = "NCBI generated")
#' @author Guy Hunt
getGeoObject <- function(geoAccessionCode, GSEMatrix=TRUE, getGPL=TRUE, platformAnnotation = "NCBI generated") {
  library(GEOquery)
  if (platformAnnotation == "Submitter supplied") {
    platformAnnotation <- FALSE
  } else if (platformAnnotation == "NCBI generated") {
    platformAnnotation <- TRUE
  } else {
    platformAnnotation <- TRUE
  }
  gset <- getGEO(geoAccessionCode, GSEMatrix=GSEMatrix, AnnotGPL=platformAnnotation) # getGPL=getGPL,
  return(gset)
}

#' A GEO Function to Obtain the Available Platforms of a GEO Object
#'
#' This function allows you to extract the platforms codes from a GEO object
#' @param gset The GEO object which can be obtained from the getGeoObject() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples platforms <- extractPlatforms(allGset)
#' @author Guy Hunt
#' @seealso [getGeoObject()] for GEO object
extractPlatforms <- function(gset) {
  library(GEOquery)
  platforms <- c()
  i <-1
  for(dataset in gset) {
    platforms[i] <- annotation(dataset)
    i <- i + 1
  }
  return(platforms)
}

#' A GEO Function to Obtain the GEO Object for a specific platform
#'
#' This function allows you to extract a platform's GEO object from a GEO object
#' @param gset The GEO object which can be obtained from the getGeoObject() function
#' @param platform The platform code
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples gsetData <- extractPlatformGset(gset, platforms[1])
#' @author Guy Hunt
#' @seealso [getGeoObject()] for GEO object
extractPlatformGset <- function(gset, platform) {
  library(GEOquery)
  if (length(gset) > 1) idx <- grep(platform[1], attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Function to Obtain the Experiment Information from a GEO Object
#'
#' This function allows you to extract experiment information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples experimentInformation <- extractExperimentInformation(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractExperimentInformation <- function(gset) {
  library(GEOquery)
  experimentalData <- experimentData(gset)
  return(experimentalData)
}

#' A GEO Function to Convert the Experiment Information Object into HTML
#'
#' This function allows you to convert experiment information into HTML
#' @param experimentData The experiment object obtained from the extractExperimentInformation() function
#' @keywords GEO
#' @export
#' @import GEOquery htmltools
#' @examples convertExperimentInformation <- convertExperimentInformation(experimentInformation)
#' @author Guy Hunt
#' @seealso [extractExperimentInformation()] for GEO object
convertExperimentInformation <- function(experimentData) {
  library(GEOquery)
  library(htmltools)
  name <- paste("<b>", "Author's Name:", "</b>", "<p>", experimentData@name, "</p>")
  lab <- paste("<b>", "Laboratory:", "</b>", "<p>", experimentData@lab, "</p>")
  contact <- paste("<b>", "Contact Information:", "</b>", "<p>", experimentData@contact, "</p>")
  title <- paste("<b>", "Paper Title:", "</b>", "<p>", experimentData@title, "</p>")
  abstract <- paste("<b>", "Abstract:", "</b>", "<p>", experimentData@abstract, "</p>")
  url <- paste("<b>", "Paper URL:", "</b>", "<p>", experimentData@url, "</p>")
  pubMedIds <- paste("<b>", "PubMed ID:", "</b>", "<p>", experimentData@pubMedIds, "</p>")
  html <- HTML(paste(title, name, lab, contact, url, pubMedIds, abstract, sep = ""))
}

#' A GEO Function to Extract Information on the Samples from a GEO object
#'
#' This function allows you to extract the studies sample information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples columnInfo <- extractSampleDetails(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractSampleDetails <- function(gset){
  library(GEOquery)
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data
  columnNames <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1")
  finalColumnNames <- c()
  i <- 1

  for (name in columnNames) {
    if (name %in% colnames(phenoData)) {
      finalColumnNames <- c(finalColumnNames,name)
    }
  }

  df <- data.frame(column=row.names(phenoData))

  for (name in finalColumnNames){
    df <- data.frame(df, phenoData[name])
  }
  return(df)
}

#' A GEO Function to Extract Expression Data from a GEO object
#'
#' This function allows you to extract the studies expression object from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples expressionData <- extractExpressionData(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractExpressionData <- function(gset) {
  library(GEOquery)
  ex <- exprs(gset)
  return(ex)}

#' A GEO Function to Extract Information on the Samples from a GEO object
#'
#' This function allows you to extract the studies sample information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @examples sampleInfo <- extractSampleInformation(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractSampleInformation <- function(gset) {
  library(GEOquery)
  sampleInfo <- pData(gset)
  return(sampleInfo)}

#' A Function to Log Transform an Expression Object
#'
#' This function allows you to log transform expression objects
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @param logTransformation Whether to auto-detect if log transformation is appropriate or to apply log transformation. Values can be "Auto-Detect" for auto detect, "Yes" to apply log transformation and "No" to not perform log transformation.
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples dataInput <- calculateLogTransformation(expressionData, "Auto-Detect")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculateLogTransformation <- function(ex, logTransformation = "Auto-Detect") {
  library(impute)
  library(limma)
  library(factoextra)
  # If log transformation is set to auto-detect
  if (logTransformation == "Auto-Detect"){
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)}
    return(ex)}

  # If log transformation is set to yes
  else if (logTransformation == "Yes") {
    #ex <- ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
    return(ex)}

  # If log transformation is set to no
  else if (logTransformation == "No") {
    return(ex)
  }
}

#' A Function to Determine if Log Transformation should Automatically be Applied
#'
#' This function allows you to determine if log transformation should be performed on an expression objects
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples autoLogInformation <- calculateAutoLogTransformApplication(expressionData)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculateAutoLogTransformApplication <- function(ex) {
  library(impute)
  library(limma)
  library(factoextra)
  # If log transformation is set to auto-detect
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    result <- "The auto-detect option applied log transformation."
  } else {
    result <- "The auto-detect option did not apply log transformation."
  }
  return(result)
}

#' A Function to Perform KNN Impute on an Expression Object
#'
#' A function to impute missing expression data, using nearest neighbor averaging.
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @param knnTransformation Whether to apply KNN impute. This can be "Yes" or "No"
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculateKnnImpute <- function(ex, knnTransformation) {
  library(impute)
  library(limma)
  library(factoextra)
  if (knnTransformation == "Yes") {
    if (ncol(ex) == 2) {
      ex <- ex[complete.cases(ex), ] # KNN does not work when there are only 2 samples
    } else {
      ex <- ex[rowSums(is.na(ex)) != ncol(ex), ] # remove rows with missing data
    }
    # remove all zeros (this was originally imported from GeoDrive but seems to be broken)
    # ex <- ex[rowSums(ex != 0) != 0,]

    # Replace missing value with calculated KNN value
    imputation <- impute.knn(ex)

    ex <- imputation$data
    return(ex)}
  else if (knnTransformation == "No") {return(ex)}}

#' A Function to Perform Principle Component Analysis on an Expression Object
#'
#' A function to perform prcomp principle component analysis on an expression object.
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples pcaDataInput <- calculatePca(knnDataInput)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculatePca <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  pca <- prcomp(ex, scale = TRUE)
  return(pca)
}

#' A Function to Perform Principle Component Analysis on an Expression Object
#'
#' A function to perform Princomp principle component analysis on an expression object.
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples pcaPrincompDataInput <- calculatePrincompPca(knnDataInput)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculatePrincompPca <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  pca <- princomp(ex, cor = TRUE)
  return(pca)
}

#' A Function to Removes Rows in an Expression Object that Contain Null Values
#'
#' A function to perform remove rows that contain a null value from an expression object.
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import impute limma factoextra
#' @examples naOmitInput <- calculateNaOmit(knnDataInput)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
calculateNaOmit <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  ex <- na.omit(ex)
  return(ex)}
