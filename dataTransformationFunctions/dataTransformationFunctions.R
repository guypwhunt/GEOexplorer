#' A Function to Log Transform an Expression Object 
#'
#' This function allows you to log transform expression objects
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @param logTransformation Whether to auto-detect if log transformation is appropriate or to apply log transformation. Values can be "Auto-Detect" for auto detect, "Yes" to apply log transformation and "No" to not perform log transformation.
#' @keywords GEO
#' @export
#' @examples dataInput <- logTransformExpressionData(expressionData, "Auto-Detect")
logTransformExpressionData <- function(ex, logTransformation) {
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
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples autoLogInformation <- isLogTransformAutoApplied(expressionData)
isLogTransformAutoApplied <- function(ex) {
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
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @param knnTransformation Whether to apply KNN impute. This can be "Yes" or "No"
#' @keywords GEO
#' @export
#' @examples knnDataInput <- knnDataTransformation(dataInput, knnTransformation)
knnDataTransformation <- function(ex, knnTransformation) {
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
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples pcaDataInput <- pcaAnalysis(knnDataInput)
pcaAnalysis <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  pca <- prcomp(ex, scale = TRUE)
  return(pca)
}

#' A Function to Perform Principle Component Analysis on an Expression Object 
#'
#' A function to perform Princomp principle component analysis on an expression object. 
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples pcaPrincompDataInput <- pcaPrincompAnalysis(knnDataInput)
pcaPrincompAnalysis <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  pca <- princomp(ex, cor = TRUE)
  return(pca)
}

#' A Function to Removes Rows in an Expression Object that Contain Null Values 
#'
#' A function to perform remove rows that contain a null value from an expression object. 
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples naOmitInput <- naOmitTransformation(knnDataInput)
naOmitTransformation <- function(ex){
  library(impute)
  library(limma)
  library(factoextra)
  ex <- na.omit(ex)
  return(ex)}