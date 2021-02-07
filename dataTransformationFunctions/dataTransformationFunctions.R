library(impute)
library(limma)

extractExpressionData <- function(gset) {
  ex <- exprs(gset)
  return(ex)}

logTransformExpressionData <- function(ex, logTransformation) {
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

knnDataTransformation <- function(ex, knnTransformation) {
  if (knnTransformation == "Yes") {
    if (ncol(ex) == 2) {
    ex <- ex[complete.cases(ex), ] # KNN does not work when there are only 2 samples
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

pcaAnalysis <- function(ex){
  pca <- prcomp(ex, scale = TRUE)
  return(pca)
}

naOmitTransformation <- function(ex){
  ex <- na.omit(ex)
  return(ex)}