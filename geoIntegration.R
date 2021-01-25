library(GEOquery)
library(limma)

getGeoData <- function(geoAccessionCode, platform) {
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

extractGeoData <- function(gset, logTransformation) {
  # If log transformation is set to auto-detect
  if (logTransformation == "Auto-Detect"){
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)}
    return(ex)} 
  
  # If log transformation is set to yes
  else if (logTransformation == "Yes") {
    ex <- exprs(gset)
    #ex <- ex[which(ex <= 0)] <- NaN
    ex <- log2(ex)
    return(ex)}
  
  # If log transformation is set to no
  else if (logTransformation == "No") {
    ex <- exprs(gset)
    return(ex)
  }
  }
  