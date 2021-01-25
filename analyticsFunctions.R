library(shiny)
library(GEOquery)
library(limma)
library(umap)
library(impute)

knnDataTransformation <- function(ex, knnTransformation) {
  if (knnTransformation == "Yes") {
    if (ncol(ex) == 2) {
    ex <- ex[complete.cases(ex), ] # KNN does not work when there are only 2 samples
    } else {
      ex <- ex[rowSums(is.na(ex)) != ncol(ex), ] # remove rows with missing data
      }
    # remove all zeros
    ex <- ex[rowSums(ex != 0) != 0,]
    
    # Replace missing value with calculated KNN value
    imputation <- impute.knn(ex)
    
    ex <- imputation$data
    return(ex)}
else if (knnTransformation == "No") {return(ex)}}