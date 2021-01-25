# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE18388", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }



# Principal Component Analysis
run_pca <- function(X, population_names, popname1, popname2, isdebug) {
  Xpca         <- prcomp(t(X))
  raw_pca      <- get_pca_raw_data(Xpca)
  pc_plot_data <- get_pc_plotdata(Xpca, population_names, popname1, popname2,
                                  isdebug)
  pca_json     <- list(pc = raw_pca, pcdata = pc_plot_data)
  return(pca_json)
}

get_pca_raw_data <- function(Xpca) {
  s              <- summary(Xpca)
  exp_var        <- s$importance[2, ] * 100 # percentages Explained Variance
  cum_var        <- s$importance[3, ] * 100 # percentage Cumulative Variance
  pcnames        <- names(exp_var)          # PC names
  names(exp_var) <- NULL
  names(cum_var) <- NULL
  pca            <- list(pcnames = pcnames, expVar = exp_var, cumVar = cum_var)
  return (pca)
}

get_pc_plotdata <- function(Xpca, populations, popname1, popname2, isdebug) {
  Xscores           <- Xpca$x
  sample_names      <- rownames(Xscores)
  rownames(Xscores) <- NULL
  cols              <- colnames(Xscores)
  # Take each column of Xscore (temp variable y) & split based on the population
  Xscores           <- lapply(1:nrow(Xscores), function(y) {
    split(Xscores[, y], populations)
  })
  names(Xscores)    <- cols
  pc                <- unlist(Xscores, recursive = FALSE)
  # Split sample names by population and add to the final list
  complete.data  <- c(split(sample_names, populations), pc)
  # add the pop_names to the json list
  pop_names      <- list(group1 = popname1, group2 = popname2)
  complete.data  <- append(complete.data, pop_names)
  if (isdebug) cat("Overview: PCA has been calculated\n")
  return(complete.data)
}


print(run_pca(ex, 'a', 'b', 'c', FALSE))