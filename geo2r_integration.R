# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

# This contains the GEO accession used to access different gene expression data
gset <- getGEO("GSE18388", GSEMatrix =TRUE, getGPL=FALSE)
print('GSET')
print(gset)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
print('GSET2')
print(gset)

ex <- exprs(gset)
print('EX')
print(ex)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }
print('EX2')
print(ex)

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), " Box and Whisker Plot", sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE18388")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)