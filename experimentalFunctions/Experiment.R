library(GEOquery)
library(gplots)


gse = getGEO('GSE976')[[1]]
heatmapPlot <- function(ex) {
sdN = 3
sds = apply(log2(ex+0.0001),1,sd)
heatmap.2(log2(ex+0.0001)[sds>sdN,],trace='none',scale='row')
}

