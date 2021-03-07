library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE18388", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "00001111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Group1","Group2"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)

adjustmentCalculation <- function(adjustment)
{
  if (adjustment == "Benjamini & Hochberg (False discovery rate)"){
    adjustment <- "fdr"
  } else if (adjustment == "Benjamini & Yekutieli"){
    adjustment <- "BY"
  } else if (adjustment == "Bonferroni"){
    adjustment <- "bonferroni"
  } else if (adjustment == "Hochberg"){
    adjustment <- "hochberg"
  } else if (adjustment == "Holm"){
    adjustment <- "holm"
  } else if (adjustment == "Hommel"){
    adjustment <- "hommel"
  } else if (adjustment == "None"){
    adjustment <- "none"
  }
  return(adjustment)
}

topDifferentiallyExpressedGenesTable <- function(fit2, adjustment) {
  tT <- topTable(fit2, adjust=adjustment, sort.by="B", number=250)
  tT["ID"] <- rownames(tT)
  tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC"))
  return(tT)
}


histogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
}

# summarize test results as "up", "down" or "not expressed"
dT <- function(fit2, adjustment) {
  dT<- decideTests(fit2, adjust.method=adjustment, p.value=0.05)
  return(dT)
}

options <- list("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Hochberg", "Holm", "Hommel", "None")

for (x in options){
  a <- adjustmentCalculation(x)
  print(a)
  b <- topDifferentiallyExpressedGenesTable(fit2, a)
  print(b)
  c <- histogramPlot(fit2, a)
  print(c)
  d <- dT(fit2, a)
  print(d)
}

