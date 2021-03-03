#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

extractColumns <- function(ex) {
  columnNames <- colnames(ex) 
  return(columnNames)
}

calculateGsms <- function(columnNames,group1, group2, group3){
  lengthOfColumns <- sum(unlist(lapply(columnNames, length)))
  gsmsList <- vector(mode = "list", length = lengthOfColumns)
  i <- 1
  
  for (column in columnNames){
    if (column %in% group1) {
      gsmsList[[i]] <- 0
      i <- i+ 1
    } else if (column %in% group2) {
      gsmsList[[i]] <- 1
      i <- i+ 1
    } else if (column %in% group3) {
      gsmsList[[i]] <- 2
      i <- i+ 1
    }
  }
  gsms <- paste(gsmsList, collapse = '')
  return(gsms)
}

differentialGeneExpression <- function(gset, gsms) {
  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))
  
  # group membership for all samples
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
  return(fit2)
}

topDifferentiallyExpressedGenesTable2 <- function(geoAccessionCode) {
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
  return(tT)
}

topDifferentiallyExpressedGenesTable <- function(geoAccessionCode) {
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, AnnotGPL=TRUE)
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
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
  
  #tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
  return(tT)
}



histogramPlot <- function(fit2) {
  tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
}

# summarize test results as "up", "down" or "not expressed"
dT <- function(fit2) {
  dT<- decideTests(fit2, adjust.method="fdr", p.value=0.05)
  return(dT)
}

# Venn diagram of results
vennDiagramPlot <- function(dT) {
  vennDiagram(dT, circle.col=palette())
}

qqPlot <- function(fit2) {
# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
}

volcanoPlot <- function(fit2, dT, ct) {
# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
}

mdPlot <- function(fit2, dT, ct) {
  # MD plot (log fold change vs mean log expression)
  # highlight statistically significant (p-adj < 0.05) probes
  plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  abline(h=0)
}