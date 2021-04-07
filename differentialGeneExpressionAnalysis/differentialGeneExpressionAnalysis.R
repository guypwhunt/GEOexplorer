#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(data.table)

extractColumns <- function(ex) {
  columnNames <- colnames(ex)
  return(columnNames)
}

calculateGsms <- function(columnNames,group1, group2){
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
    } else {
      gsmsList[[i]] <- "X"
      i <- i+ 1
    }
  }
  gsms <- paste(gsmsList, collapse = '')
  return(gsms)
}

calculateFit2 <- function(geoAccessionCode, platform, gsms, logTransformation, limmaPrecisionWeights, forceNormalization, knnTransformation, gset, platformAnnotation = "NCBI generated"){
  if (platformAnnotation == "Submitter supplied") {
    platformAnnotation <- FALSE
  } else if (platformAnnotation == "NCBI generated") {
    platformAnnotation <- TRUE
  } else {
    platformAnnotation <- TRUE
  }
  #gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, AnnotGPL=platformAnnotation)
  #gset <- getPlatformGset(gset, platform)

  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  # group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]

  ex <- extractExpressionData(gset)

  ex <- logTransformExpressionData(ex, logTransformation)

  ex <- knnDataTransformation(ex, knnTransformation)

  sel <- which(sml != "X")
  sml <- sml[sel]
  gset <- gset[ ,sel]
  ex <- ex[ ,sel]
  exprs(gset) <- ex

  if(forceNormalization == "Yes"){
    exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
  }



  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("Group1","Group2"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)

  if (limmaPrecisionWeights == "Yes"){
    nall <- nrow(gset)
    gset <- gset[complete.cases(exprs(gset)), ]

    # calculate precision weights and show plot of mean-variance trend
    v <- vooma(gset, design, plot=T)
    # OR weights by group
    # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
    v$genes <- fData(gset) # attach gene annotations

    # fit linear model
    fit  <- lmFit(v)
  } else if (limmaPrecisionWeights == "No"){
    fit <- lmFit(gset, design)  # fit linear model
  }

  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)

  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  return(fit2)
}

differentialGeneExpression <- function(gset, ex, gsms, limmaPrecisionWeights, forceNormalization) {
  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  # group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]

  # filter out excluded samples (marked as "X")
  sel <- which(sml != "X")
  sml <- sml[sel]
  gset <- gset[ ,sel]
  ex <- ex[ ,sel]
  exprs(gset) <- ex

  if(forceNormalization == "Yes"){
    exprs(gset) <- normalizeBetweenArrays(exprs(gset)) # normalize data
  }

  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("Group1","Group2"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)

  if (limmaPrecisionWeights == "Yes"){
    nall <- nrow(gset)
    gset <- gset[complete.cases(exprs(gset)), ]

    # calculate precision weights and show plot of mean-variance trend
    v <- vooma(gset, design, plot=T)
    # OR weights by group
    # v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))
    v$genes <- fData(gset) # attach gene annotations

    # fit linear model
    fit  <- lmFit(v)
  } else if (limmaPrecisionWeights == "No"){
    fit <- lmFit(gset, design)  # fit linear model
  }

  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)

  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  return(fit2)
}

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
  columnNamesList <- c()
  optionalColumnNamesList <- c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title", "Gene.ID", "GB_LIST", "SPOT_ID", "RANGE_GB", "RANGE_STRAND", "RANGE_START", "GB_ACC", "GB_RANGE", "SEQUENCE")
  tTColumnNames <- colnames(tT)
  for(columnName in optionalColumnNamesList)
  {
    if (columnName %in% tTColumnNames)
    {
      columnNamesList <- c(columnNamesList, columnName)
    }
  }
  tT["ID"] <- rownames(tT)
  tT <- subset(tT, select=columnNamesList)
  return(tT)
}


histogramPlot <- function(fit2, adjustment) {
  tT2 <- topTable(fit2, adjust=adjustment, sort.by="B", number=Inf)
  hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
       ylab = "Number of genes", main = "P-adj value distribution")
}

# summarize test results as "up", "down" or "not expressed"
calculateDT <- function(fit2, adjustment, significanceLevelCutOff) {
  dT <- decideTests(fit2, adjust.method=adjustment, p.value=significanceLevelCutOff)
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

exclusiveColumns <- function(columns, inputColumns) {
  columns1Input <- c()
  for (value in columns) {
    if(value %in% inputColumns) {

    } else {
      columns1Input = c(columns1Input, value)
    }
    }
    return(columns1Input)
}

heatmapPlot <- function(fit2, ex) {
  full_results <- topTable(fit2, number=Inf)
  topN <- 20
  ##
  ids_of_interest <- mutate(full_results, Rank = 1:n()) %>%
    filter(Rank < topN) %>%
    pull(ID)

  gene_names <- mutate(full_results, Rank = 1:n()) %>%
    filter(Rank < topN) %>%
    pull(ID)

  ## Get the rows corresponding to ids_of_interest and all columns
  gene_matrix <- ex[ids_of_interest,]

  pheatmap(gene_matrix)

  pheatmap(gene_matrix,
           scale="row")
}

