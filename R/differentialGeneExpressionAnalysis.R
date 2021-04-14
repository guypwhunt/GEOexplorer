#' A Function to Extract Columns Names from an Expression Object
#'
#' This function extracts the column names from an expression object
#' @param ex A GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples columnNames <- extractColumns(expressionData)
#' @author Guy Hunt
extractColumns <- function(ex) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  columnNames <- colnames(ex)
  return(columnNames)
}

#' A Function to Calculate the GSMS Object
#'
#' This function calculates the GSMS object for differential expression from the sample names and samples in each group
#' @param columnNames All the sample names in the expression object which can be obtained from the extractColumns() function
#' @param group1 The sample names in group 1. This must not contain any sample names that are in group2
#' @param group2 The sample names in group 2. This must not contain any sample names that are in group1
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples gsms <- calculateGsms(columnNames,c("GSM455528", "GSM455541", "GSM455542", "GSM455543", "GSM455578", "GSM455610", "GSM455782"), c("GSM455783", "GSM455784", "GSM455785", "GSM455786", "GSM455787"))
#' @author Guy Hunt
calculateGsms <- function(columnNames, group1, group2){
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
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

#' A Function to Calculate the Differential Gene EXpression between two groups
#'
#' This function calculates the differential expression for two groups
#' @param gsms TBD
#' @param logTransformation Whether to auto-detect if log transformation is appropriate or to apply log transformation. Values can be "Auto-Detect" for auto detect, "Yes" to apply log transformation and "No" to not perform log transformation.
#' @param limmaPrecisionWeights Whether to apply limma precision weights (vooma)
#' @param forceNormalization Whether to force normalization
#' @param knnTransformation Whether to fill in missing values using Knn
#' @param gset TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples fit2 <- calculateFit2(gsms, logTransformation, limmaPrecisionWeights, forceNormalization, knnTransformation)
#' @author Guy Hunt
calculateFit2 <- function(gsms, logTransformation, limmaPrecisionWeights, forceNormalization, knnTransformation, gset){
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  # group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]

  # This might need to be looked into
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

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
differentialGeneExpression <- function(gset, ex, gsms, limmaPrecisionWeights, forceNormalization) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
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

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
adjustmentCalculation <- function(adjustment){
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
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

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
topDifferentiallyExpressedGenesTable <- function(fit2, adjustment) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
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

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
exclusiveColumns <- function(columns, inputColumns) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  columns1Input <- c()
  for (value in columns) {
    if(value %in% inputColumns) {

    } else {
      columns1Input = c(columns1Input, value)
    }
    }
    return(columns1Input)
}

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
calculateDT <- function(fit2, adjustment, significanceLevelCutOff) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  dT <- decideTests(fit2, adjust.method=adjustment, p.value=significanceLevelCutOff)
  return(dT)
}

#' TBD
#'
#' TBD
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples TBD
#' @author Guy Hunt
vennDiagramPlot <- function(dT) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  vennDiagram(dT, circle.col=palette())
}
