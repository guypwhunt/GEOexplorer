#' A Function to Extract the Sample/Columns Names from an Expression Object
#'
#' This function extracts the sample/column names from an expression object
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @examples columnNames <- extractSampleNames(expressionData)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object
extractSampleNames <- function(ex) {
  columnNames <- colnames(ex)
  return(columnNames)
}

#' A Function to Calculate the Samples Selected in Each Group
#'
#' This function calculates the GSMS object for differential expression from the sample names and samples in each group
#' @param columnNames All the sample names in the expression object which can be obtained from the extractSampleNames() function
#' @param group1 The sample names in group 1. This must not contain any sample names that are in group2
#' @param group2 The sample names in group 2. This must not contain any sample names that are in group1
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples gsms <- calculateEachGroupsSamples(columnNames,c("GSM455528", "GSM455541", "GSM455542", "GSM455543", "GSM455578", "GSM455610", "GSM455782"), c("GSM455783", "GSM455784", "GSM455785", "GSM455786", "GSM455787"))
#' @author Guy Hunt
#' @seealso [extractSampleNames()] for all the sample names
calculateEachGroupsSamples <- function(columnNames, group1, group2){
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
#' @param gsms A string of integers indicating which group a sample belongs to, which can be calculated form the calculateEachGroupsSamples() function
#' @param limmaPrecisionWeights Whether to apply limma precision weights (vooma)
#' @param forceNormalization Whether to force normalization
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gset, ex)
#' @author Guy Hunt
#' @seealso [extractExpressionData()] for expression object, [extractPlatformGset()] for GEO object, [calculateEachGroupsSamples()] for the string of integers indicating which group a sample belongs to
calculateDifferentialGeneExpression <- function(gsms, limmaPrecisionWeights, forceNormalization, gset, ex){
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  # make proper column names to match toptable
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  # group membership for all samples
  sml <- strsplit(gsms, split="")[[1]]

  # Reduce the dimensionality of gset to that of ex
  gset <- gset[row.names(gset) %in% row.names(ex), ]
  gset <- gset[,colnames(gset) %in% colnames(ex)]

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

#' A Function to Convert the UI P-Value Adjustment into the Backend P-Value Adjustment
#'
#' This function converts the P-value adjustment value from the UI into the value required by the backend
#' @param adjustment A string character containing the adjustment to the P-value. The values can be: "Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Hochberg", "Holm", "Hommel"or "None"
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples adjustment <- convertAdjustment("Benjamini & Hochberg (False discovery rate)")
#' @author Guy Hunt
convertAdjustment <- function(adjustment){
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
  } else {
    adjustment <- "none"
  }
  return(adjustment)
}

#' A Function to Create a Table of the Top Differentially Expressed Genes
#'
#' This function creates a table of the top differentially expressed genes
#' @param fit2 An object containing the differentially expressed genes analysis that can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A string character containing the adjustment to the P-value. The values can be: "fdr", "BY", "bonferroni", "hochberg", "holm", "hommel" or "none"
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples tT <- calculateTopDifferentiallyExpressedGenes(fit2, "fdr")
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
calculateTopDifferentiallyExpressedGenes <- function(fit2, adjustment) {
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

#' A Function to Create a List of Columns that Have Not Been Selected
#'
#' This function creates a list of columns that have not been selected in the UI
#' @param columns A list of all columns within the expression object
#' @param inputColumns A list of the columns that have been selected
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples column2 <- calculateExclusiveColumns(c("GSM455528", "GSM455541", "GSM455542", "GSM455543"), c("GSM455528", "GSM455541"))
#' @author Guy Hunt
calculateExclusiveColumns <- function(columns, inputColumns) {
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

#' A Function to Create an Object Containing if Each Gene is Unregulated, Down Regulated or has a Similar Level of Expression between the Groups
#'
#' This function creates an object containing if each gene is unreguate, downregulated or not
#' @param fit2 An object containing the differentially expressed genes analysis that can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A string character containing the adjustment to the P-value. The values can be: "fdr", "BY", "bonferroni", "hochberg", "holm", "hommel" or "none"
#' @param significanceLevelCutOff A float indicating the P-value cutoff. The values can be between 0 and 1
#' @keywords GEO
#' @export
#' @import GEOquery limma umap data.table
#' @examples dT <- calculateDifferentialGeneExpressionSummary(fit2, "fdr", 0.05)
#' @author Guy Hunt
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
calculateDifferentialGeneExpressionSummary <- function(fit2, adjustment, significanceLevelCutOff) {
  library(GEOquery)
  library(limma)
  library(umap)
  library(data.table)
  dT <- decideTests(fit2, adjust.method=adjustment, p.value=significanceLevelCutOff)
  return(dT)
}
