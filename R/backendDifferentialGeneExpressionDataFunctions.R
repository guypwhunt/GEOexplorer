#' A Function to Extract the Sample/Columns Names from an Expression Object
#'
#' This function extracts the sample/column names from an expression object
#' @param ex The GEO expression object which can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' @author Guy Hunt
#' @keywords internal
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
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#' group1 <- c()
#' group2 <- c()
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#'     group1 <- c(group1, name)
#'         i <- i +1
#'           } else {
#'           group2 <- c(group2, name)
#'           i <- i +1
#'           }
#'           }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' @author Guy Hunt
#' @keywords internal
#' @seealso [extractSampleNames()] for all the sample names
calculateEachGroupsSamples <- function(columnNames, group1, group2){
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
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#' group1 <- c()
#' group2 <- c()
#'
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#' group1 <- c(group1, name)
#'  i <- i +1
#'  } else {
#'  group2 <- c(group2, name)
#'  i <- i +1
#'  }
#'  }
#'
#'  # Select columns in group2
#'  column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#'  # Calculate gsms
#'  gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#'  # Convert P value adjustment
#'  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#'  adjustment <- convertAdjustment(pValueAdjustment)
#'  # Get fit 2
#'  limmaPrecisionWeights <- "Yes"
#'  forceNormalization <- "Yes"
#'  fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' @author Guy Hunt
#' @keywords internal
#' @seealso [extractExpressionData()] for expression object, [extractPlatformGset()] for GEO object, [calculateEachGroupsSamples()] for the string of integers indicating which group a sample belongs to
calculateDifferentialGeneExpression <- function(gsms, limmaPrecisionWeights, forceNormalization, gset, ex){
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
    v <- vooma(gset, design, plot=FALSE)
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
#' @examples
#'  # Convert P value adjustment
#'  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#'  adjustment <- convertAdjustment(pValueAdjustment)
#'
#' @author Guy Hunt
#' @keywords internal
convertAdjustment <- function(adjustment){
  # List of UI P value adjustments
  uiAdjustment <- c('Benjamini & Hochberg (False discovery rate)','Benjamini & Yekutieli','Bonferroni','Hochberg','Holm','None')

  # List of Backend P value adjustments
  backendAdjustment <- c("fdr", "BY", "bonferroni",'hochberg','holm','none')

  # Lookup table of UI and Backend P value adjustments
  adjustments <- backendAdjustment
  names(adjustments) <- uiAdjustment

  # Adjustment value
  adjustmentResult <- unname(adjustments[adjustment])

  return(adjustmentResult)
}

#' A Function to Create a Table of the Top Differentially Expressed Genes
#'
#' This function creates a table of the top differentially expressed genes
#' @param fit2 An object containing the differentially expressed genes analysis that can be obtained from the calculateDifferentialGeneExpression() function
#' @param adjustment A string character containing the adjustment to the P-value. The values can be: "fdr", "BY", "bonferroni", "hochberg", "holm", "hommel" or "none"
#' @keywords GEO
#' @import limma
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#' group1 <- c()
#' group2 <- c()
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#' group1 <- c(group1, name)
#' i <- i +1
#' } else {
#' group2 <- c(group2, name)
#' i <- i +1
#' }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Print Top deferentially expressed genes
#' tT <- calculateTopDifferentiallyExpressedGenes(fit2, adjustment)
#'
#' @author Guy Hunt
#' @keywords internal
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
calculateTopDifferentiallyExpressedGenes <- function(fit2, adjustment) {
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
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#' group1 <- c()
#' group2 <- c()
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#' group1 <- c(group1, name)
#' i <- i +1
#' } else {
#' group2 <- c(group2, name)
#' i <- i +1
#' }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#'
#' @author Guy Hunt
#' @keywords internal
calculateExclusiveColumns <- function(columns, inputColumns) {
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
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' # Extract platforms
#' platforms <- extractPlatforms(allGset)
#' platform <- platforms[1]
#'
#' # Extract the GEO2R data from the specified platform
#' gsetData <- extractPlatformGset(allGset, platform)
#'
#' # Extract expression data
#' expressionData <- extractExpressionData(gsetData)
#'
#' # Apply log transformation to expression data if necessary
#' logTransformation <- "Auto-Detect"
#' dataInput <- calculateLogTransformation(expressionData, logTransformation)
#'
#' # Perform KNN transformation on log expression data if necessary
#' knnDataInput <- calculateKnnImpute(dataInput, "Yes")
#'
#' # Extract experimental condition/sample names
#' columnNames <- extractSampleNames(expressionData)
#'
#' # Define Groups
#' numberOfColumns <- length(columnNames)
#' numberOfColumns <- numberOfColumns + 1
#' halfNumberOfColumns <- ceiling(numberOfColumns/2)
#' i <- 0
#' group1 <- c()
#' group2 <- c()
#' for (name in columnNames) {
#' if (i < halfNumberOfColumns) {
#' group1 <- c(group1, name)
#' i <- i +1
#' } else {
#' group2 <- c(group2, name)
#' i <- i +1
#' }
#' }
#'
#' # Select columns in group2
#' column2 <- calculateExclusiveColumns(columnNames, group1)
#' # Calculate gsms
#' gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
#'
#' # Convert P value adjustment
#' pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#' adjustment <- convertAdjustment(pValueAdjustment)
#'
#' # Get fit 2
#' limmaPrecisionWeights <- "Yes"
#' forceNormalization <- "Yes"
#' fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
#'
#' @author Guy Hunt
#' @keywords internal
#' @import limma
#' @seealso [calculateDifferentialGeneExpression()] for differential gene expression object
calculateDifferentialGeneExpressionSummary <- function(fit2, adjustment, significanceLevelCutOff) {
  dT <- decideTests(fit2, adjust.method=adjustment, p.value=significanceLevelCutOff)
  return(dT)
}

#' A Function to Calculate the Samples Selected in Each Group
#'
#' This function calculates the GSMS object for differential expression from a dataframe of the groups
#' @param groupDataFrame A dataframe containing the row number and group selection
#' @keywords GEO
#' @importFrom stringr str_remove_all
#' @author Guy Hunt
#' @keywords internal
calculateEachGroupsSamplesFromDataFrame <- function(groupDataFrame) {
  # Convert the input to a dataframe
  groupDataFrame <- as.data.frame(groupDataFrame)

  # For each row convert the UI codes to backend codes
  for (val in seq_len(nrow(groupDataFrame))) {
    if (groupDataFrame[val,1] == "N/A")
    {
      groupDataFrame[val,1] <- "X"
    } else if (groupDataFrame[val,1] == "Group 1"){
      groupDataFrame[val,1] <- 0
    } else if (groupDataFrame[val,1] == "Group 2"){
      groupDataFrame[val,1] <- 1
    }
  }

  # Convert the outputs to a string
  stringGroup <- toString(groupDataFrame[,1])

  # Remove all commas and spaces in the string
  stringGroup <- str_remove_all(str_remove_all(stringGroup, ",")," ")

  return(stringGroup)
}
