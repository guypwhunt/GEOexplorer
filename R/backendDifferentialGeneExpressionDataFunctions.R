#' A Function to Extract the Sample/Columns Names from an Expression Object
#'
#' This function extracts the sample/column names
#' from an expression object
#' @param ex The GEO expression object which
#' can be obtained from the extractExpressionData() function
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
#' @noRd
#' @seealso [extractExpressionData()]
#' for expression object
extractSampleNames <- function(ex) {
  columnNames <- colnames(ex)
  return(columnNames)
}

#' A Function to Calculate the Samples Selected in Each Group
#'
#' This function calculates the GSMS object
#' for differential expression from the sample
#' names and samples in each group
#' @param columnNames All the sample names in the
#' expression object which can be obtained from the extractSampleNames()
#' function
#' @param group1 The sample names in group 1.
#' This must not contain any sample names that are in group2
#' @param group2 The sample names in group 2.
#' This must not contain any sample names that are in group1
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
#' @noRd
#' @seealso [extractSampleNames()]
#' for all the sample names
calculateEachGroupsSamples <-
  function(columnNames, group1, group2) {
    lengthOfColumns <- sum(unlist(lapply(columnNames, length)))
    gsmsList <- vector(mode = "list", length = lengthOfColumns)
    i <- 1

    for (column in columnNames) {
      if (column %in% group1) {
        gsmsList[[i]] <- 0
        i <- i + 1
      } else if (column %in% group2) {
        gsmsList[[i]] <- 1
        i <- i + 1
      } else {
        gsmsList[[i]] <- "X"
        i <- i + 1
      }
    }
    gsms <- paste(gsmsList, collapse = '')
    return(gsms)
  }

#' A Function to Calculate the Differential Gene
#' EXpression between two groups
#'
#' This function calculates the differential
#' expression for two groups
#' @param gsms A string of integers indicating
#' which group a sample belongs to, which can be calculated form the
#' calculateEachGroupsSamples() function
#' @param limmaPrecisionWeights Whether to apply
#' limma precision weights (vooma)
#' @param forceNormalization Whether to force normalization
#' @param gset The GEO object which can be
#' obtained from the extractPlatformGset() function
#' @param ex The GEO expression object which
#' can be obtained from the extractExpressionData() function
#' @keywords GEO
#' @import limma
#' @importFrom edgeR DGEList as.matrix.DGEList calcNormFactors
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data
#' # if necessary
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
#'  gsms <- calculateEachGroupsSamples(columnNames,group1,
#'  group2)
#'
#'  # Convert P value adjustment
#'  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#'  adjustment <- convertAdjustment(pValueAdjustment)
#'  # Get fit 2
#'  limmaPrecisionWeights <- "Yes"
#'  forceNormalization <- "Yes"
#'  fit2 <- calculateDifferentialGeneExpression(gsms,
#'  limmaPrecisionWeights, forceNormalization, gsetData,
#'  knnDataInput)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()]
#' for expression object, [extractPlatformGset()]
#' for GEO object, [calculateEachGroupsSamples()]
#' for the string of integers indicating which group a
#' sample belongs to
calculateDifferentialGeneExpression <-
  function(gsms,
           input,
           all) {
    # Define results variable
    results <- NULL
    
    dataSource <- input$dataSource
    
    if (is.null(all$gsetData)){
      dataSource <- "Upload"
    }

    if (input$dataSetType == "Single") {
      if (dataSource == "GEO" & all$typeOfData == "Microarray") {
        # make proper column names to match toptable
        fvarLabels(all$gsetData) <- make.names(fvarLabels(all$gsetData))

        # Reduce the dimensionality of all$gsetData to that of ex
        all$gsetData <- all$gsetData[row.names(all$gsetData) %in% row.names(all$knnDataInput),]
        all$gsetData <- all$gsetData[, colnames(all$gsetData) %in% colnames(all$knnDataInput)]
      }

      # group membership for all samples
      sml <- strsplit(gsms, split = "")[[1]]
      sel <- which(sml != "X")
      sml <- sml[sel]
      all$knnDataInput <- all$knnDataInput[, sel]

      if (dataSource == "GEO" & all$typeOfData == "Microarray") {
        # Update gset data
        all$gsetData <- all$gsetData[, sel]
        exprs(all$gsetData) <- all$knnDataInput
      } 
      else if (all$typeOfData == "RNA Sequencing") {
        all$knnDataInput <- DGEList(all$knnDataInput, group = sml)
      }

      if (input$forceNormalization == "Yes") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          # normalize data
          exprs(all$gsetData) <- normalizeBetweenArrays(all$knnDataInput)
        } 
        else if (all$typeOfData == "RNA Sequencing") {
          all$knnDataInput = calcNormFactors(all$knnDataInput, method = "TMM")
        } else if (all$typeOfData == "Microarray") {
          all$knnDataInput <- normalizeBetweenArrays(all$knnDataInput)
        }
      }

      # assign samples to groups and set up design matrix
      gs <- factor(sml)
      groups <- make.names(c("Group1", "Group2"))
      levels(gs) <- groups

      if (dataSource == "GEO" & all$typeOfData == "Microarray") {
        # Update gset data
        all$gsetData$group <- gs

        # Create design
        design <- model.matrix(~ group + 0, all$gsetData)
      } 
      else if (all$typeOfData == "Microarray") {
        # Convert all$knnDataInput to expression dataset
        all$knnDataInput <- ExpressionSet(all$knnDataInput)
        all$knnDataInput$group <- gs
        # Create design
        design <- model.matrix(~ group + 0, all$knnDataInput)
      } else if (all$typeOfData == "RNA Sequencing") {
        all$knnDataInput$samples$group <- gs
        # Create design
        design <- model.matrix(~ group + 0, all$knnDataInput$samples)
      }

      colnames(design) <- levels(gs)

      if (input$limmaPrecisionWeights == "Yes") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          all$gsetData <- all$gsetData[complete.cases(exprs(all$gsetData)), ]

          # calculate precision weights and show plot of
          # mean-variance trend
          v <- vooma(all$gsetData, design, plot = FALSE)
          # attach gene annotations
          v$genes <- fData(all$gsetData)
        } 
        else if (all$typeOfData == "Microarray") {
          # Convert all$knnDataInput to matrix
          all$knnDataInput <- as.matrix(all$knnDataInput)
          all$knnDataInput <- all$knnDataInput[complete.cases(all$knnDataInput), ]
          # calculate precision weights
          v <- vooma(all$knnDataInput, design, plot = FALSE)
          # Add gene information
          v$genes <- as.data.frame(row.names(all$knnDataInput))
          colnames(v$genes) <- list("ID")
        }
        else if (all$typeOfData == "RNA Sequencing") {
          # calculate precision weights
          v <- voom(all$knnDataInput, design, plot = FALSE)
        }

        # fit linear model
        fit  <- lmFit(v)

        # Update results
        results$ex <- v

      } else if (input$limmaPrecisionWeights == "No") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          # fit linear model
          fit <- lmFit(all$gsetData, design)
          # Update results
          results$ex <- exprs(all$gsetData)
        } else if (all$typeOfData == "Microarray") {
          # fit linear model
          fit <- lmFit(all$knnDataInput, design)
          # attach gene annotations
          fit$genes <- as.data.frame(row.names(all$knnDataInput))
          # Update column name
          colnames(fit$genes) <- list("ID")

          # Update results as a matrix
          results$ex <- as.matrix(all$knnDataInput)
        } else if (all$typeOfData == "RNA Sequencing") {
          # fit linear model
          fit <- lmFit(as.matrix.DGEList(all$knnDataInput), design)
          # Update results
          results$ex <- all$knnDataInput$counts
        }
      }

      # set up contrasts of interest and recalculate
      # model coefficients
      cts <- paste(groups[1], groups[2], sep = "-")
      cont.matrix <- makeContrasts(contrasts = cts,
                                   levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)

      # compute statistics and table of top significant genes
      fit2 <- eBayes(fit2, 0.01)

      # Update results
      results$fit2 <- fit2

    } else if (input$dataSetType == "Combine")
    {
      if (dataSource == "GEO" & all$typeOfData == "Microarray") {
        # make proper column names to match toptable
        fvarLabels(all$gsetData) <- make.names(fvarLabels(all$gsetData))

        # Reduce the dimensionality of all$gsetData to that of all$knnDataInput
        all$gsetData <- all$gsetData[row.names(all$gsetData) %in% row.names(all$knnDataInput),]
        all$gsetData <- all$gsetData[, colnames(all$gsetData) %in% colnames(all$knnDataInput)]
      }

      # group membership for all samples
      sml <- strsplit(gsms, split = "")[[1]]
      sel <- which(sml != "X")
      sml <- sml[sel]
      all$knnDataInput <- all$knnDataInput[, sel]

      if (dataSource == "Upload") {
        if (all$typeOfData == "RNA Sequencing") {
          all$knnDataInput = DGEList(all$knnDataInput, group = sml)
        }
      }

      if (input$forceNormalization == "Yes") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          # normalize data
          all$knnDataInput <- normalizeBetweenArrays(all$knnDataInput)
        } else if (all$typeOfData == "RNA Sequencing") {
          all$knnDataInput = calcNormFactors(all$knnDataInput, method = "TMM")
        } else if (all$typeOfData == "Microarray") {
          all$knnDataInput <- normalizeBetweenArrays(all$knnDataInput)
        }
      }

      # assign samples to groups and set up design matrix
      gs <- factor(sml)
      groups <- make.names(c("Group1", "Group2"))
      levels(gs) <- groups

      if (dataSource == "GEO" & all$typeOfData == "Microarray") {
        # Update all$gsetData data
        all$knnDataInput <- ExpressionSet(all$knnDataInput)
        all$knnDataInput$group <- gs

        # Create design
        design <- model.matrix(~ group + 0, all$knnDataInput)
      } else if (all$typeOfData == "Microarray" ) {
        # Convert all$knnDataInput to expression dataset
        all$knnDataInput <- ExpressionSet(all$knnDataInput)
        all$knnDataInput$group <- gs
        # Create design
        design <- model.matrix(~ group + 0, all$knnDataInput)
      } else if (all$typeOfData == "RNA Sequencing") {
        all$knnDataInput$samples$group <- gs
        # Create design
        design <- model.matrix(~ group + 0, all$knnDataInput$samples)
      }

      colnames(design) <- levels(gs)

      if (input$limmaPrecisionWeights == "Yes") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          # Convert to matrix
          all$knnDataInput <- as.matrix(all$knnDataInput)
          all$knnDataInput <- all$knnDataInput[complete.cases(all$knnDataInput), ]

          # calculate precision weights and show plot of
          # mean-variance trend
          v <- vooma(all$knnDataInput, design, plot = FALSE)
          # attach gene annotations
          v$genes <- fData(all$gsetData)
        } else if (all$typeOfData == "Microarray") {
          # Convert all$knnDataInput to matrix
          all$knnDataInput <- as.matrix(all$knnDataInput)
          all$knnDataInput <- all$knnDataInput[complete.cases(all$knnDataInput), ]
          # calculate precision weights
          v <- vooma(all$knnDataInput, design, plot = FALSE)
          v$genes <- as.data.frame(row.names(all$knnDataInput))
          colnames(v$genes) <- list("ID")
        }
        else if (all$typeOfData == "RNA Sequencing") {
          # calculate precision weights
          v <- voom(all$knnDataInput, design, plot = FALSE)
          # fit linear model
          fit  <- lmFit(v)
          # Udate results
          results$ex <- v
        }

        # fit linear model
        fit  <- lmFit(v)

        # Update results
        results$ex <- v

      } else if (input$limmaPrecisionWeights == "No") {
        if (dataSource == "GEO" & all$typeOfData == "Microarray") {
          # fit linear model
          fit <- lmFit(all$knnDataInput, design)
          # Update gene information
          fit$genes <- fData(all$gsetData)
          # Update results
          results$ex <- as.matrix(all$knnDataInput)
        } else if (all$typeOfData == "Microarray") {
          # fit linear model
          fit <- lmFit(all$knnDataInput, design)
          # attach gene annotations
          fit$genes <- as.data.frame(row.names(all$knnDataInput))
          # Update column name
          colnames(fit$genes) <- list("ID")

          # Update results as a matrix
          results$ex <- as.matrix(all$knnDataInput)
        } else if (input$dataSetType == "RNA Sequencing") {
          # fit linear model
          fit <- lmFit(as.matrix.DGEList(all$knnDataInput), design)
          # Update results
          results$ex <- all$knnDataInput$counts
        }
      }

      # set up contrasts of interest and recalculate
      # model coefficients
      cts <- paste(groups[1], groups[2], sep = "-")
      cont.matrix <- makeContrasts(contrasts = cts,
                                   levels = design)
      fit2 <- contrasts.fit(fit, cont.matrix)

      # compute statistics and table of top significant genes
      fit2 <- eBayes(fit2, 0.01)

      # Update results
      results$fit2 <- fit2
    }

    return(results)
  }

#' A Function to Convert the UI P-Value Adjustment
#' into the Backend P-Value Adjustment
#'
#' This function converts the P-value adjustment
#' value from the UI into the value required by the backend
#' @param adjustment A string character containing the
#' adjustment to the P-value.
#' The values can be: "Benjamini & Hochberg (False discovery rate)",
#' "Benjamini & Yekutieli", "Bonferroni", "Hochberg",
#' "Holm", "Hommel"or "None"
#' @keywords GEO
#' @examples
#'  # Convert P value adjustment
#'  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
#'  adjustment <- convertAdjustment(pValueAdjustment)
#'
#' @author Guy Hunt
#' @noRd
convertAdjustment <- function(adjustment) {
  # List of UI P value adjustments
  uiAdjustment <-
    c(
      'Benjamini & Hochberg (False discovery rate)',
      'Benjamini & Yekutieli',
      'Bonferroni',
      'Hochberg',
      'Holm',
      'None'
    )

  # List of Backend P value adjustments
  backendAdjustment <-
    c("fdr", "BY", "bonferroni", 'hochberg', 'holm', 'none')

  # Lookup table of UI and Backend P value adjustments
  adjustments <- backendAdjustment
  names(adjustments) <- uiAdjustment

  # Adjustment value
  adjustmentResult <- unname(adjustments[adjustment])

  return(adjustmentResult)
}

#' A Function to Create a Table of the Top
#' Differentially Expressed Genes
#'
#' This function creates a table of the top
#' differentially expressed genes
#' @param fit2 An object containing the differentially
#' expressed genes analysis that can be obtained from
#' the calculateDifferentialGeneExpression() function
#' @param adjustment A string character containing
#' the adjustment to the P-value. The values can be:
#' "fdr", "BY", "bonferroni", "hochberg", "holm", "hommel"
#' or "none"
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if
#' # necessary
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Print Top deferentially expressed genes
#' tT <- calculateTopDifferentiallyExpressedGenes(fit2, adjustment)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object
calculateTopDifferentiallyExpressedGenes <-
  function(fit2, adjustment, numberOfGenes = 250) {
    tT <- topTable(fit2,
                   adjust.method = adjustment,
                   sort.by = "B",
                   number = numberOfGenes)
    columnNamesList <- c()
    optionalColumnNamesList <-
      c(
        "ID",
        "adj.P.Val",
        "P.Value",
        "t",
        "B",
        "logFC",
        "Gene.symbol",
        "Gene.title",
        "Gene.ID",
        "GB_LIST",
        "SPOT_ID",
        "RANGE_GB",
        "RANGE_STRAND",
        "RANGE_START",
        "GB_ACC",
        "GB_RANGE",
        "SEQUENCE"
      )
    tTColumnNames <- colnames(tT)
    for (columnName in optionalColumnNamesList)
    {
      if (columnName %in% tTColumnNames)
      {
        columnNamesList <- c(columnNamesList, columnName)
      }
    }
    tT["ID"] <- rownames(tT)
    tT <- subset(tT, select = columnNamesList)
    return(tT)
  }

#' A Function to Create a List of Columns that
#' Have Not Been Selected
#'
#' This function creates a list of columns that have not
#' been selected in the UI
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if
#' # necessary
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
#' @noRd
calculateExclusiveColumns <- function(columns, inputColumns) {
  columns1Input <- c()
  for (value in columns) {
    if (value %in% inputColumns) {

    } else {
      columns1Input = c(columns1Input, value)
    }
  }
  return(columns1Input)
}

#' A Function to Create an Object Containing if Each Gene
#' is Unregulated, Down Regulated or has a Similar Level of
#' Expression between the Groups
#'
#' This function creates an object containing if each gene
#' is unreguate, downregulated or not
#' @param fit2 An object containing the differentially
#' expressed genes analysis that can be obtained from the
#' calculateDifferentialGeneExpression() function
#' @param adjustment A string character containing the
#' adjustment to the P-value. The values can be:
#' "fdr", "BY", "bonferroni", "hochberg", "holm", "hommel"
#' or "none"
#' @param significanceLevelCutOff A float indicating the
#' P-value cutoff. The values can be between 0 and 1
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
#' dataInput <- calculateLogTransformation(expressionData,
#' logTransformation)
#'
#' # Perform KNN transformation on log expression data if
#' # necessary
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
#' fit2 <- calculateDifferentialGeneExpression(gsms,
#' limmaPrecisionWeights, forceNormalization, gsetData,
#' knnDataInput)
#'
#' # Summarize test results as "up", "down" or "not expressed"
#' significanceLevelCutOff <- 0.05
#' dT <- calculateDifferentialGeneExpressionSummary(fit2,
#' adjustment, significanceLevelCutOff)
#'
#' @author Guy Hunt
#' @noRd
#' @import limma
#' @seealso [calculateDifferentialGeneExpression()]
#' for differential gene expression object
calculateDifferentialGeneExpressionSummary <-
  function(fit2,
           adjustment,
           significanceLevelCutOff) {
    dT <-
      decideTests(fit2, adjust.method = adjustment,
                  p.value = significanceLevelCutOff)
    return(dT)
  }

#' A Function to Calculate the Samples Selected in Each Group
#'
#' This function calculates the GSMS object for
#' differential expression from a dataframe of the groups
#' @param groupDataFrame A dataframe containing
#' the row number and group selection
#' @keywords GEO
#' @importFrom stringr str_remove_all
#' @author Guy Hunt
#' @noRd
calculateEachGroupsSamplesFromDataFrame <-
  function(groupDataFrame) {
    # Convert the input to a dataframe
    groupDataFrame <- as.data.frame(groupDataFrame)

    # For each row convert the UI codes to backend codes
    for (val in seq_len(nrow(groupDataFrame))) {
      if (groupDataFrame[val, 1] == "N/A")
      {
        groupDataFrame[val, 1] <- "X"
      } else if (groupDataFrame[val, 1] == "Group 1") {
        groupDataFrame[val, 1] <- 0
      } else if (groupDataFrame[val, 1] == "Group 2") {
        groupDataFrame[val, 1] <- 1
      }
    }

    # Convert the outputs to a string
    stringGroup <- toString(groupDataFrame[, 1])

    # Remove all commas and spaces in the string
    stringGroup <-
      str_remove_all(str_remove_all(stringGroup, ","), " ")

    return(stringGroup)
  }

#' A Function to Calculate the Samples Selected in Each Group
#'
#' This function calculates the GSMS object
#' @author Guy Hunt
#' @noRd
calculateEachGroupsSamplesGsms <-
  function(columnInfo, groupOne, groupTwo) {
    # Update all rows to X
    columnInfo$group <- "X"

    # Add group 1 columns
    columnInfo[groupOne,]$group <- "0"

    # Add group 2 columns
    columnInfo[groupTwo,]$group <- "1"

    # colapse columns to character
    gsms <- paste(columnInfo$group, collapse = '')

    return(gsms)
  }
