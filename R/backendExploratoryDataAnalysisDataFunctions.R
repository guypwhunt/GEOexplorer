#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object
#' from GEO when you know the specific Geo Accession
#' code and platform
#' @param geoAccessionCode A character string representing
#' a GEO object for download and parsing
#' @param platform The study's platform
#' @keywords GEO
#' @importFrom GEOquery getGEO
#' @examples extractGeoData("GSE18380", "GPL4694")
#' @author Guy Hunt
#' @noRd
extractGeoData <- function(geoAccessionCode, platform) {
  gset <- getGEO(geoAccessionCode, GSEMatrix = TRUE,
                 getGPL = FALSE)
  if (length(gset) > 1)
    idx <- grep(platform, attr(gset, "names"))
  else
    idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object
#' from GEO when you only know the Geo Accession code
#' @param geoAccessionCode A character string representing a
#' GEO object for download and parsing
#' @param GSEMatrix A boolean telling GEOquery whether or not
#' to use GSE Series Matrix files from GEO. The parsing of
#' these files can be many orders-of-magnitude faster than
#' parsing the GSE SOFT format files. Defaults to TRUE,
#' meaning that the SOFT format parsing will not occur;
#' set to FALSE if you for some reason need other columns from
#' the GSE records.
#' @param getGPL A boolean defaulting to TRUE as to whether
#' or not to download and include GPL information when
#' getting a GSEMatrix file. You may want to set this to
#' FALSE if you know that you are going to annotate your
#' featureData using Bioconductor tools rather than relying
#' on information provided through NCBI GEO. Download times
#' can also be greatly reduced by specifying FALSE.
#' @param platformAnnotation A string defaulting to
#' "NCBI generated" meaning true as to whether or not to
#' use the Annotation GPL information. These files are nice
#' to use because they contain up-to-date information remapped
#' from Entrez Gene on a regular basis. However, they do not
#' exist for all GPLs; in general, they are only available for
#' GPLs referenced by a GDS. Input "Submitter supplied" for
#' FALSE
#' @keywords GEO
#' @importFrom GEOquery getGEO
#' @examples
#' # Get the GEO data for all platforms
#' geoAccessionCode <- "GSE18388"
#' allGset <- getGeoObject(geoAccessionCode)
#'
#' @author Guy Hunt
#' @noRd
getGeoObject <-
  function(geoAccessionCode,
           GSEMatrix = TRUE,
           getGPL = TRUE,
           platformAnnotation = "NCBI generated") {
    if (platformAnnotation == "Submitter supplied") {
      platformAnnotation <- FALSE
    } else if (platformAnnotation == "NCBI generated") {
      platformAnnotation <- TRUE
    } else {
      platformAnnotation <- TRUE
    }
    # Remove white space from geoAccessionCode
    gset <- tryCatch({
      getGEO(geoAccessionCode,
             GSEMatrix = GSEMatrix,
             AnnotGPL = platformAnnotation)
    }, error = function(x) {
      backupGset <-
        getGEO(geoAccessionCode,
               GSEMatrix = GSEMatrix,
               getGPL = FALSE)
      return(backupGset)
    })
    return(gset)
  }

#' A GEO Function to Obtain the Available Platforms of a
#' GEO Object
#'
#' This function allows you to extract the platforms codes
#' from a GEO object
#' @param gset The GEO object which can be obtained from
#' the getGeoObject() function
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
#' @author Guy Hunt
#' @noRd
#' @seealso [getGeoObject()] for GEO object
extractPlatforms <- function(gset) {
  platforms <- c()
  i <- 1
  for (dataset in gset) {
    platforms[i] <- annotation(dataset)
    i <- i + 1
  }
  return(platforms)
}

#' A GEO Function to Obtain the GEO Object for a specific
#' platform
#'
#' This function allows you to extract a platform's GEO
#' object from a GEO object
#' @param gset The GEO object which can be obtained from
#' the getGeoObject() function
#' @param platform The platform code
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
#' @author Guy Hunt
#' @noRd
#' @seealso [getGeoObject()] for GEO object
extractPlatformGset <- function(gset, platform) {
  if (length(gset) > 1)
    idx <- grep(platform[1], attr(gset, "names"))
  else
    idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Function to Obtain the Experiment Information from
#' a GEO Object
#'
#' This function allows you to extract experiment information
#' from a GEO object
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' # Extract the experiment information
#' experimentInformation <- extractExperimentInformation(gsetData)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractPlatformGset()] for GEO object
extractExperimentInformation <- function(gset) {
  experimentalData <- experimentData(gset)
  return(experimentalData)
}

#' A GEO Function to Convert the Experiment Information
#' Object into HTML
#'
#' This function allows you to convert experiment information
#' into HTML
#' @param experimentData The experiment object obtained from
#' the extractExperimentInformation() function
#' @keywords GEO
#' @importFrom htmltools HTML
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
#' # Extract the experiment information
#' experimentInformation <- extractExperimentInformation(gsetData)
#'
#' # Extract Experiment Information
#' extractedExperimentInformation <- convertExperimentInformation(
#' experimentInformation)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExperimentInformation()] for GEO object
convertExperimentInformation <- function(experimentData) {
  name <-
    paste("<b>",
          "Author's Name:",
          "</b>",
          "<p>",
          experimentData@name,
          "</p>")
  lab <-
    paste("<b>",
          "Laboratory:",
          "</b>",
          "<p>",
          experimentData@lab,
          "</p>")
  contact <-
    paste("<b>",
          "Contact Information:",
          "</b>",
          "<p>",
          experimentData@contact,
          "</p>")
  title <-
    paste("<b>",
          "Paper Title:",
          "</b>",
          "<p>",
          experimentData@title,
          "</p>")
  abstract <-
    paste("<b>",
          "Abstract:",
          "</b>",
          "<p>",
          experimentData@abstract,
          "</p>")
  url <-
    paste("<b>", "Paper URL:", "</b>", "<p>",
          experimentData@url, "</p>")
  pubMedIds <-
    paste("<b>",
          "PubMed ID:",
          "</b>",
          "<p>",
          experimentData@pubMedIds,
          "</p>")
  html <-
    HTML(paste(title, name, lab, contact, url,
               pubMedIds, abstract, sep = ""))
}

#' A GEO Function to Extract Information on the Samples
#'  from a GEO object
#'
#' This function allows you to extract the studies sample
#' information from a GEO object
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' # Get column Details
#' columnInfo <- extractSampleDetails(gsetData)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractPlatformGset()] for GEO object
extractSampleDetails <- function(gset) {
  # Extract the experimental conditions
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data

  # Define the desired column names and variable
  columnNames <-
    c("title",
      "source_name_ch1",
      "characteristics_ch1",
      "characteristics_ch1.1")
  finalColumnNames <- c()
  rejectedColumnNames <- c()
  i <- 1

  # Identify the columns that are available and not available
  for (name in columnNames) {
    if (name %in% colnames(phenoData)) {
      finalColumnNames <- c(finalColumnNames, name)
    } else {
      rejectedColumnNames <- c(rejectedColumnNames, name)
    }
  }

  # Create the dataframe with the available columns
  df <- data.frame(column = row.names(phenoData))

  for (name in finalColumnNames) {
    df <- data.frame(df, phenoData[name])
  }

  # Create the dataframe with the non-available columns
  if (!is.null(rejectedColumnNames)) {
    rejectedDf <- data.frame(column = row.names(phenoData))
    row.names(rejectedDf) <- row.names(phenoData)
    x <- 2
    for (name in rejectedColumnNames) {
      rejectedDf[x] <- NA
      x <- x +1
    }

    # Remove the column column
    rejectedDf <- subset(rejectedDf, select = -column)

    # Update the column names
    colnames(rejectedDf) <- rejectedColumnNames

    # Bind the two dataframes
    df <- cbind(df, rejectedDf)
  }
  return(df)
}

#' A GEO Function to Extract Expression Data from a GEO object
#'
#' This function allows you to extract the studies expression
#' object from a GEO object
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' @author Guy Hunt
#' @noRd
#' @seealso [extractPlatformGset()] for GEO object
extractExpressionData <- function(gset) {
  ex <- exprs(gset)

  # Delete columns that are all na unless there is only one column
  if (ncol(ex) > 1) {
    # Deletes columns for which all values are na
    try(ex <- ex[, colSums(is.na(ex)) < nrow(ex)])

    # Deletes rows for which all values are na
    try(ex <- ex[rowSums(is.na(ex)) < nrow(ex), ])
  }

  # Ensure the expression data contains numerical data
  if (typeof(ex) != "double") {
    try(ex2 <- apply(ex, 2, as.double))
    try(row.names(ex2) <- row.names(ex))
    try(ex <- ex2)
  }

  return(ex)
}

#' A GEO Function to Extract Information on the Samples from
#' a GEO object
#'
#' This function allows you to extract the studies sample
#' information from a GEO object
#' @param gset The GEO object which can be obtained from the
#' extractPlatformGset() function
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
#' #' # Extract Sample Information
#' sampleInfo <- extractSampleInformation(gsetData)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractPlatformGset()] for GEO object
extractSampleInformation <- function(gset) {
  sampleInfo <- pData(gset)
  return(sampleInfo)
}

#' A Function to Log Transform an Expression Object
#'
#' This function allows you to log transform expression objects
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
#' @param logTransformation Whether to auto-detect if log
#' transformation is appropriate or to apply log transformation. Values can be
#' "Auto-Detect" for auto detect, "Yes" to apply log transformation and "No" to
#' not perform log transformation.
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
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculateLogTransformation <-
  function(ex, logTransformation = "Auto-Detect") {
    # If log transformation is set to auto-detect
    if (logTransformation == "Auto-Detect") {
      # log2 transform
      qx <- as.numeric(# Define the quantile values
        quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
                 na.rm = TRUE))
      LogC <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0)
      if (LogC) {
        ex[which(ex <= 0)] <- NaN
        ex <- log2(ex)
      }
      return(ex)
    }

    # If log transformation is set to yes
    else if (logTransformation == "Yes") {
      #ex <- ex[which(ex <= 0)] <- NaN
      ex <- log2(ex)
      return(ex)
    }

    # If log transformation is set to no
    else if (logTransformation == "No") {
      return(ex)
    }
  }

#' A Function to Determine if Log Transformation should
#' Automatically be Applied
#'
#' This function allows you to determine if log transformation
#' should be performed on an expression objects
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Is log transformation auto applied
#' autoLogInformation <- calculateAutoLogTransformApplication(
#' expressionData
#' )
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculateAutoLogTransformApplication <- function(ex) {
  # If log transformation is set to auto-detect
  qx <- as.numeric(# Define the quantile values
    quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    result <- "The auto-detect option applied log transformation."
  } else {
    result <- "The auto-detect option did not apply log transformation."
  }
  return(result)
}

#' A Function to Perform KNN Impute on an Expression Object
#'
#' A function to impute missing expression data, using
#' nearest neighbor averaging.
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
#' @param knnTransformation Whether to apply KNN impute.
#' This can be "Yes" or "No"
#' @keywords GEO
#' @import impute
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
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculateKnnImpute <- function(ex, knnTransformation) {
  if (knnTransformation == "Yes") {
    # Check if there are less than 3 samples
    if (ncol(ex) > 3) {
      # If there are less than three samples remove rows with
      # blank values
      try(ex <- ex[rowSums(is.na(ex)) != ncol(ex),])

      # Replace missing value with calculated KNN value
      imputation <- impute.knn(ex)

      # Extract data from imputation object
      ex <- imputation$data
    } else {
      stop("At least 3 columns are required for KNN imputation")
    }
  }
  return(ex)
}

#' A Function to Perform Principle Component Analysis on an
#' Expression Object
#'
#' A function to perform prcomp principle component analysis
#' on an expression object.
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Prcomp PCA analysis on KNN transformation
#' # expression data
#' pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculatePrcompPca <- function(ex) {
  pca <- prcomp(ex, scale = TRUE)
  return(pca)
}

#' A Function to Perform Principle Component Analysis on an
#' Expression Object
#'
#' A function to perform Princomp principle component analysis
#' on an expression object.
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' # Perform Princomp PCA analysis on KNN transformation
#' # expression data
#' pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculatePrincompPca <- function(ex) {
  pca <- princomp(ex, cor = TRUE)
  return(pca)
}

#' A Function to Removes Rows in an Expression Object that
#' Contain Null Values
#'
#' A function to perform remove rows that contain a null
#' value from an expression object.
#' @param ex The GEO expression object which can be obtained
#' from the extractExpressionData() function
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
#' # Remove all incomplete rows
#' naOmitInput <- calculateNaOmit(knnDataInput)
#'
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExpressionData()] for expression object
calculateNaOmit <- function(ex) {
  ex <- na.omit(ex)
  return(ex)
}

#' A function to read the input CSVs
#'
#' This function allows you to read the input CSVs
#' @param file A object CSV object
#' @keywords csv
#' @importFrom utils read.csv
#' @examples # Extract CSV
#' rnaExpressionData <- readCsvFile("17cbedyM0aXEJD47wKoL6HhMhsg0w2Vop")
#' @author Guy Hunt
#' @noRd
readCsvFileFromGoogleDocs <- function(id) {
  # read the file
  dF <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download",
                         id))

  return(dF)
}

#' A function to read the input CSVs
#'
#' This function allows you to read the input CSVs
#' @param file A object CSV object
#' @keywords csv
#' @importFrom utils read.csv
#' @examples # Extract CSV
#' rnaExpressionData <- readCsvFile("C:/Users/guypw/OneDrive/Documents/
#' GEOexplorer/R/testScripts/exampleGeneExpressionCsv.csv")
#' @author Guy Hunt
#' @noRd
readCsvFile <- function(file) {
  # read the file
  dF <- read.csv(file = file, check.names = FALSE)

  return(dF)
}

#' A function to preprocess user uploaded gene expression data
#'
#' This function allows you to preprocess user uploaded gene expression data
#' @param file A object CSV object
#' @keywords csv
#' @import Biobase
#' @examples # Extract CSV
#' rnaExpressionData <- readCsvFile("C:/Users/guypw/OneDrive/Documents/
#' GEOexplorer/R/testScripts/exampleGeneExpressionCsv.csv")
#'
#' # Get a list of all the columns
#' columns <- extractSampleNames(rnaExpressionData)
#' @author Guy Hunt
#' @importFrom stringr str_remove_all
#' @noRd
preProcessGeneExpressionData <- function(expressionData) {
  # Remove quote marks from column names
  colnames(expressionData) <- str_remove_all(colnames(expressionData), "'")
  colnames(expressionData) <- str_remove_all(colnames(expressionData), '"')

  # Replace row names with gene names
  row.names(expressionData) <- expressionData$Genes

  # Remove quote marks from row names
  row.names(expressionData) <- str_remove_all(row.names(expressionData), "'")
  row.names(expressionData) <- str_remove_all(row.names(expressionData), '"')


  # Delete gene column
  expressionData <- expressionData[ , !(colnames(expressionData) == "Genes")]

  # Convert to matrix
  expressionData <- data.matrix(expressionData)

  # Replace NULL values with NA
  expressionData[expressionData == "NULL"] <- NA

  return(expressionData)
}

#' A function to preprocess user uploaded gene expression data
#'
#' This function allows you to convert expression data to counts per million
#' @param expressionData A object containing the gene expression data
#' @keywords rnaSeq
#' @importFrom edgeR cpm
#' @examples # Define Variables
#' geoAccessionCode <- "GSE63310"
#'
#' # Source and extract expression data
#' rnaExpressionData <- extractGeoSupFiles(geoAccessionCode)
#'
#' # Update Sample Names
#' columns <- calculateSampleNames(rnaExpressionData)
#'
#' # Raw counts are converted to counts-per-million (CPM)
#' cpm <- cpm(rnaExpressionData, "Yes")
#' @author Guy Hunt
#' @noRd
calculateCountsPerMillion <- function(rnaExpressionData, applyCpm) {
  if (applyCpm == "Yes"){
  # Raw counts are converted to counts-per-million (CPM)
  cpm <- cpm(rnaExpressionData)} else {
    cpm <- rnaExpressionData
  }
  return(cpm)
}


#' A function to search for datasets in GEO
#'
#' This function allows you to search for datasets in GEO
#' @param expressionData A object containing the gene expression data
#' @keywords rnaSeq
#' @importFrom htmltools HTML
#' @importFrom stringr str_replace
#' @importFrom httr GET
#' @importFrom XML xmlParse xmlToList xmlTreeParse
#' @author Guy Hunt
#' @noRd
searchGeo <- function(searchTerm, resultsLimit) {
  # Convert resultsLimit to a string
  resultsLimit <- as.character(resultsLimit)

  # Define the two URLs
  firstURL <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?
  db=gds&term=searchTerm&retmax=resultsLimit&usehistory=y'
  secondURL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?
  db=gds&version=2.0&query_key=X&WebEnv=ENTER_WEBENV_PARAMETER_HERE"

  # Update the first url with the search parameters
  firstURL <- str_replace(firstURL, "searchTerm", searchTerm)
  firstURL <- str_replace(firstURL, "resultsLimit", resultsLimit)

  # Perform the API call
  firstResults <- GET(firstURL)

  # Parse the response
  firstResultsData <- xmlParse(firstResults)
  firstResultsData <- xmlToList(firstResultsData)

  # Update the second url with the search parameters
  secondURL <- str_replace(secondURL, "X", firstResultsData$QueryKey)
  secondURL <- str_replace(secondURL, "ENTER_WEBENV_PARAMETER_HERE",
                           firstResultsData$WebEnv)

  # Perform the API call
  secondResults <- GET(secondURL)

  # Parse the response
  #secondResultsData <- xmlTreeParse(secondResults)

  # Convert from list to character
  #secondResultsData <- unlist(secondResultsData)

  # Convert to HTML
  #secondResultsData <- HTML(secondResultsData)

  return(secondResults)
  #return(secondResultsData$doc$children)
}

#' A function to apply batch correction
#'
#' This function allows you to apply batch correction
#' @keywords geneExpression
#' @import limma
#' @importFrom sva ComBat
#' @author Guy Hunt
#' @noRd
calculateBatchCorrection <- function(expressionData, expressionData2,
                            combinedExpressionData, batchCorrection) {
  if (batchCorrection != "None"){
    # Define the batches
    batchOne <- replicate(ncol(expressionData), "Batch1")
    batchTwo <- replicate(ncol(expressionData2), "Batch2")
    finalBatch <- c(batchOne, batchTwo)

    if (batchCorrection == "Empirical Bayes") {
      # Remove null values
      combinedExpressionData <- na.omit(combinedExpressionData)
      # Remove batch effect
      combinedExpressionDataBatchRemoved <- ComBat(combinedExpressionData,
                                                   batch = finalBatch)
    } else if (batchCorrection == "Linear Model") {
      # Remove batch effect using Limma
      # (should not be done prior to fitting a linear mode)
      combinedExpressionDataBatchRemoved <- removeBatchEffect(
        combinedExpressionData,
        batch = finalBatch)
    }
  } else {
    combinedExpressionDataBatchRemoved <- combinedExpressionData
  }
  return(combinedExpressionDataBatchRemoved)
}

#' A GEO Function to Convert Two Experiment Information
#' Object into HTML
#'
#' This function allows you to convert two experiment information
#' into HTML
#' @importFrom htmltools HTML
#' @author Guy Hunt
#' @noRd
#' @seealso [extractExperimentInformation()] for GEO object
convertTwoExperimentInformation <- function(experimentData, experimentData2) {
  experimentHtml <- HTML(
    paste(
      "<b>First Experiment</b> <br>",
      experimentData,
      " <br><b>Second Experiment</b><br>",
      experimentData2
    )
  )

  return(experimentHtml)
}

#' A Function to Convert the Column Names to Experimental Information Table
#'
#' This function allows you to Columns to an Experimental Information Table
#' @author Guy Hunt
#' @noRd
convertExpressionDataToExperimentInformation <- function(expressionData) {
  # Convert expression data columns to data frame
  column <- colnames(expressionData)
  # Define other columns
  title <-
    source_name_ch1 <-
    characteristics_ch1 <-
    characteristics_ch1.1 <- NA
  # Create dataframe
  columnInfo <-
    data.frame(column, title, source_name_ch1, characteristics_ch1,
               characteristics_ch1.1)
  # Update rownames
  rownames(columnInfo) <- columnInfo[, 1]
  # Update colnames
  colnames(columnInfo, do.NULL = FALSE)
  colnames(columnInfo) <- c("column", "title", "source_name_ch1",
                            "characteristics_ch1",
                            "characteristics_ch1.1")

  return(columnInfo)
}

#' A Function to combine expression datasets
#'
#' This function allows you to combine expression datasets
#' @author Guy Hunt
#' @noRd
combineExpressionData <- function(expressionData1, expressionData2) {
  # Convert everything to a dataframe
  expressionData1 <- as.data.frame(expressionData1)
  expressionData2 <- as.data.frame(expressionData2)
  
  # Identify the rows in expressionData1 that are not in expressionData2
  expressionData2RowNamesToAdd <- expressionData1[
    !(rownames(expressionData1) %in% rownames(expressionData2)),]
  
  # Identify the rows in expressionData1 that are not in expressionData2
  expressionDataRowNamesToAdd <-
    expressionData2[
      !(rownames(expressionData2) %in% rownames(expressionData1)),]
  
  if (nrow(expressionData2RowNamesToAdd) > 0) {
    # Extract the rownames
    expressionData2RowNamesToAdd <- rownames(expressionData2RowNamesToAdd)
    
    # Convert to a dataframe
    expressionData2RowNamesToAdd <- as.data.frame(expressionData2RowNamesToAdd)
    
    # Update the rownames
    rownames(expressionData2RowNamesToAdd) <- expressionData2RowNamesToAdd[,1]
    
    # Update the colnames
    colnames(expressionData2RowNamesToAdd) <- c("rowname")
    
    # Add the new columns with NA data
    expressionData2RowNamesToAdd[colnames(expressionData2)] <- NA
    
    #x <- data.frame(NA, row.names = row.names(expressionData2RowNamesToAdd))
    
    # Remove the rownames column
    expressionData2RowNamesToAdd <- subset(expressionData2RowNamesToAdd,
                                           select = -rowname)
    
    # Add the new rows to expressionData2
    expressionData2 <- rbind(expressionData2, expressionData2RowNamesToAdd)
  }
  
  if (nrow(expressionDataRowNamesToAdd) > 0) {
    # Extract the rownames
    expressionDataRowNamesToAdd <- rownames(expressionDataRowNamesToAdd)
    
    # Convert to a dataframe
    expressionDataRowNamesToAdd <- as.data.frame(expressionDataRowNamesToAdd)
    
    # Update the rownames
    rownames(expressionDataRowNamesToAdd) <- expressionDataRowNamesToAdd[,1]
    
    # Update the colnames
    colnames(expressionDataRowNamesToAdd) <- c("rowname")
    
    # Add the new columns with NA data
    expressionDataRowNamesToAdd[,colnames(expressionData1)] <- NA
    
    # Remove the rownames column
    expressionDataRowNamesToAdd <- subset(expressionDataRowNamesToAdd,
                                          select = -rowname)
    
    # Add the new rows to expressionData1
    expressionData1 <- rbind(expressionData1, expressionDataRowNamesToAdd)
  }
  
  mergedExpressionData <- cbind(expressionData1, expressionData2)
  
  return(mergedExpressionData)
}

#' A Function to reset all errorChecks
#'
#' This function allows you to reset all errorChecks
#' @author Guy Hunt
#' @noRd
resetErrorChecks <- function(errorChecks) {
  errorChecks$continueWorkflow <- TRUE
  errorChecks$geoAccessionCode <- TRUE
  errorChecks$geoMicroarrayAccessionCode <- TRUE
  errorChecks$geoPlatform <- TRUE
  errorChecks$expressionData <- TRUE
  errorChecks$dataInput <- TRUE
  errorChecks$knnDataInput <- TRUE
  errorChecks$pcaPrcompDataInput <- TRUE
  errorChecks$expressionDataOverTwoColumns <- TRUE
  errorChecks$expressionDataOverOneColumns <- TRUE
  errorChecks$differentialGeneExpression <- TRUE
  errorChecks$differentialGeneExpressionGroup <- TRUE
  errorChecks$uploadFile <- TRUE
  errorChecks$uploadFileExtension <- TRUE
  errorChecks$uploadLogData <- TRUE
  errorChecks$continueWorkflow2 <- TRUE
  errorChecks$geoAccessionCode2 <- TRUE
  errorChecks$geoMicroarrayAccessionCode2 <- TRUE
  errorChecks$geoPlatform2 <- TRUE
  errorChecks$expressionData2 <- TRUE
  errorChecks$dataInput2 <- TRUE
  errorChecks$knnDataInput2 <- TRUE
  errorChecks$pcaPrcompDataInput2 <- TRUE
  errorChecks$expressionDataOverTwoColumns2 <- TRUE
  errorChecks$expressionDataOverOneColumns2 <- TRUE
  errorChecks$differentialGeneExpression2 <- TRUE
  errorChecks$differentialGeneExpressionGroup2 <- TRUE
  errorChecks$uploadFile2 <- TRUE
  errorChecks$uploadFileExtension2 <- TRUE
  errorChecks$uploadLogData2 <- TRUE

  return(errorChecks)

}

#' A Function to Enable users to download dataframes as CSVs
#' @importFrom utils write.csv
#' @import shiny
#' @author Guy Hunt
#' @noRd
dowloadFile <- function(fileName, dataFrame) {
  downloadHandler(
    filename = function() {
      fileName
    },
    content = function(file) {
      write.csv(dataFrame, file,
                row.names = FALSE)
    }
  )
}

#' A Function to remove the lowly expressed genes
#' @importFrom edgeR filterByExpr
#' @author Guy Hunt
#' @noRd
removeLowlyExpressedGenes <- function(ex) {
  keep.exprs <- filterByExpr(rnaExpressionData, group=group)
  ex <- ex[keep.exprs,, keep.lib.sizes=FALSE]

  return(ex)

}

#' A function to download the GEO supplementary files
#'
#' @keywords GEO rnaSeq
#' @importFrom GEOquery getGEOSuppFiles
#' @author Guy Hunt
#' @noRd
downloadGeoSupFiles <- function(geoAccessionCode, directoryPath) {
  geoAccessionCode <- 
  # DOWNLOAD TAR FILE
  geoTarFile <- getGEOSuppFiles(geoAccessionCode, baseDir = directoryPath)
}

#' A function to extract the GEO supplementary files and extract the expression data
#'
#' @keywords GEO rnaSeq
#' @importFrom utils untar
#' @importFrom R.utils gunzip
#' @author Guy Hunt
#' @noRd
extractGeoSupFiles <- function(geoAccessionCode, filePath, directoryPath) {
  tarFileName <- paste0(geoAccessionCode, "_RAW.tar")
  
  setwd(directoryPath)

  # Untar the file
  untar(filePath, exdir = ".")

  # Obtain list of file names in the directory
  listOfFileInDir <- list.files(path = directoryPath)

  # Remove tar file name from list
  listOfFileInDir <- listOfFileInDir[listOfFileInDir != tarFileName]

  # For each file save as a gz file
  for(i in listOfFileInDir){
    try(gunzip(i, overwrite=TRUE))
  }

  return(tarFileName)
}

#' A function to extract the expression data from raw files
#'
#' @author Guy Hunt
#' @importFrom edgeR readDGE
#' @noRd
extractExpressionDataFromGeoSupRawFiles <- function(directoryPath, tarFileName, 
                                                    geneNamesCol, countsCol) {
  setwd(directoryPath)
  
  # Define file names
  files <- list.files(".")
  files <- files[files != tarFileName]

  # Combine the text files into a matrix of counts using edgeR
  expressionData <- tryCatch({
    readDGE(files, columns=c(geneNamesCol,countsCol))},
    error = function(e) {
      # return a safeError if a parsing error occurs
      return(NULL)
    })

  # Set Working directory
  setwd("..")

  return(expressionData)
}

#' A function to update the expression data sample names
#'
#' This function allows you to update the expression data sample names
#' @param rnaExpressionData A object containing the RNA seq expression data
#' @keywords rnaSeq
#' @importFrom stringr str_locate
#' @examples # Define Variables
#' geoAccessionCode <- "GSE63310"
#'
#' # Update Sample Names
#' rnaExpressionData <- calculateSampleNames(rnaExpressionData)
#'
#' # Update Sample Names
#' rnaExpressionData <- calculateSampleNames(rnaExpressionData)
#' @author Guy Hunt
#' @noRd
calculateSampleNames <- function(rnaExpressionData) {
  # Extract the sample information
  samplenames <- substring(colnames(rnaExpressionData), 0,
                           (str_locate(colnames(
                             rnaExpressionData), "_")-1)[,1])
  
  # Update the colnames with the sample names
  colnames(rnaExpressionData) <- samplenames
  
  return(rnaExpressionData)
}

#' Create a directory for GEO supplementary files
#' @author Guy Hunt
#' @noRd
createGeoSupplementaryFilesDirectory <- function() {
  directoryPath <- file.path(getwd(),"geoSupFiles")
  
  if(!dir.exists(directoryPath)){
    try(dir.create(directoryPath, showWarnings = FALSE))
    try(Sys.chmod(directoryPath, "777"))
  }
  
  return(directoryPath)
}

#' Delete directory for GEO supplementary files
#' @author Guy Hunt
#' @noRd
deleteGeoSupplementaryFilesDirectory <- function(directoryPath) {
  unlink(directoryPath, recursive = TRUE, force = TRUE)
}

#' Delete directory for GEO supplementary files
#' @author Guy Hunt
#' @noRd
calculateGeoAccessionDirectory <- function(directoryPath, geoAccessionCode) {
  directoryPath <- file.path(directoryPath,geoAccessionCode)
  return(directoryPath)
}

#' Extract file extensions
#' @importFrom xfun file_ext
#' @author Guy Hunt
#' @noRd
extractFileExtensions <- function(geoTarFile) {
  fileExtensions <- c()
  
  for (file in row.names(geoTarFile)) {
    fileExtensions <-append(fileExtensions, file_ext(file))
  }
  return(fileExtensions)
}

#' Extract expression data from excel
#' @importFrom readxl read_excel
#' @author Guy Hunt
#' @noRd
extractExpressionExcel <- function(filePath) {
  expressionData <- read_excel(filePath)
  
  expressionData <- as.data.frame(expressionData)
  
  return(expressionData)
}

#' Extract expression data from csv
#' @importFrom utils read.csv
#' @author Guy Hunt
#' @noRd
extractExpressionCsv <- function(filePath) {
  expressionData <- read.csv(filePath)
  
  expressionData <- as.data.frame(expressionData)
  
  return(expressionData)
}


#' Deduplicated expression data 
#' @importFrom readxl read_excel
#' @author Guy Hunt
#' @noRd
deduplicatedExpressiondata <- function(expressionData) {
  expressionData <- expressionData[!duplicated(expressionData[,1]),] 
  
  return(expressionData)
}

#' Remove non numeric columns from expression data 
#' @importFrom readxl read_excel
#' @author Guy Hunt
#' @noRd
removeNonNumericColumnsFromExpressiondata <- function(expressionData) {
  for (column in seq_len(ncol(expressionData))) {
    try({expressionData[,column] <- as.numeric(expressionData[,column])})
  }
  
  expressionData <- expressionData[,colSums(is.na(expressionData))<
                                     nrow(expressionData)]  
  return(expressionData)
}

#' Remove White Space from expression data 
#' @importFrom stringr str_trim
#' @author Guy Hunt
#' @noRd
removeWhiteSpace <- function(stringText) {
  stringText <- str_trim(stringText, side = "both")
  return(stringText)
}

#' Converts NAs to 0s in expression data 
#' @author Guy Hunt
#' @noRd
convertNaToZero <- function(expressionData) {
  expressionData[is.na(expressionData)] <- 0
  return(expressionData)
}