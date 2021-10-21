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
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data
  columnNames <-
    c("title",
      "source_name_ch1",
      "characteristics_ch1",
      "characteristics_ch1.1")
  finalColumnNames <- c()
  i <- 1

  for (name in columnNames) {
    if (name %in% colnames(phenoData)) {
      finalColumnNames <- c(finalColumnNames, name)
    }
  }

  df <- data.frame(column = row.names(phenoData))

  for (name in finalColumnNames) {
    df <- data.frame(df, phenoData[name])
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

#' A function to obtain the supplementary files and extract the expression data
#'
#' This function allows you to source a the GEO supplementary files
#' and extract the expression data from the files
#' @param geoAccessionCode A character string representing
#' a GEO object for download and parsing
#' @keywords GEO rnaSeq
#' @importFrom GEOquery getGEOSuppFiles
#' @importFrom utils untar
#' @importFrom R.utils gunzip
#' @examples # Define Variables
#' geoAccessionCode <- "GSE63310"
#'
#' # Source and extract expression data
#' rnaExpressionData <- extractGeoSupFiles(geoAccessionCode)
#' @author Guy Hunt
#' @noRd
extractGeoSupFiles <- function(geoAccessionCode, geneNamesCol, countsCol) {
  # Define Variables
  path <- paste0(".\\", geoAccessionCode)
  tarFileName <- paste0(geoAccessionCode, "_RAW.tar")

  # DOWNLOAD TAR FILE
  geoTarFile <- getGEOSuppFiles(geoAccessionCode)

  # Set working directory
  setwd(path)

  # Untar the file (THIS IS ERORRING) (CONTINUE HERE)
  untar(tarFileName, exdir = ".")

  # Obtain list of file names in the directory
  listOfFileInDir <- list.files(path = path)

  # Remove tar file name from list
  listOfFileInDir <- listOfFileInDir[listOfFileInDir != tarFileName]

  # For each file save as a gz file
  for(i in listOfFileInDir)
    gunzip(i, overwrite=TRUE)

  # Define file names
  files <- list.files(".")
  files <- files[files != tarFileName]

  # Combine the text files into a matrix of counts using edgeR
  expressionData <- readDGE(files, columns=c(geneNamesCol,countsCol))

  # Delete files in directory
  for (fileName in list.files(".")){
    file.remove(fileName)
  }

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
  samplenames

  # Update the colnames with the sample names
  colnames(rnaExpressionData) <- samplenames

  return(rnaExpressionData)
}

#' A function to convert expression data to counts per million
#'
#' This function allows you to convert expression data to counts per million
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
#'
#' # Raw counts are converted to counts-per-million (CPM)
#' cpm <- cpm(rnaExpressionData)
#' @author Guy Hunt
#' @noRd
calculateCountsPerMillion <- function(rnaExpressionData) {
  # Raw counts are converted to counts-per-million (CPM)
  cpm <- cpm(rnaExpressionData)

  return(cpm)
}
