#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object from GEO when you know the specific Geo Accession code and platform
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param platform The study's platform
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples extractGeoData("GSE18380", "GPL4694")
#' @author Guy Hunt
extractGeoData <- function(geoAccessionCode, platform) {
  library(GEOquery)
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Data Sourcing Function
#'
#' This function allows you to source a GEO Object from GEO when you only know the Geo Accession code
#' @param geoAccessionCode A character string representing a GEO object for download and parsing
#' @param GSEMatrix A boolean telling GEOquery whether or not to use GSE Series Matrix files from GEO. The parsing of these files can be many orders-of-magnitude faster than parsing the GSE SOFT format files. Defaults to TRUE, meaning that the SOFT format parsing will not occur; set to FALSE if you for some reason need other columns from the GSE records.
#' @param getGPL A boolean defaulting to TRUE as to whether or not to download and include GPL information when getting a GSEMatrix file. You may want to set this to FALSE if you know that you are going to annotate your featureData using Bioconductor tools rather than relying on information provided through NCBI GEO. Download times can also be greatly reduced by specifying FALSE.
#' @param platformAnnotation A string defaulting to "NCBI generated" meaning true as to whether or not to use the Annotation GPL information. These files are nice to use because they contain up-to-date information remapped from Entrez Gene on a regular basis. However, they do not exist for all GPLs; in general, they are only available for GPLs referenced by a GDS. Input "Submitter supplied" for FALSE
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples allGset <- getGeoObject("GSE18380", GSEMatrix=TRUE, getGPL=TRUE, platformAnnotation = "NCBI generated")
#' @author Guy Hunt
getGeoObject <- function(geoAccessionCode, GSEMatrix=TRUE, getGPL=TRUE, platformAnnotation = "NCBI generated") {
  library(GEOquery)
  if (platformAnnotation == "Submitter supplied") {
    platformAnnotation <- FALSE
  } else if (platformAnnotation == "NCBI generated") {
    platformAnnotation <- TRUE
  } else {
    platformAnnotation <- TRUE
  }
  gset <- getGEO(geoAccessionCode, GSEMatrix=GSEMatrix, AnnotGPL=platformAnnotation) # getGPL=getGPL,
  return(gset)
}

#' A GEO Function to Obtain the Available Platforms of a GEO Object
#'
#' This function allows you to extract the platforms codes from a GEO object
#' @param gset The GEO object which can be obtained from the getGeoObject() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples platforms <- extractPlatforms(allGset)
#' @author Guy Hunt
#' @seealso [getGeoObject()] for GEO object
extractPlatforms <- function(gset) {
  library(GEOquery)
  platforms <- c()
  i <-1
  for(dataset in gset) {
    platforms[i] <- annotation(dataset)
    i <- i + 1
  }
  return(platforms)
}

#' A GEO Function to Obtain the GEO Object for a specific platform
#'
#' This function allows you to extract a platform's GEO object from a GEO object
#' @param gset The GEO object which can be obtained from the getGeoObject() function
#' @param platform The platform code
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples gsetData <- extractPlatformGset(gset, platforms[1])
#' @author Guy Hunt
#' @seealso [getGeoObject()] for GEO object
extractPlatformGset <- function(gset, platform) {
  library(GEOquery)
  if (length(gset) > 1) idx <- grep(platform[1], attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

#' A GEO Function to Obtain the Experiment Information from a GEO Object
#'
#' This function allows you to extract experiment information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples experimentInformation <- extractExperimentInformation(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractExperimentInformation <- function(gset) {
  library(GEOquery)
  experimentalData <- experimentData(gset)
  return(experimentalData)
}

#' A GEO Function to Convert the Experiment Information Object into HTML
#'
#' This function allows you to convert experiment information into HTML
#' @param experimentData The experiment object obtained from the extractExperimentInformation() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples convertExperimentInformation <- convertExperimentInformation(experimentInformation)
#' @author Guy Hunt
#' @seealso [extractExperimentInformation()] for GEO object
convertExperimentInformation <- function(experimentData) {
  library(GEOquery)
  name <- paste("<b>", "Author's Name:", "</b>", "<p>", experimentData@name, "</p>")
  lab <- paste("<b>", "Laboratory:", "</b>", "<p>", experimentData@lab, "</p>")
  contact <- paste("<b>", "Contact Information:", "</b>", "<p>", experimentData@contact, "</p>")
  title <- paste("<b>", "Paper Title:", "</b>", "<p>", experimentData@title, "</p>")
  abstract <- paste("<b>", "Abstract:", "</b>", "<p>", experimentData@abstract, "</p>")
  url <- paste("<b>", "Paper URL:", "</b>", "<p>", experimentData@url, "</p>")
  pubMedIds <- paste("<b>", "PubMed ID:", "</b>", "<p>", experimentData@pubMedIds, "</p>")
  html <- HTML(paste(title, name, lab, contact, url, pubMedIds, abstract, sep = ""))
}

#' A GEO Function to Extract Information on the Samples from a GEO object
#'
#' This function allows you to extract the studies sample information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples columnInfo <- extractSampleDetails(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractSampleDetails <- function(gset){
  library(GEOquery)
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data
  columnNames <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1")
  finalColumnNames <- c()
  i <- 1

  for (name in columnNames) {
    if (name %in% colnames(phenoData)) {
      finalColumnNames <- c(finalColumnNames,name)
    }
  }

  df <- data.frame(column=row.names(phenoData))

  for (name in finalColumnNames){
    df <- data.frame(df, phenoData[name])
  }
  return(df)
}

#' A GEO Function to Extract Expression Data from a GEO object
#'
#' This function allows you to extract the studies expression object from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @import GEOquery
#' @examples expressionData <- extractExpressionData(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractExpressionData <- function(gset) {
  library(GEOquery)
  ex <- exprs(gset)
  return(ex)}

#' A GEO Function to Extract Information on the Samples from a GEO object
#'
#' This function allows you to extract the studies sample information from a GEO object
#' @param gset The GEO object which can be obtained from the extractPlatformGset() function
#' @keywords GEO
#' @export
#' @examples sampleInfo <- extractSampleInformation(gsetData)
#' @author Guy Hunt
#' @seealso [extractPlatformGset()] for GEO object
extractSampleInformation <- function(gset) {
  library(GEOquery)
  sampleInfo <- pData(gset)
  return(sampleInfo)}
