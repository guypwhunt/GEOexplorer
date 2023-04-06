#' A Function to create a link to additional study information
#'
#' @keywords GEO
#' @author Guy Hunt
#' @noRd
createStudyLink <- function(website, Id) {
  if (website == "GEO") {
    link <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",Id)
  } else {
    link <- paste0("https://pubmed.ncbi.nlm.nih.gov/",Id)
  }
  linkHtml <- sprintf(
  paste0('<a href=', link, ' target="_blank" class="btn btn-link">',website
         ,'</a>'))

  return(linkHtml)
}


#' A Function to Search the GEO database for studies
#'
#' @keywords GEO
#' @importFrom httr GET
#' @importFrom stringr str_trim str_replace str_replace_all
#' @importFrom xml2 read_xml as_list
#' @importFrom XML xmlParse xmlToList
#' @author Guy Hunt
#' @noRd
searchGeo <- function(searchTerm,
                      firstResultNumber = "0",
                      lastResultNumber = "10") {
  # Trim the search term and replace whitespace
  searchTerm <- str_trim(searchTerm, side = "both")
  searchTerm <-
    str_replace_all(searchTerm, pattern = " ", replacement = "+")
  
  # Define the two URLs
  eSearchUrl <-
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=searchTerm&retmax=lastResultNumber&retstart=firstResultNumber&usehistory=y&tool=GEOexplorer&email=guy.hunt@kcl.ac.uk'
  eSummaryUrl <-
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&version=2.0&retmax=lastResultNumber&retstart=firstResultNumber&usehistory=y&query_key=X&WebEnv=ENTER_WEBENV_PARAMETER_HERE&tool=GEOexplorer&email=guy.hunt@kcl.ac.uk"
  
  
  # Update the first url with the search parameters
  eSearchUrl <- str_replace(eSearchUrl, "searchTerm", searchTerm)
  eSearchUrl <- str_replace(eSearchUrl, "lastResultNumber", "0")
  eSearchUrl <-
    str_replace(eSearchUrl, "firstResultNumber", firstResultNumber)
  
  # Perform the API call
  eSearchResponse <- GET(eSearchUrl)
  
  # Parse the response
  eSearchData <- xmlParse(eSearchResponse)
  eSearchData <- xmlToList(eSearchData)
  
  # Update the second url with the search parameters
  eSummaryUrl <- str_replace(eSummaryUrl, "X", eSearchData$QueryKey)
  eSummaryUrl <-
    str_replace(eSummaryUrl,
                "ENTER_WEBENV_PARAMETER_HERE",
                eSearchData$WebEnv)
  eSummaryUrl <-
    str_replace(eSummaryUrl, "lastResultNumber", lastResultNumber)
  eSummaryUrl <-
    str_replace(eSummaryUrl, "firstResultNumber", firstResultNumber)
  
  # Perform the API call
  eSummaryResponse <- GET(eSummaryUrl)
  
  # Read XML object
  eSummaryXml <- read_xml(eSummaryResponse)
  
  # Convert XML to List
  eSummaryList <- as_list(eSummaryXml)
  
  # Extra the useful list elements
  eSummaryList <- eSummaryList$eSummaryResult$DocumentSummarySet
  eSummaryList <- eSummaryList[2:length(eSummaryList)]
  
  # Define list variables
  geoAccessionCode <- c()
  studyTitle <- c()
  platform <- c()
  pubMedLink <- c()
  geoLink <- c()
  studyOrganism <- c()
  studyPublishedDate <- c()
  microarrayDataSet <- c()
  eSummaryListLength <- length(eSummaryList)
  dataFormat <- c()
  processable <- c()
  
  # For each results update the variables above
  for (result in seq_len(eSummaryListLength)) {
    # Extract the relevant result and the level containing the required
    # information
    eSummaryDataLevel <- eSummaryList[result]$DocumentSummary
    
    # Update all the variables
    geoAccessionCode <- tryCatch({
      append(geoAccessionCode, eSummaryDataLevel$Accession[[1]])
    },
    error = function(cond) {
      append(geoAccessionCode, NA)
    })
    studyTitle <- tryCatch({
      append(studyTitle, eSummaryDataLevel$title[[1]])
    },
    error = function(cond) {
      append(studyTitle, NA)
    })
    platform <- tryCatch({
      append(platform, paste0("GPL", eSummaryDataLevel$GPL) )
    },
    error = function(cond) {
      append(platform, NA)
    })
    pubMedLink <- tryCatch({
      append(pubMedLink, createStudyLink("PubMed",
                                         eSummaryDataLevel$PubMedIds$int[[1]])
      )
    },
    error = function(cond) {
      append(pubMedLink, NA)
    })
    geoLink <- tryCatch({
      append(geoLink, createStudyLink("GEO",
                                      eSummaryDataLevel$Accession[[1]])
      )
    },
    error = function(cond) {
      append(geoLink, NA)
    })
    studyOrganism <- tryCatch({
      append(studyOrganism, eSummaryDataLevel$taxon[[1]])
    },
    error = function(cond) {
      append(studyOrganism, NA)
    })
    studyPublishedDate <- tryCatch({
      append(studyPublishedDate, eSummaryDataLevel$PDAT[[1]])
    },
    error = function(cond) {
      append(studyPublishedDate, NA)
    })
    microarrayDataSet <- tryCatch({
      if (eSummaryDataLevel$GEO2R[[1]] == "yes") {
        append(microarrayDataSet, "Microarray")
      } else {
        append(microarrayDataSet, "RNA Sequencing")
      }
    },
    error = function(cond) {
      append(microarrayDataSet, NA)
    })
    dataFormat <- tryCatch({
      append(dataFormat, eSummaryDataLevel$suppFile[[1]])
    },
    error = function(cond) {
      append(dataFormat, NA)
    })
    processable <- tryCatch({
      if (eSummaryDataLevel$GEO2R[[1]] == "yes") {
        append(processable, "Yes")
      } else if (
        grepl("XLSX", eSummaryDataLevel$suppFile[[1]], fixed = TRUE) |
        grepl("TXT", eSummaryDataLevel$suppFile[[1]], fixed = TRUE) |
        grepl("CSV", eSummaryDataLevel$suppFile[[1]], fixed = TRUE) |
        grepl("TSV", eSummaryDataLevel$suppFile[[1]], fixed = TRUE) |
        grepl("TAR", eSummaryDataLevel$suppFile[[1]], fixed = TRUE) |
        grepl("XLS", eSummaryDataLevel$suppFile[[1]], fixed = TRUE)) {
        append(processable, "Potentially, if not, please download the file(s) 
               from the GEO link and format them into a count matrix as per 
               the template in the 'Example Datasets' tab")
      } else {
        append(processable, "If the file(s) are downloaded from the GEO link 
               and formatted into a count matrix as per the template in 
               the 'Example Datasets' tab")
      }
    },
    error = function(cond) {
      append(processable, NA)
    })
  }
  
  # Define the search results table
  searchResultsTable <- matrix(ncol = 11, nrow = eSummaryListLength)
  
  # Update the search results table column names
  colnames(searchResultsTable) <- c(
    "GEO Accession Code",
    "Title",
    "Organism",
    "Published Date",
    "Platform",
    "Data Type",
    "Data Format",
    "Processable",
    "GEO Link",
    "PubMed Link",
    "Load as First Dataset"
    #, "Load as Second Dataset"
  )
  
  # Update the results table with values
  searchResultsTable[, "GEO Accession Code"] <- geoAccessionCode
  searchResultsTable[, "Title"] <- studyTitle
  searchResultsTable[, "Organism"] <- studyOrganism
  searchResultsTable[, "Published Date"] <- studyPublishedDate
  searchResultsTable[, "Platform"] <- platform
  searchResultsTable[, "GEO Link"] <- geoLink
  searchResultsTable[, "PubMed Link"] <- pubMedLink
  searchResultsTable[, "Data Type"] <- microarrayDataSet
  searchResultsTable[, "Data Format"] <- dataFormat
  searchResultsTable[, "Processable"] <- processable
  
  
  # Add action buttons
  searchResultsTable[, "Load as First Dataset"] <-
    addActionButtonToTable(actionButton, eSummaryListLength, "button_",
                           label = "Load",
                           onclick =
                             'Shiny.onInputChange(\"loadGeoSearchAsFirstDataset\",
                           this.id)' )
  
  geoSearchResults <- NULL
  geoSearchResults$searchResultsTable <- searchResultsTable
  geoSearchResults$totalResults <- eSearchData$Count
  
  return(geoSearchResults)
}

#' A Function to add action buttons to a table
#'
#' @keywords GEO
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @author Guy Hunt
#' @noRd
addActionButtonToTable <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))
  }
  return(inputs)
}
