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
  #studySummary <- c()
  studyOrganism <- c()
  studyPublishedDate <- c()
  microarrayDataSet <- c()
  eSummaryListLength <- length(eSummaryList)

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
      return(NA)
    })

    studyTitle <- tryCatch({
      append(studyTitle, eSummaryDataLevel$title[[1]])
    },
    error = function(cond) {
      return(NA)
    })
    #studySummary <- append(studySummary, eSummaryDataLevel$summary[[1]])
    studyOrganism <- tryCatch({
      append(studyOrganism, eSummaryDataLevel$taxon[[1]])
    },
    error = function(cond) {
      return(NA)
    })
    studyPublishedDate <- tryCatch({
      append(studyPublishedDate, eSummaryDataLevel$PDAT[[1]])
    },
    error = function(cond) {
      return(NA)
    })
    microarrayDataSet <- tryCatch({
      append(microarrayDataSet, eSummaryDataLevel$GEO2R[[1]])
    },
    error = function(cond) {
      return(NA)
    })
  }

  # Define the search results table
  searchResultsTable <- matrix(ncol = 6, nrow = eSummaryListLength)

  # Update the search results table column names
  colnames(searchResultsTable) <- c(
    "GEO Accession Code",
    "Title",
    "Organism",
    "Published Date",
    #"studySummary",
    "Microarray DataSet",
    "Actions"
  )

  # Update the rtesults table with values
  searchResultsTable[, "GEO Accession Code"] <- geoAccessionCode
  searchResultsTable[, "Title"] <- studyTitle
  searchResultsTable[, "Organism"] <- studyOrganism
  searchResultsTable[, "Published Date"] <- studyPublishedDate
  #searchResultsTable[,"studySummary"] <- studySummary
  searchResultsTable[, "Microarray DataSet"] <- microarrayDataSet

  # Add action buttons
  searchResultsTable[, "Actions"] <-
    addActionButtonToTable(actionButton, eSummaryListLength, "button_",
                           label = "Load Dataset",
                           onclick =
                             'Shiny.onInputChange(\"loadGeoSearchAccession\",
                           this.id)' )

  return(searchResultsTable)
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
