library(GEOquery)

getGeoData <- function(geoAccessionCode, platform) {
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

getGset <- function(geoAccessionCode, platformAnnotation) {
  if (platformAnnotation == "Submitter supplied") {
    platformAnnotation == FALSE
  } else if (platformAnnotation == "NCBI generated") {
    platformAnnotation == TRUE
  } 
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE, AnnotGPL=platformAnnotation)
  return(gset)
}

getPlatforms <- function(gset) {
  platforms <- list()
  i <-1
  for(dataset in gset) {
    platforms[[i]] <- annotation(dataset)
    i <- i + 1
  }
  return(platforms)
  }

getPlatformGset <- function(gset, platform) {
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

getExperimentInformation <- function(gset) {
  experimentalData <- experimentData(gset)
  return(experimentalData)
}

extractExperimentInformation <- function(experimentData) {
  name <- paste("<b>", "Author's Name:", "</b>", "<p>", experimentData@name, "</p>")
  lab <- paste("<b>", "Laboratory:", "</b>", "<p>", experimentData@lab, "</p>")
  contact <- paste("<b>", "Contact Information:", "</b>", "<p>", experimentData@contact, "</p>")
  title <- paste("<b>", "Paper Title:", "</b>", "<p>", experimentData@title, "</p>")
  abstract <- paste("<b>", "Abstract:", "</b>", "<p>", experimentData@abstract, "</p>")
  url <- paste("<b>", "Paper URL:", "</b>", "<p>", experimentData@url, "</p>")
  pubMedIds <- paste("<b>", "PubMed ID:", "</b>", "<p>", experimentData@pubMedIds, "</p>")
  html <- HTML(paste(title, name, lab, contact, url, pubMedIds, abstract, sep = ""))
}

getColumnDetails <- function(gset){
  phenoDataset <- phenoData(gset)
  phenoData <- phenoDataset@data
  df <- data.frame(
    column=row.names(phenoData),
    title=phenoData["title"],
    source=phenoData["source_name_ch1"], 
    characteristic1=phenoData["characteristics_ch1"], 
    characteristic2=phenoData["characteristics_ch1.1"]) 
  return(df)
}