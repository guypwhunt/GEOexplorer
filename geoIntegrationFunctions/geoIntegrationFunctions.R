library(GEOquery)

getGeoData <- function(geoAccessionCode, platform) {
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  return(gset)
}

getGset <- function(geoAccessionCode) {
  gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
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
  experimentData <- experimentData(gset)
  return(experimentData)
}

extractExperimentInformation <- function(experimentData) {
  name <- paste("<b>", "Experimenter's Name:", "</b>", experimentData@name)
  lab <- paste("<b>", "Laboratory:", "</b>", experimentData@lab)
  contact <- paste("<b>", "Contact Information:", "</b>", experimentData@contact)
  title <- paste("<b>", "Paper Title:", "</b>", experimentData@title)
  abstract <- paste("<b>", "Abstract:", "</b>", experimentData@abstract)
  url <- paste("<b>", "Paper URL:", "</b>", experimentData@url)
  pubMedIds <- paste("<b>", "PubMed ID:", "</b>", experimentData@pubMedIds)
  html <- HTML(paste(title, name, lab, contact, url, pubMedIds, abstract, sep = "<br/><br/>"))
}