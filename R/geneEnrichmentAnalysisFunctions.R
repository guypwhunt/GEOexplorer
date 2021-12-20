#' A Function to extract the differentially expressed genes
#' @author Guy Hunt
#' @noRd
extractDifferenitallyExpressedGenes <- function(tT, dT) {
  # Merge the top differentially expressed genes with of they are upregulated
  # or not
  differentiallyExpressedGenes <- merge(tT, dT, by = "row.names")

  # Add rownames
  rownames(differentiallyExpressedGenes) <-
    differentiallyExpressedGenes$Row.names

  # Remove non-differentially expressed genes
  differentiallyExpressedGenes <-
    differentiallyExpressedGenes[!(
      differentiallyExpressedGenes[, "Group1-Group2"] ==
                                     0),]

  return(differentiallyExpressedGenes)
}

#' A Function to extract the upregulated genes
#' @author Guy Hunt
#' @noRd
extractUpregulatedGenes <- function(differentiallyExpressedGenes) {
  # Remove downregulated genes
  upregulatedGenes <-
    differentiallyExpressedGenes[!(
      differentiallyExpressedGenes[, "Group1-Group2"] ==
                                     -1),]

  return(upregulatedGenes)
}

#' A Function to extract the downregulated genes
#' @author Guy Hunt
#' @noRd
extractdowregulatedGenes <- function(differentiallyExpressedGenes) {
  # Remove upregulated genes
  downregulatedGenes <-
    differentiallyExpressedGenes[!(
      differentiallyExpressedGenes[, "Group1-Group2"] ==
                                     1),]

  return(downregulatedGenes)
}

#' A Function to extract the downregulated genes
#' @author Guy Hunt
#' @noRd
extractGeneSymbols <- function(differentiallyExpressedGenes,
                   geneSymbolColumnName = "Gene.symbol") {
  # Extract all differentially expressed gene symbols
  geneSymbols <-
    differentiallyExpressedGenes[, geneSymbolColumnName]

  return(geneSymbols)
}

#' A Function to extract the available databases from enrich R
#' @importFrom enrichR listEnrichrDbs
#' @author Guy Hunt
#' @noRd
extractDatabaseNamesFromEnrichR <- function() {

  databases <- listEnrichrDbs()
  databaseNames <- databases$libraryName
  databaseNames <- as.list(databaseNames)
  return(databaseNames)
}

#' A Function to enrich the gene using enrich R
#' @importFrom enrichR enrichr
#' @author Guy Hunt
#' @noRd
enrichGenes <- function(geneSymbols,
             databaseNames =
               "GO_Biological_Process_2015") {

  # Query the enrich R database
  enrichedGeneInformation <- enrichr(geneSymbols, databaseNames)

  # Extract the relevant information
  enrichedGeneInformation <- enrichedGeneInformation[[1]]
  return(enrichedGeneInformation)
}

#' A Function to plot the gene enrichment information
#' @importFrom enrichR plotEnrich
#' @author Guy Hunt
#' @noRd
plotGeneEnrichmentinformation  <- function(geneEnrichmentTable) {

  # Query the enrich R database
  plot <- plotEnrich(geneEnrichmentTable, showTerms = 20, numChar = 40,
             y = "Count", orderBy = "P.value")


  return(plot)
}
