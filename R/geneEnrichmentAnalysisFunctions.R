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

#' A Function to calculate the -log(P value)
#' @author Guy Hunt
#' @noRd
calculateLogPValue  <- function(geneEnrichmentTable) {

  geneEnrichmentTable[,"Minus.Log.P.Value"] <-
    0-log10(geneEnrichmentTable[, "P.value"])

  return(geneEnrichmentTable)
}

#' A Function to calculate the overlap fraction values
#' @importFrom stringr str_split
#' @author Guy Hunt
#' @noRd
calculateOverlapFractions  <- function(geneEnrichmentTable) {

  # Update this to divide each number
  overlapFractions <- str_split(geneEnrichmentTable[,"Overlap"],"/")

  fractionList <- c()

  for (fractionValue in overlapFractions) {
    try(numerator <- as.numeric(fractionValue[1]))
    try(denominator <- as.numeric(fractionValue[2]))
    try(fraction <- numerator/denominator)

    fractionList <- tryCatch({
      append(fractionList, fraction)
    }, error = function(e) {
      # return a safeError if a parsing error occurs
      append(fractionList, NULL)
    })
  }

  geneEnrichmentTable[,"Overlap.Value"] <- fractionList

  return(geneEnrichmentTable)
}

#' A Function to create an interactive gene enrichment Manhattan Plot
#' @import plotly
#' @author Guy Hunt
#' @noRd
interactiveGeneEnrichmentManhattanPlot  <- function(geneEnrichmentTable,
                                                    columnToSort) {

  selectedColumn <- geneEnrichmentColumnSelection(geneEnrichmentTable,
                                                  columnToSort)

  fig19 <- plot_ly(geneEnrichmentTable,
                 y = selectedColumn,
                 x = ~Term,
                 color = selectedColumn,
                 type = "scatter",
                 mode = 'markers',
                 text = geneEnrichmentPlotText(geneEnrichmentTable)
  )

  fig19 <- fig19 %>% layout(xaxis = list(categoryorder = "array",
                                     categoryarray = selectedColumn,
                                     showticklabels = FALSE,
                                     title = ""
                                     ),
                        yaxis = list(title = columnToSort))
  
  try(fig19 <- toWebGL(fig19))
  
  return(fig19)
}

#' A Function to create an interactive gene enrichment volcano plot
#' @import plotly
#' @author Guy Hunt
#' @noRd
interactiveGeneEnrichmentVolcanoPlot  <- function(geneEnrichmentTable) {
  fig20 <- plot_ly(
    data = geneEnrichmentTable,
    x = ~ Odds.Ratio,
    y = ~ Minus.Log.P.Value,
    color = ~ Combined.Score,
    text = geneEnrichmentPlotText(geneEnrichmentTable),
    type = 'scatter',
    mode = 'markers')

  fig20 <- fig20 %>% layout(
    xaxis = list(title = "Odds Ratio"),
    yaxis = list(title = "-log10(P-value)")
  )
  
  try(fig20 <- toWebGL(fig20))
  
  return(fig20)
}

#' A Function to create an interactive gene enrichment bar plot
#' @import plotly
#' @author Guy Hunt
#' @noRd
interactiveGeneEnrichmentBarPlot  <- function(geneEnrichmentTable,
                                              columnToSort) {
  selectedColumn <- geneEnrichmentColumnSelection(geneEnrichmentTable,
                                                  columnToSort)

  fig21 <- plot_ly(geneEnrichmentTable,
                 x = selectedColumn,
                 color = selectedColumn,
                 text = geneEnrichmentPlotText(geneEnrichmentTable),
                 y = ~Term, type = 'bar', orientation = 'h')

  fig21 <- fig21 %>% layout(xaxis = list(title = columnToSort),
                        yaxis = list(title = "Term",
                                     categoryorder = "array",
                                     categoryarray = selectedColumn))
  
  return(fig21)
}

#' A Function to return the gene enrichment column selection logic
#' @import plotly
#' @author Guy Hunt
#' @noRd
geneEnrichmentColumnSelection <- function(geneEnrichmentTable, columnToSort) {
  columnSelection <- if(columnToSort == "Odds.Ratio") {
    ~Odds.Ratio} else
    if(columnToSort == "P.value") {
      ~P.value} else
      if(columnToSort == "Adjusted.P.value") {
        ~Adjusted.P.value} else
        if(columnToSort == "Odds.Ratio") {
          ~Odds.Ratio} else
          if(columnToSort == "Combined.Score") {
            ~Combined.Score} else
            if(columnToSort == "Minus.Log.P.Value") {
              ~Minus.Log.P.Value} else
              if(columnToSort == "Overlap.Value") {
                ~Overlap.Value}
  return(columnSelection)
}

#' A Function to sort the gene enrichment text for plots
#' @import plotly
#' @author Guy Hunt
#' @noRd
geneEnrichmentPlotText <- function(geneEnrichmentTable) {
  plotText <- {~ paste(
    'Term: ',
    Term,
    '<br>',
    'Genes: ',
    Genes ,
    '<br>',
    'Overlap: ',
    Overlap ,
    '<br>',
    'P.value: ',
    P.value ,
    '<br>',
    'Adjusted P-value: ',
    Adjusted.P.value ,
    '<br>',
    'Odds Ratio: ',
    Odds.Ratio,
    '<br>',
    'Combined.Score: ',
    Combined.Score ,
    '<br></br>'
  )
    }

    return(plotText)
}


#' A Function to gene enrichment table
#' @author Guy Hunt
#' @noRd
sortGeneEnrichmentTable <- function(geneEnrichmentTable,
                                    columnToSort = "Adjusted.P.value",
                                    sortDecreasingly = FALSE) {

  geneEnrichmentTable <-
    geneEnrichmentTable[order(geneEnrichmentTable[,columnToSort], decreasing
                              = sortDecreasingly),]

  return(geneEnrichmentTable)
}

#' A Function to select the top gene enrichment records
#' @author Guy Hunt
#' @noRd
selectTopGeneEnrichmentRecords <- function(geneEnrichmentTable,
                                           recordsToDisplay = 10) {

  geneEnrichmentTable <- geneEnrichmentTable[seq_len(recordsToDisplay),]


  return(geneEnrichmentTable)
}


#' A Function to convert ascending/descebding to false/true
#' @author Guy Hunt
#' @noRd
convertUiSortingMethod <- function(sortDecreasingly) {
  if (sortDecreasingly == "Ascendingly") {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

