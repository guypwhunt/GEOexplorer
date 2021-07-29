library(GEOexplorer)
context("calculateTopDifferentiallyExpressedGenes")

test_that("Microarray GSE without missing values and NCBI submitted column names is handled correctly with calculateTopDifferentiallyExpressedGenes", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  columnNames <- extractSampleNames(expressionData)

  numberOfColumns <- length(columnNames)
  numberOfColumns <- numberOfColumns + 1
  halfNumberOfColumns <- ceiling(numberOfColumns/2)
  i <- 0

  group1 <- c()
  group2 <- c()

  for (name in columnNames) {
    if (i < halfNumberOfColumns) {
      group1 <- c(group1, name)
      i <- i +1
    } else {
      group2 <- c(group2, name)
      i <- i +1
    }
  }

  column2 <- calculateExclusiveColumns(columnNames, group1)

  gsms <- calculateEachGroupsSamples(columnNames,group1, group2)

  limmaPrecisionWeights <- "Yes"
  forceNormalization <- "Yes"
  fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, expressionData)

  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
  adjustment <- convertAdjustment(pValueAdjustment)
  tT <- calculateTopDifferentiallyExpressedGenes(fit2, adjustment)

  expect_type(tT, "list")
  expect_type(tT$ID, "character")
  expect_type(tT$t, "double")
  expect_type(tT$Gene.symbol, "character")
  expect_type(tT$adj.P.Val, "double")
  expect_type(tT$B, "double")
  expect_type(tT$Gene.title, "character")
  expect_type(tT$P.Value, "double")
  expect_type(tT$logFC, "double")
  expect_type(tT$Gene.ID, "character")
})
