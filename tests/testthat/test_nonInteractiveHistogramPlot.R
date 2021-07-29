library(GEOexplorer)
context("nonInteractiveHistogramPlot")

test_that("Microarray GSE without missing values is handled correctly with nonInteractiveHistogramPlot", {
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

  fig <- nonInteractiveHistogramPlot(fit2, adjustment)

  expect_type(fig, "list")
  expect_type(fig$breaks, "double")
  expect_type(fig$counts, "integer")
  expect_type(fig$density, "double")
  expect_type(fig$mids, "double")
  expect_type(fig$xname, "character")
  expect_type(fig$equidist, "logical")
})
