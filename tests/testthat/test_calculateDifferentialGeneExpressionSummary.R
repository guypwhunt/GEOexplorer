library(GEOexplorer)
context("calculateDifferentialGeneExpressionSummary")

test_that("Microarray GSE without missing values is handled correctly with calculateDifferentialGeneExpressionSummary", {
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

  significanceLevelCutOff <- 0.05
  dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)

  expect_type(dT, 'double')
  expect_equal(ncol(dT),1)
  expect_equal(nrow(dT),35557)
})

test_that("Microarray GSE with missing values is handled correctly with calculateDifferentialGeneExpressionSummary", {
  allGset <- getGeoObject('GSE18380')

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

  significanceLevelCutOff <- 0.05
  dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)

  expect_type(dT, 'double')
  expect_equal(ncol(dT),1)
  expect_equal(nrow(dT),3234)
})

