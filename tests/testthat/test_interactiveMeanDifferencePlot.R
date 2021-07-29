library(GEOexplorer)
context("interactiveMeanDifferencePlot")

test_that("Microarray GSE without missing values is handled correctly with interactiveMeanDifferencePlot", {
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

  ct <- 1

  fig <- interactiveMeanDifferencePlot(fit2, dT, ct)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
})
