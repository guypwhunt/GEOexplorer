library(GEOexplorer)
context("calculateLogTransformation")

test_that("Microarray GSE is handled correctly with Auto-Detect log transformation", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  logTransformation <- "Auto-Detect"  #
  dataInput <- calculateLogTransformation(expressionData, logTransformation)

  expect_type(dataInput, 'double')
  expect_equal(ncol(dataInput),8)
  expect_equal(nrow(dataInput),35557)
  expect_equal(dataInput[1,1],11.711505)
})
