library(GEOexplorer)
context("extractSampleNames")

test_that("Microarray GSE is handled correctly for expression data", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  columns <- extractSampleNames(expressionData)

  expect_type(columns, 'character')
  expect_equal(columns[1],"GSM458594")
})
