library(GEOexplorer)
context("extractExpressionData")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  expect_type(expressionData, 'double')
  expect_equal(ncol(expressionData),8)
  expect_equal(nrow(expressionData),35557)
})
