library(GEOexplorer)
context("extractSampleDetails")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  columnInfo <- extractSampleDetails(gsetData)

  expect_type(columnInfo, 'list')
  expect_equal(ncol(columnInfo),5)
  expect_equal(nrow(columnInfo),8)
})
