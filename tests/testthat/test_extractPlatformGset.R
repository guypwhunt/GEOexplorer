library(GEOexplorer)
context("extractPlatformGset")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expect_type(gsetData, 'S4')
  expect_s4_class(gsetData, 'ExpressionSet')
  expect_equal(nrow(pData(gsetData)),8)
  expect_equal(nrow(fData(gsetData)),35557)
})
