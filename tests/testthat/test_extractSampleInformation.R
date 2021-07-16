library(GEOexplorer)
context("extractSampleInformation")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  sampleInfo <- extractSampleInformation(gsetData)

  expect_type(sampleInfo, 'list')
  expect_equal(nrow(sampleInfo),8)
  expect_equal(ncol(sampleInfo),35)
})
