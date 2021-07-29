library(GEOexplorer)
context("extractPlatforms")

test_that("Microarray GSE with one platform is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  expect_type(platforms, 'character')
  expect_type(platform, 'character')
  expect_equal(platform,"GPL6246")
})
