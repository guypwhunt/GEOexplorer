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

test_that("Microarray GSE with multiple platforms is handled correctly", {
  geoAccessionCode <- "GSE1838"
  allGset <- getGeoObject(geoAccessionCode)

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  expect_type(platforms, 'character')
  expect_type(platform, 'character')
  expect_equal(platform,"GPL1444")
  expect_equal(platforms[2],"GPL1445")
  expect_equal(platforms[3],"GPL1446")
})
