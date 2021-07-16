library(GEOexplorer)
context("nonInteractiveUmapPlot")

test_that("Microarray GSE without missing values is handled correctly with nonInteractiveUmapPlot", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  fig <- nonInteractiveUmapPlot(naOmitInput, 2, geoAccessionCode)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$x, 'double')
  expect_type(fig$y, 'double')
})

test_that("Microarray GSE with missing values is handled correctly with nonInteractiveUmapPlot", {
  allGset <- getGeoObject('GSE18380')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  fig <- nonInteractiveUmapPlot(naOmitInput, 2, geoAccessionCode)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$x, 'double')
  expect_type(fig$y, 'double')
})
