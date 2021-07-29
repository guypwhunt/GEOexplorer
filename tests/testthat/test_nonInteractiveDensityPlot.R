library(GEOexplorer)
context("nonInteractiveDensityPlot")

test_that("Microarray GSE without missing values is handled correctly with nonInteractiveDensityPlot", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  fig <- nonInteractiveDensityPlot(ex = naOmitInput, geoAccessionCode = 'GSE18388', platform = platform)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$X, 'double')
  expect_type(fig$Y, 'double')
})
