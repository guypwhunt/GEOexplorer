library(GEOexplorer)
context("interactiveThreeDDensityPlot")

test_that("Microarray GSE without missing values is handled correctly by the 3D interactive density plot", {
  geoAccessionCode <- 'GSE18388'
  allGset <- getGeoObject(geoAccessionCode)

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  fig <- interactiveThreeDDensityPlot(naOmitInput, geoAccessionCode, platform)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
})
