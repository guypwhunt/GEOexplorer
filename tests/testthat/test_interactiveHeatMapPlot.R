library(GEOexplorer)
context("interactiveHeatMapPlot")

test_that("Microarray GSE without missing values is handled correctly by interactiveHeatMapPlot", {
  geoAccessionCode <- 'GSE18388'
  allGset <- getGeoObject(geoAccessionCode)

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  logTransformation <- "Auto-Detect"
  dataInput <- calculateLogTransformation(expressionData, logTransformation)

  knnDataInput <- calculateKnnImpute(dataInput, "No")

  naOmitInput <- calculateNaOmit(knnDataInput)

  fig <- interactiveHeatMapPlot(naOmitInput)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
})

test_that("Microarray GSE with missing values  is handled correctly by interactiveHeatMapPlot", {
  geoAccessionCode <- 'GSE18380'
  allGset <- getGeoObject(geoAccessionCode)

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  logTransformation <- "Auto-Detect"
  dataInput <- calculateLogTransformation(expressionData, logTransformation)

  knnDataInput <- calculateKnnImpute(dataInput, "No")

  naOmitInput <- calculateNaOmit(knnDataInput)

  fig <- interactiveHeatMapPlot(naOmitInput)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
})
