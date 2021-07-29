library(GEOexplorer)
context("interactiveMeanVariancePlot")

test_that("Microarray GSE without missing values is handled correctly by the interactive Mean Variance plot", {
  geoAccessionCode <- 'GSE18388'
  allGset <- getGeoObject(geoAccessionCode)

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  logTransformation <- "Auto-Detect"
  dataInput <- calculateLogTransformation(expressionData, logTransformation)

  knnDataInput <- calculateKnnImpute(dataInput, "Yes")

  naOmitInput <- calculateNaOmit(knnDataInput)

  fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode, gsetData)

  expect_type(fig, 'list')
  expect_type(fig[1], 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
})
