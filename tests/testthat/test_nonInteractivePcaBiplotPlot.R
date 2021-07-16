library(GEOexplorer)
context("nonInteractivePcaBiplotPlot")

test_that("Microarray GSE without missing values is handled correctly with nonInteractivePcaBiplotPlot", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)

  fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)

  expect_type(fig$data, "list")
  expect_type(fig$layers, "list")
  expect_type(fig$scales, "environment")
  expect_type(fig$mapping, "list")
  expect_type(fig$theme, "list")
  expect_type(fig$coordinates, "environment")
  expect_type(fig$facet, "environment")
  expect_type(fig$plot_env, "environment")
  expect_type(fig$labels, "list")
})

test_that("Microarray GSE with missing values is handled correctly with nonInteractivePcaBiplotPlot", {
  allGset <- getGeoObject('GSE18380')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)

  fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)

  expect_type(fig$data, "list")
  expect_type(fig$layers, "list")
  expect_type(fig$scales, "environment")
  expect_type(fig$mapping, "list")
  expect_type(fig$theme, "list")
  expect_type(fig$coordinates, "environment")
  expect_type(fig$facet, "environment")
  expect_type(fig$plot_env, "environment")
  expect_type(fig$labels, "list")
})
