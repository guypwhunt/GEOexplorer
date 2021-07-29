library(GEOexplorer)
context("interactiveBoxAndWhiskerPlot")

test_that("Microarray GSE without missing values is handled correctly with interactiveBoxAndWhiskerPlot", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  fig <- interactiveBoxAndWhiskerPlot(expressionData, 'GSE18388', platform)
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
