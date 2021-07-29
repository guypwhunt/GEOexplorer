library(GEOexplorer)
context("nonInteractiveBoxAndWhiskerPlot")

test_that("Microarray GSE without missing values is handled correctly with nonInteractiveBoxAndWhiskerPlot", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  fig <- nonInteractiveBoxAndWhiskerPlot(ex = expressionData, geoAccessionCode = 'GSE18388', platform = platform)

  expect_type(fig, 'list')
  expect_type(fig$stats, 'double')
  expect_type(fig$n, 'double')
  expect_type(fig$conf, 'double')
  expect_type(fig$out, 'double')
  expect_type(fig$group, 'double')
  expect_type(fig$names, 'character')
})
