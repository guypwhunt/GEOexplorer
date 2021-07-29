library(GEOexplorer)
context("calculateLogTransformation")

test_that("Microarray GSE without missing values is handled correctly with NA omit", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  expect_type(naOmitInput, 'double')
  expect_equal(ncol(naOmitInput),8)
  expect_equal(nrow(naOmitInput),35557)
  expect_equal(naOmitInput[1,1],11.711505)
})
