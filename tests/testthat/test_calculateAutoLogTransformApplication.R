library(GEOexplorer)
context("calculateAutoLogTransformApplication")

test_that("Microarray GSE is handled correctly when it should not apply log transofrmation", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  autoLogInformation <- calculateAutoLogTransformApplication(expressionData)

  expect_type(autoLogInformation, 'character')
  expect_equal(autoLogInformation,"The auto-detect option did not apply log transformation.")
})
