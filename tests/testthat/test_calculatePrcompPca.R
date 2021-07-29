library(GEOexplorer)
context("calculatePrcompPca")

test_that("Microarray GSE without missing values is handled correctly with Prcomp PCA", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)

  expect_type(pcaPrcompDataInput, 'list')
  expect_s3_class(pcaPrcompDataInput, 'prcomp')
})
