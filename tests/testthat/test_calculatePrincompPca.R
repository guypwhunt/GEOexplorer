library(GEOexplorer)
context("calculatePrincompPca")

test_that("Microarray GSE without missing values is handled correctly with Princomp PCA", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)

  expect_type(pcaPrincompDataInput, 'list')
  expect_s3_class(pcaPrincompDataInput, 'princomp')
})

test_that("Microarray GSE with missing values is handled correctly with Princomp PCA", {
  allGset <- getGeoObject('GSE18380')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  naOmitInput <- calculateNaOmit(expressionData)

  pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)

  expect_type(pcaPrincompDataInput, 'list')
  expect_s3_class(pcaPrincompDataInput, 'princomp')
})
