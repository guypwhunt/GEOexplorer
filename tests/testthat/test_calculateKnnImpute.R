library(GEOexplorer)
context("calculateKnnImpute")

test_that("Microarray GSE without missing values is handled correctly with KNN imputation", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  logTransformation <- "Auto-Detect"  #
  dataInput <- calculateLogTransformation(expressionData, logTransformation)

  knnTransformation <- "Yes"
  knnDataInput <- calculateKnnImpute(dataInput, knnTransformation)

  expect_type(knnDataInput, 'double')
  expect_equal(ncol(knnDataInput),8)
  expect_equal(nrow(knnDataInput),35557)
  expect_equal(knnDataInput[1,1],11.711505)
})
