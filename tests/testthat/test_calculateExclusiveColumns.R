library(GEOexplorer)
context("calculateExclusiveColumns")

test_that("Microarray GSE without missing values is handled correctly with calculateExclusiveColumns", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  columnNames <- extractSampleNames(expressionData)

  numberOfColumns <- length(columnNames)
  numberOfColumns <- numberOfColumns + 1
  halfNumberOfColumns <- ceiling(numberOfColumns/2)
  i <- 0

  group1 <- c()
  group2 <- c()

  for (name in columnNames) {
    if (i < halfNumberOfColumns) {
      group1 <- c(group1, name)
      i <- i +1
    } else {
      group2 <- c(group2, name)
      i <- i +1
    }
  }

  column2 <- calculateExclusiveColumns(columnNames, group1)

  expect_type(column2, "character")
  expect_equal(column2[1], "GSM458599")
  expect_equal(column2[2], "GSM458600")
  expect_equal(column2[3], "GSM458601")
  expect_equal(column2[4], "NA")
  expect_equal(length(column2), 3)

})

test_that("Microarray GSE with missing values is handled correctly with calculateExclusiveColumns", {
  allGset <- getGeoObject('GSE18380')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  expressionData <- extractExpressionData(gsetData)

  columnNames <- extractSampleNames(expressionData)

  numberOfColumns <- length(columnNames)
  numberOfColumns <- numberOfColumns + 1
  halfNumberOfColumns <- ceiling(numberOfColumns/2)
  i <- 0

  group1 <- c()
  group2 <- c()

  for (name in columnNames) {
    if (i < halfNumberOfColumns) {
      group1 <- c(group1, name)
      i <- i +1
    } else {
      group2 <- c(group2, name)
      i <- i +1
    }
  }

  column2 <- calculateExclusiveColumns(columnNames, group1)

  expect_type(column2, "character")
  expect_equal(column2[1], "GSM455783")
  expect_equal(column2[2], "GSM455784")
  expect_equal(column2[3], "GSM455785")
  expect_equal(column2[4], "GSM455786")
  expect_equal(column2[5], "GSM455787")
  expect_equal(column2[6], "NA")
  expect_equal(length(column2), 5)
})

