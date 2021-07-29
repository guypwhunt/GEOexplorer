library(GEOexplorer)
context("calculateDifferentialGeneExpression")

test_that("Microarray GSE without missing values is handled correctly with calculateDifferentialGeneExpression", {
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

  gsms <- calculateEachGroupsSamples(columnNames,group1, group2)

  limmaPrecisionWeights <- "Yes"
  forceNormalization <- "Yes"
  fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, expressionData)

  expect_type(fit2, "list")
  expect_type(fit2$coefficients, "double")
  expect_type(fit2$sigma, "double")
  expect_type(fit2$cov.coefficients, "double")
  expect_type(fit2$rank, "integer")
  expect_type(fit2$Amean, "double")
  expect_type(fit2$design, "double")
  expect_type(fit2$df.prior, "double")
  expect_type(fit2$var.prior, "double")
  expect_type(fit2$s2.post, "double")
  expect_type(fit2$df.total, "double")
  expect_type(fit2$lods, "double")
  expect_type(fit2$F.p.value, "double")
  expect_type(fit2$stdev.unscaled, "double")
  expect_type(fit2$df.residual, "double")
  expect_type(fit2$pivot, "integer")
  expect_type(fit2$genes, "list")
  expect_type(fit2$method, "character")
  expect_type(fit2$contrasts, "double")
  expect_type(fit2$s2.prior, "double")
  expect_type(fit2$proportion, "double")
  expect_type(fit2$t, "double")
  expect_type(fit2$p.value, "double")
  expect_type(fit2$F, "double")
})
