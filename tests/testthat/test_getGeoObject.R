library(GEOexplorer)
context("getGeoObject")

test_that("Empty GSE is handled correctly", {
  allGset <- getGeoObject('GSE11413')

  expect_type(allGset, 'list')
  expect_type(allGset[[1]], 'S4')
  expect_s4_class(allGset[[1]], 'ExpressionSet')
  expect_equal(nrow(pData(allGset[[1]])),12)
  expect_equal(nrow(fData(allGset[[1]])),0)
})


test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  expect_type(allGset, 'list')
  expect_type(allGset[[1]], 'S4')
  expect_s4_class(allGset[[1]], 'ExpressionSet')
  expect_equal(nrow(pData(allGset[[1]])),8)
  expect_equal(nrow(fData(allGset[[1]])),35557)
})

test_that("Microarray GSE has populated experimentData", {
  allGset <- getGeoObject("GSE18388")

  ed <- experimentData(allGset[[1]])
  expect_equal(pubMedIds(ed), "20213684")

  ei <- expinfo(ed)
  expect_equal(ei[1], "Ty,W,Lebsack", ignore_attr = TRUE)
  expect_equal(ei[2], "", ignore_attr = TRUE) #lab
  expect_equal(ei[3], "lebsack@email.arizona.edu", ignore_attr = TRUE)
  expect_equal(ei[4], "Microarray Analysis of Space-flown Murine Thymus Tissue", ignore_attr = TRUE)
  expect_equal(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18388", ignore_attr = TRUE) #url
})
