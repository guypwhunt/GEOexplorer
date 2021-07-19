library(GEOexplorer)
context("convertAdjustment")

test_that("convertAdjustment can handle all adjustment values", {
  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "fdr")

  pValueAdjustment <- "Benjamini & Yekutieli"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "BY")

  pValueAdjustment <- "Bonferroni"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "bonferroni")

  pValueAdjustment <- "Holm"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "holm")

  pValueAdjustment <- "None"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "none")
})
