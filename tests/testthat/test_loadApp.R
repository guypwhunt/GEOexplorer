library(GEOexplorer)
context("loadApp")

test_that("loadApp function works", {
  app <- loadApp()

  expect_type(app, "list")
  expect_type(app$httpHandler, "closure")
  expect_type(app$serverFuncSource, "closure")
  expect_type(app$onStart, "NULL")
  expect_type(app$options, "list")
  expect_type(app$appOptions, "list")
})
