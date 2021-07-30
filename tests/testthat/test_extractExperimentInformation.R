library(GEOexplorer)
context("extractExperimentInformation")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18388')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  experimentInformation <- extractExperimentInformation(gsetData)

  expect_type(experimentInformation, 'S4')
  expect_s4_class(experimentInformation, 'MIAME')
  expect_equal(experimentInformation@name,"Ty,W,Lebsack")
  expect_equal(experimentInformation@lab, "")
  expect_equal(experimentInformation@contact, "lebsack@email.arizona.edu")
  expect_equal(experimentInformation@title, "Microarray Analysis of Space-flown Murine Thymus Tissue")
  expect_equal(nchar(experimentInformation@title), 55)
  expect_equal(experimentInformation@url, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18388")
  expect_equal(experimentInformation@pubMedIds, "20213684")
})
