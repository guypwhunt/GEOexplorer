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
  expect_equal(experimentInformation@abstract, "Microarray Analysis of Space-flown Murine Thymus Tissue Reveals Changes in Gene Expression Regulating Stress and Glucocorticoid Receptors. We used microarrays to detail the gene expression of space-flown thymic tissue and identified distinct classes of up-regulated genes during this process. We report here microarray gene expression analysis in young adult C57BL/6NTac mice at 8 weeks of age after exposure to spaceflight aboard the space shuttle (STS-118) for a period of 13 days. Upon conclusion of the mission, thymus lobes were extracted from space flown mice (FLT) as well as age- and sex-matched ground control mice similarly housed in animal enclosure modules (AEM). mRNA was extracted and an automated array analysis for gene expression was performed. Examination of the microarray data revealed 970 individual probes that had a 1.5 fold or greater change. When these data were averaged (n=4), we identified 12 genes that were significantly up- or down-regulated by at least 1.5 fold after spaceflight (pâ‰¤0.05). Together, these data demonstrate that spaceflight induces significant changes in the thymic mRNA expression of genes that regulate stress, glucocorticoid receptor metabolism, and T cell signaling activity. These data explain, in part, the reported systemic compromise of the immune system after exposure to the microgravity of space.")
  expect_equal(experimentInformation@url, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18388")
  expect_equal(experimentInformation@pubMedIds, "20213684")
})
