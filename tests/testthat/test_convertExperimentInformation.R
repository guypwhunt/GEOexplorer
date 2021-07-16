library(GEOexplorer)
context("convertExperimentInformation")

test_that("Microarray GSE is handled correctly", {
  allGset <- getGeoObject('GSE18380')

  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]

  gsetData <- extractPlatformGset(allGset, platform)

  experimentInformation <- extractExperimentInformation(gsetData)

  extractedExperimentInformation <- convertExperimentInformation(experimentInformation)

  expect_type(extractedExperimentInformation, 'character')
  expect_equal(extractedExperimentInformation[1], "<b> Paper Title: </b> <p> Autonomous plasmid-like replication of a conjugative transposon </p><b> Author's Name: </b> <p> Alan,D.,Grossman </p><b> Laboratory: </b> <p>  </p><b> Contact Information: </b> <p> clee2@mit.edu, adg@mit.edu </p><b> Paper URL: </b> <p> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18380 </p><b> PubMed ID: </b> <p> 19943900 </p><b> Abstract: </b> <p> Integrative and conjugative elements (ICEs), a.k.a., conjugative transposons, are mobile genetic elements involved in many biological processes, including the spread of antibiotic resistance. Unlike conjugative plasmids that are extra-chromosomal and replicate autonomously, ICEs are integrated in the chromosome and replicate passively during chromosomal replication. It is generally thought that ICEs do not replicate autonomously. We found that when induced, Bacillus subtilis ICEBs1 replicates as a plasmid. The ICEBs1 origin of transfer (oriT) served as the origin of replication and the conjugal DNA relaxase served as the replication initiation protein. Autonomous replication of ICEBs1 conferred genetic stability to the excised element, but was not required for mating. The B. subtilis helicase PcrA that mediates unwinding and replication of Gram-positive rolling circle replicating plasmids was required for ICEBs1 replication and mating. Nicking of oriT by the relaxase and unwinding by PcrA likely directs transfer of a single-strand of ICEBs1 into recipient cells.\n\nThis SuperSeries is composed of the SubSeries listed below. </p>")
})
