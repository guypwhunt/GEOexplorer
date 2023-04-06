library(GEOexplorer)

test_that("GEO search works",
          {
            searchTerm <- "TXT"
            firstResultNumber <- "0"
            lastResultNumber <- "1"
            
            x <- searchGeo(searchTerm, firstResultNumber, lastResultNumber)
          })