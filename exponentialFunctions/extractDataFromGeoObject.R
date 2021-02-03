# Test to extract data from 
library(GEOquery)
geo <- "GSE18380"
geoFile <- getGEO(GEO = geo)
geoFile <- geoFile[[1]]
geoFile <- attr(geoFile, 'experimentData')
otherInformation <- attr(geoFile, 'other')
platformId <- attr(otherInformation, 'platform_id')
author <- attr(geoFile, 'name')
attributes(geoFile)
print(author)