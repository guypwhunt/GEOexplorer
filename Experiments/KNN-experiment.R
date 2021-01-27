source("geoIntegration.R")
source("geo2rDataVisualisation.R")
source("analyticsFunctions.R")
library(impute)

data <-getGeoData("GSE18380", "GPL4694")
data <-extractGeoData(data, "Auto-Detect")

# str(data)
# head(data)
# summary(data)

imputation <- impute.knn(data)
ex <- imputation$data
#summary(ex)

ex2 <- knnDataTransformation(data, 'Yes')

print('Head Data')
head(data)


print('Head EX')
head(ex)

print('Head EX2')
head(ex2)
