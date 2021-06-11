# Download necessary libraries
#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

# Wrap shiny app in a function

# Create a description file
usethis::use_description()

# Load all files
devtools::load_all()

# Remove all source()

# Move all the code to the R/ folder

# Import necessary libraries
library(devtools)
library(roxygen2)

# Create documentation
#setwd("./documents")
setwd("./geo2rShinyApp")
document()

# Build the package

# Install the package
setwd("..")
install("geo2rShinyApp")

# Test cats worked
library(GEOexplorer)
pkgload::load_all(".")
GEOexplorer::loadApp()

# Mack the package a github repo
#setwd("./cats")
#install_github('cats','github_username')

# Push to github

# Install from GitHub
install.packages("devtools")
library(devtools)
install_github("guypwhunt/r_shiny_geo2r_visulisation", force = TRUE, ref = "main")


# Test library
install_github("guypwhunt/GEOexplorer", force = TRUE, ref = "master")
library(GEOexplorer)
geo2rShinyApp::loadApp()
help(package = "GEOexplorer")

# Test functions
geoAccessionCode <- "GSE18388"
logTransformation <- "Auto-Detect"
knnTransformation <- "Yes"

allGset <- getGeoObject(geoAccessionCode)
platforms <- extractPlatforms(allGset)
gsetData <- extractPlatformGset(allGset, platforms[1])
expressionData <- extractExpressionData(gsetData)
logExpressionData <- calculateLogTransformation(expressionData, logTransformation)
knnDataInput <- calculateKnnImpute(logExpressionData, knnTransformation)
pcaPrincompExpressionData <- calculatePrincompPca(knnDataInput)
fig <- interactivePrincompPcaScreePlot(pcaPrincompExpressionData, geoAccessionCode)
fig



