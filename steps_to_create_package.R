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
library(geo2rShinyApp)
pkgload::load_all(".")
geo2rShinyApp::myApp()
geo2rShinyApp::myApp()

# Mack the package a github repo
#setwd("./cats")
#install_github('cats','github_username')
