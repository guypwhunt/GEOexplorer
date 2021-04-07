# Download necessary libraries
#install.packages("devtools")
#devtools::install_github("klutometis/roxygen")

# Import necessary libraries
library(devtools)
library(roxygen2)

# Create cats package
setwd("..")
create("geo2rShinyAppPackage")

# Add files to R folder

# Create documentation
setwd("./geo2rShinyAppPackage")
document()

# Install the package
setwd("..")
install("geo2rShinyAppPackage")

# Test cats worked
library(geo2rShinyAppPackage)
?geo2rShinyAppPackage

# Macke the package a github repo
setwd("./cats")
install_github('cats','github_username')
