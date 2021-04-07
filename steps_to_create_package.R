# Download neccessary libraries
install.packages("devtools")
library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)

# Create cats package
setwd("..")
create("cats")

# Create documentation
setwd("./cats")
document()

# Install the package
setwd("..")
install("cats")

# Test cats worked
?cat_function

# Macke the package a github repo
setwd("./cats")
install_github('cats','github_username')
