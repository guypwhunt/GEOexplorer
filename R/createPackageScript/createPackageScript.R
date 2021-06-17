# Create a description file
usethis::use_description()

# Load all files
devtools::load_all()

# Import necessary libraries
library(devtools)
library(roxygen2)

# Set working directory
setwd("./geo2rShinyApp")

# Create documentation
document()

rm(list = c("nonInteractivePcaBiplotPlot", "nonInteractivePcaIndividualsPlot", "nonInteractivePcaScreePlot", "nonInteractivePcaVariablesPlot", "nonInteractiveQQPlot", "nonInteractiveUmapPlot", "nonInteractiveVennDiagramPlot", "nonInteractiveVolcanoPlot"))
