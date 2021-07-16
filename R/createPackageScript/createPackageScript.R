# Create a description file
usethis::use_description()

# Load all files
devtools::load_all()

# Import necessary libraries
library(devtools)
library(roxygen2)

# Set working directory
setwd("C:/Users/guypw/OneDrive/Documents/GEOexplorer")

# Remove documentation
rm(list = c("calculateAutoLogTransformApplication", "calculateDifferentialGeneExpression", "calculateDifferentialGeneExpressionSummary", "calculateEachGroupsSamples", "calculateEachGroupsSamplesFromDataFrame", "calculateExclusiveColumns", "calculateKnnImpute", "calculateLogTransformation", "calculateNaOmit", "calculatePca", "calculatePrincompPca", "calculateTopDifferentiallyExpressedGenes", "convertAdjustment", "convertExperimentInformation", "extractExperimentInformation", "extractExpressionData", "extractGeoData", "extractPlatformGset", "extractPlatforms", "extractSampleDetails", "extractSampleInformation", "extractSampleNames", "getGeoObject", "interactiveBoxAndWhiskerPlot", "interactiveDensityPlot", "interactiveHeatMapPlot", "interactiveHistogramPlot", "interactiveMeanDifferencePlot", "interactiveMeanVariancePlot", "interactivePcaScreePlot", "interactivePrincompPcaIndividualsPlot", "interactivePrincompPcaScreePlot", "interactivePrincompPcaVariablesPlot", "interactiveQQPlot", "interactiveThreeDDesnityPlot", "interactiveUmapPlot", "interactiveVolcanoPlot", "loadApp", "nonInteractiveBoxAndWhiskerPlot", "nonInteractiveCorrelationMatrixPlot", "nonInteractiveDesnityPlot", "nonInteractiveHistogramPlot", "noninteractiveMeanDifferencePlot", "nonInteractiveMeanVariancePlot", "nonInteractivePcaBiplotPlot", "nonInteractivePcaIndividualsPlot", "nonInteractivePcaScreePlot", "nonInteractivePcaVariablesPlot", "nonInteractiveQQPlot", "nonInteractiveUmapPlot", "nonInteractiveVennDiagramPlot", "nonInteractiveVolcanoPlot"))
rm(list = c("datasetInformationUi", "differentialGeneExpressionAnalysisUi", "exploratoryDataAnalysisUi", "sideBarUi", "sourceDatasetInformationUi", "sourceDifferentialGeneExpressionAnalysisUi", "sourceExploratoryDataAnalysisUi", "sourceSideBarUi"))
rm(list = c("calculatePrcompPca", "interactivePrcompPcaIndividualsPlot", "interactivePrcompPcaScreePlot", "interactivePrcompPcaVariablesPlot", "sourceServer"))
rm(list = c("interactiveThreeDDensityPlot", "nonInteractiveDensityPlot"))

# Create documentation
document()
