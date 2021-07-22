# Source the files
try(source("R/app.R"))
try(source("R/backendDifferentialGeneExpressionDataFunctions.R"))
try(source("R/backendExploratoryDataAnalysisDataFunctions.R"))
try(source("R/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R"))
try(source("R/interactiveExploratoryDataAnalysisDataVisualizationFunctions.R"))
try(source("R/nonInteractiveDifferentialGeneExpressionDataVisualizationFunctions.R"))
try(source("R/nonInteractiveExploratoryDataAnalysisDataVisualizationFunctions.R"))
try(source("R/appUiComponents.R"))
try(source("R/appServerComponents.R"))


# Launch app
GEOexplorer::loadApp()
