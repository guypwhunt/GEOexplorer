# Source the files
try(source("R/app.R"))
try(source("R/backendDifferentialGeneExpressionDataFunctions.R"))
try(source("R/backendExploratoryDataAnalysisDataFunctions.R"))
try(source("R/interactiveDifferentialGeneExpressionDataVisualizationFunctions.R"))
try(source("R/interactiveExploratoryDataAnalysisDataVisualizationFunctions.R"))
try(source("R/nonInteractiveDifferentialGeneExpressionDataVisualizationFunctions.R"))
try(source("R/nonInteractiveExploratoryDataAnalysisDataVisualizationFunctions.R"))
try(source("R/appSideBarUi.R"))
try(source("R/appDatasetInformationUi.R"))
try(source("R/appExploratoryDataAnalysisUi.R"))
try(source("R/appDifferentialGeneExpressionAnalysisUi.R"))


# Launch app
loadApp()
