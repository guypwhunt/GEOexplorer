if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

gc()


devtools::install_github("guypwhunt/GEOexplorer")

gc()

GEOexplorer::loadApp()
