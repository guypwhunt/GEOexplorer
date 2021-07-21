if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

gc()


devtools::install_github("guypwhunt/GEOexplorer", ref = "master")

gc()

GEOexplorer::loadApp()
