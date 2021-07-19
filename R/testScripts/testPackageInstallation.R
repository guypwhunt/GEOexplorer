if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

gc()


devtools::install_github("guypwhunt/GEOexplorer", ref = "master")

gc()

#git clone https://github.com/guypwhunt/GEOexplorer/tree/master



GEOexplorer::loadApp()
