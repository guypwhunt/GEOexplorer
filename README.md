# GEOexplorer

GEOexplorer contains a Shiny app that enables exploratory data analysis and differential gene expression of gene expression analysis to be performed on microarray data held in the GEO database. The results are displayed as interactive graphs.

## Acknowledgements

The development of GEOexplorer was made possible because of the excellent code provided by GEO2R
(https://www.ncbi.nlm.nih.gov/geo/geo2r/).

## Installation

Use the devtools install_github function for installation e.g.

```R
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("guypwhunt/GEOexplorer")
```

Or install from BioConductor

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOexplorer")
```

## Usage

```R
library(GEOexplorer)
loadApp()
```

## Contributing
Pull requests are welcome. 

For any changes or issues, please open an issue to discuss what you would like to change.

## License
[GPL-3](https://choosealicense.com/licenses/gpl-3.0/)

## Authors
Guy Hunt with contributions from Fabrizio Smeraldi, Michael Barnes and Rafael Henkin
