# GEOexplorer

GEOexplorer enables exploratory data analysis and differential gene expression of gene expression analysis to be performed on data held in the GEO database. The outputs are interactive graphs. GEOexplorer also contains a shiny app that can be used to perform exploratory data analysis and differential gene expression of gene expression analysis.

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
GEOexplorer::loadApp()
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[GPL-3](https://choosealicense.com/licenses/gpl-3.0/)

## Authors
Guy Hunt with contributions from Fabrizio Smeraldi, Michael Barnes and Rafael Henkin
