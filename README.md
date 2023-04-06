# GEOexplorer

GEOexplorer is an R/Bioconductor package and web application
that enables users to perform gene expression analysis. 

## Acknowledgements

The development of GEOexplorer was made possible because of the
excellent code provided by GEO2R
(https://www.ncbi.nlm.nih.gov/geo/geo2r/).

## Web application
GEOexplorer is also available as a website on the following link:
https://geoexplorer.rosalind.kcl.ac.uk/

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

For any changes or issues, please open an issue to discuss what you
would like to change.

## License
[GPL-3](https://choosealicense.com/licenses/gpl-3.0/)

## Authors
Guy Hunt with contributions from Fabrizio Smeraldi, Michael Barnes 
and Rafael Henkin.

## Future Updates
Please see our Google [doc](https://docs.google.com/spreadsheets/d/14UW_-9pLbqJsh5pYUcSy2hE9CZCJ1nw5o5izuqo-fOw/edit?usp=sharing) for the planned feature updates.