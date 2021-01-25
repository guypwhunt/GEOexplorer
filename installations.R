#!/usr/bin/Rsript
# ---------------------------------------------------------
# File name      : installations.R
# Authors       : Guy Hunt
# Description   : Pre-Installation libraries
# ---------------------------------------------------------

#############################################################################
#       Load necessary dependencies, if not previously installed            #
#############################################################################

install.packages("pacman")
source("http://bioconductor.org/biocLite.R")

pacman::p_load("GEOquery", "gage", "gageData", "GO.db", "pathview", "limma", "impute")
# Species database
pacman::p_load("org.Ag.eg.db", "org.At.tair.db", "org.Bt.eg.db", "org.Ce.eg.db", "org.Cf.eg.db", "org.Dm.eg.db", "org.Dr.eg.db", "org.EcK12.eg.db", "org.EcSakai.eg.db", "org.Gg.eg.db", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Pf.plasmo.db", "org.Pt.eg.db", "org.Rn.eg.db", "org.Sc.sgd.db", "org.Ss.eg.db", "org.Xl.eg.db")
pacman::p_load("argparser", "Cairo", "dendextend", "DMwR", "ggplot2", "gplots", "jsonlite", "pheatmap", "plyr", "RColorBrewer", "reshape2", "squash")

#############################################################################
#            For additional information on required packages                #
#############################################################################

# GEOquery      - http://bioconductor.org/packages/release/bioc/html/GEOquery.html
# gage          - http://bioconductor.org/packages/release/bioc/html/gage.html
# gageData      - https://bioconductor.org/packages/release/data/experiment/html/gageData.html
# GO.db         - https://bioconductor.org/packages/release/data/annotation/html/GO.db.html
# pathview      - http://bioconductor.org/packages/release/bioc/html/pathview.html
# limma         - https://bioconductor.org/packages/release/bioc/html/limma.html
# org.Mm.eg.db  - http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
# argparser     - https://cran.r-project.org/web/packages/argparser/index.html
# Cairo         - https://cran.r-project.org/web/packages/Cairo/index.html
# dendextend    - https://cran.r-project.org/web/packages/dendextend/index.html
# impute        - http://svitsrv25.epfl.ch/R-doc/library/impute/html/impute.knn.html
# DMwR          - https://cran.r-project.org/web/packages/DMwR/index.html
# ggplot2       - https://cran.r-project.org/web/packages/ggplot2/index.html
# gplots        - https://cran.r-project.org/web/packages/gplots/index.html
# jsonlite      - https://cran.r-project.org/web/packages/jsonlite/index.html
# pheatmap      - https://cran.r-project.org/web/packages/pheatmap/index.html
# plyr          - https://cran.r-project.org/web/packages/plyr/index.html
# RColorBrewer  - https://cran.r-project.org/web/packages/RColorBrewer/index.html
# reshape2      - https://cran.r-project.org/web/packages/reshape2/index.html
# squash        - https://cran.r-project.org/web/packages/squash/index.html
