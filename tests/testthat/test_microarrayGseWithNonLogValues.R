library(GEOexplorer)
context("Missing Values")

test_that("Microarray GSE with missing values is handled correctly by all functions", {
  # Input Values
  logTransformation <- "Auto-Detect"
  knnTransformation <- "Yes"
  knn <- 2
  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
  limmaPrecisionWeights <- "Yes"
  forceNormalization <- "Yes"
  platformAnnotation <- "NCBI generated"
  significanceLevelCutOff <- 0.05

  # Get the GEO data for all platforms
  geoAccessionCode <- "GSE2"
  allGset <- getGeoObject(geoAccessionCode)
  ed <- experimentData(allGset[[1]])
  expect_equal(pubMedIds(ed), "")
  ei <- expinfo(ed)
  expect_equal(ei[1], "Yoshihiro,,Kagami", ignore_attr = TRUE)
  expect_equal(ei[2], "", ignore_attr = TRUE) #lab
  expect_equal(ei[3], "ykagami@brain.riken.go.jp", ignore_attr = TRUE)
  expect_equal(ei[4], "Cerebellar development", ignore_attr = TRUE)
  expect_equal(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2", ignore_attr = TRUE)

  # Extract platforms
  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]
  expect_type(platforms, 'character')
  expect_type(platform, 'character')
  expect_equal(platform,"GPL8")

  # Extract the GEO2R data from the specified platform
  gsetData <- extractPlatformGset(allGset, platform)
  expect_type(gsetData, 'S4')
  expect_s4_class(gsetData, 'ExpressionSet')
  expect_equal(nrow(pData(gsetData)),5)
  expect_equal(nrow(fData(gsetData)),897)

  # Extract the experiment information
  experimentInformation <- extractExperimentInformation(gsetData)
  expect_type(experimentInformation, 'S4')
  expect_s4_class(experimentInformation, 'MIAME')
  expect_equal(experimentInformation@name,"Yoshihiro,,Kagami")
  expect_equal(experimentInformation@lab, "")
  expect_equal(experimentInformation@contact, "ykagami@brain.riken.go.jp")
  expect_equal(experimentInformation@title, "Cerebellar development")
  expect_equal(nchar(experimentInformation@title), 22)
  expect_equal(experimentInformation@url, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2")
  expect_equal(experimentInformation@pubMedIds, "")

  # Extract Sample Information
  sampleInfo <- extractSampleInformation(gsetData)
  expect_type(sampleInfo, 'list')
  expect_equal(nrow(sampleInfo),5)
  expect_equal(ncol(sampleInfo),30)

  # Extract expression data
  expressionData <- extractExpressionData(gsetData)
  expect_type(expressionData, 'double')
  expect_equal(ncol(expressionData),5)
  expect_equal(nrow(expressionData),897)

  # Get column Details
  columnInfo <- extractSampleDetails(gsetData)
  expect_type(columnInfo, 'list')
  expect_equal(ncol(columnInfo),3)
  expect_equal(nrow(columnInfo),5)

  # Is log transformation auto applied
  autoLogInformation <- calculateAutoLogTransformApplication(expressionData)
  expect_type(autoLogInformation, 'character')
  expect_equal(autoLogInformation,"The auto-detect option applied log transformation.")

  # Get a list of all the columns
  columns <- extractSampleNames(expressionData)
  expect_type(columns, 'character')
  expect_equal(columns[1],"GSM50")

  # Apply log transformation to expression data if necessary
  dataInput <- calculateLogTransformation(expressionData, logTransformation)
  expect_type(dataInput, 'double')
  expect_equal(ncol(dataInput),5)
  expect_equal(nrow(dataInput),897)
  expect_equal(dataInput[1,1],6.5833085)

  # Perform KNN transformation on log expression data if necessary
  knnDataInput <- calculateKnnImpute(dataInput, "Yes")
  expect_type(knnDataInput, 'double')
  expect_equal(ncol(knnDataInput),5)
  expect_equal(nrow(knnDataInput),897)
  expect_equal(knnDataInput[1,1],6.5833085)

  # Get a list of all the columns in the KNN output
  knnColumns <- extractSampleNames(knnDataInput)

  # Get knn output column Details
  knnColumnInfo <- extractSampleDetails(gsetData)
  knnColumnInfo <- knnColumnInfo[knnColumns,]

  # Remove all incomplete rows
  naOmitInput <- calculateNaOmit(knnDataInput)
  expect_type(naOmitInput, 'double')
  expect_equal(ncol(naOmitInput),5)
  expect_equal(nrow(naOmitInput),897)
  expect_equal(naOmitInput[1,1],6.5833085)

  # Perform Princomp PCA analysis on KNN transformation expression data
  pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
  expect_type(pcaPrincompDataInput, 'list')
  expect_s3_class(pcaPrincompDataInput, 'princomp')

  # Perform Prcomp PCA analysis on KNN transformation expression data
  pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
  expect_type(pcaPrcompDataInput, 'list')
  expect_s3_class(pcaPrcompDataInput, 'prcomp')

  # Extract Experiment Information
  extractedExperimentInformation <- convertExperimentInformation(experimentInformation)
  expect_type(extractedExperimentInformation, 'character')
  expect_equal(extractedExperimentInformation[1], "<b> Paper Title: </b> <p> Cerebellar development </p><b> Author's Name: </b> <p> Yoshihiro,,Kagami </p><b> Laboratory: </b> <p>  </p><b> Contact Information: </b> <p> ykagami@brain.riken.go.jp </p><b> Paper URL: </b> <p> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2 </p><b> PubMed ID: </b> <p>  </p><b> Abstract: </b> <p> The Affymetrix GeneChip Mu11K was used to analyze the gene expression profile in developing mouse cerebellum (two GeneChips per E18, P7, P14, P21, and P56) to assist in the understanding of the genetic basis of cerebellar development in mice.\nThe analysis showed 81.6% (10,321/12,654) of the genes represented on the GeneChip were expressed in the postnatal cerebellum, and among those, 8.7% (897/10,321) were differentially expressed with more than a two-fold change in their maximum and minimum expression levels during the developmental time course.\nThe expression data (mean signal in relative unit) of all of these 897 differentially expressed genes were listed in GSM50(for E18), GSM51(for P7), GSM52(for P14), GSM53(for P21), and GSM54(for P56)  as well as our homepage at http://www.brain.riken.go.jp/labs/lm\nKeywords = mouse cerebellum development\nKeywords: time-course </p>")

  # Non-Interactive Box-and-Whisker Plot
  fig <- nonInteractiveBoxAndWhiskerPlot(ex = knnDataInput, geoAccessionCode = geoAccessionCode, platform = platform)
  expect_type(fig, 'list')
  expect_type(fig$stats, 'double')
  expect_type(fig$n, 'double')
  expect_type(fig$conf, 'double')
  expect_type(fig$out, 'double')
  expect_type(fig$group, 'double')
  expect_type(fig$names, 'character')

  # Interactive Box-and-Whisker Plot
  fig <- interactiveBoxAndWhiskerPlot(knnDataInput, geoAccessionCode, platform)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Non-Interactive Density Plot
  fig <- nonInteractiveDensityPlot(ex = naOmitInput, geoAccessionCode = geoAccessionCode, platform = platform)
  expect_type(fig, 'list')
  expect_type(fig$X, 'double')
  expect_type(fig$Y, 'double')

  # Interactive Density Plot
  fig <- interactiveDensityPlot(naOmitInput, geoAccessionCode, platform)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # 3D Interactive Density Plot
  fig <- interactiveThreeDDensityPlot(naOmitInput, geoAccessionCode, platform)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive UMAP
  fig <- interactiveUmapPlot(naOmitInput, knn, geoAccessionCode)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Mean Variance Plot
  fig <- interactiveMeanVariancePlot(naOmitInput, geoAccessionCode, gsetData)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Princomp PCA Scree Plot
  fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput, geoAccessionCode)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Princomp PCA Individual Plot
  fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput, geoAccessionCode, gsetData)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Princomp PCA Variables Plot
  fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput, geoAccessionCode)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Prcomp PCA Scree Plot
  fig <- interactivePrcompPcaScreePlot(pcaPrcompDataInput, geoAccessionCode)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Prcomp PCA Individual Plot
  fig <- interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput, geoAccessionCode, gsetData)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Interactive Prcomp PCA Variables Plot
  fig <- interactivePrcompPcaVariablesPlot(pcaPrcompDataInput, geoAccessionCode)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Correlation Matrix of samples
  fig <- interactiveHeatMapPlot(naOmitInput)
  fig
  expect_type(fig, 'list')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')

  # Non-Interactive UMAP
  fig <- nonInteractiveUmapPlot(naOmitInput, knn, geoAccessionCode)
  expect_type(fig, 'list')
  expect_type(fig$x, 'double')
  expect_type(fig$y, 'double')

  # Non-Interactive Mean Variance Plot
  fig <- nonInteractiveMeanVariancePlot(naOmitInput, geoAccessionCode)

  # Non-Interactive Princomp PCA Scree Plot
  fig <- nonInteractivePcaScreePlot(pcaPrincompDataInput)
  fig
  expect_type(fig, 'list')
  expect_type(fig$data, 'list')
  expect_type(fig$layers, 'list')
  expect_type(fig$scales, 'environment')
  expect_type(fig$mapping, 'list')
  expect_type(fig$theme, 'list')
  expect_type(fig$coordinates, 'environment')
  expect_type(fig$plot_env, 'environment')
  expect_type(fig$labels, 'list')

  # Non-Interactive Princomp PCA Individual Plot
  fig <- nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
  fig
  expect_type(fig$data, "list")
  expect_type(fig$layers, "list")
  expect_type(fig$scales, "environment")
  expect_type(fig$mapping, "list")
  expect_type(fig$theme, "list")
  expect_type(fig$coordinates, "environment")
  expect_type(fig$facet, "environment")
  expect_type(fig$plot_env, "environment")
  expect_type(fig$labels, "list")

  # Non-Interactive Princomp PCA Variables Plot
  fig <- nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
  fig
  expect_type(fig$data, "list")
  expect_type(fig$layers, "list")
  expect_type(fig$scales, "environment")
  expect_type(fig$mapping, "list")
  expect_type(fig$theme, "list")
  expect_type(fig$coordinates, "environment")
  expect_type(fig$facet, "environment")
  expect_type(fig$plot_env, "environment")
  expect_type(fig$labels, "list")

  # Non-Interactive Princomp PCA Individual and Variables Bilot
  fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
  fig

  # Differential gene expression analysis functions
  # Get column names
  columnNames <- extractSampleNames(expressionData)

  # Define Groups
  numberOfColumns <- length(columnNames)
  numberOfColumns <- numberOfColumns + 1
  halfNumberOfColumns <- ceiling(numberOfColumns/2)
  i <- 0

  group1 <- c()
  group2 <- c()

  for (name in columnNames) {
    if (i < halfNumberOfColumns) {
      group1 <- c(group1, name)
      i <- i +1
    } else {
      group2 <- c(group2, name)
      i <- i +1
    }
  }

  # Select columns in group2
  column2 <- calculateExclusiveColumns(columnNames, group1)
  expect_type(column2, "character")
  expect_equal(column2[1], "GSM53")
  expect_equal(column2[2], "GSM54")
  expect_equal(column2[3], "NA")
  expect_equal(length(column2), 2)

  # Calculate gsms
  gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
  expect_type(gsms, "character")
  expect_equal(gsms, "00011")
  expect_equal(nchar(gsms), 5)

  # Convert P value adjustment
  pValueAdjustment <- "Benjamini & Hochberg (False discovery rate)"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "fdr")

  pValueAdjustment <- "Benjamini & Yekutieli"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "BY")

  pValueAdjustment <- "Bonferroni"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "bonferroni")

  pValueAdjustment <- "Holm"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "holm")

  pValueAdjustment <- "None"
  adjustment <- convertAdjustment(pValueAdjustment)
  expect_type(adjustment, "character")
  expect_equal(adjustment, "none")

  adjustment <- convertAdjustment(pValueAdjustment)

  # Get fit 2
  fit2 <- calculateDifferentialGeneExpression(gsms, limmaPrecisionWeights, forceNormalization, gsetData, expressionData)
  expect_type(fit2, "list")
  expect_type(fit2$coefficients, "double")
  expect_type(fit2$sigma, "double")
  expect_type(fit2$cov.coefficients, "double")
  expect_type(fit2$rank, "integer")
  expect_type(fit2$Amean, "double")
  expect_type(fit2$design, "double")
  expect_type(fit2$df.prior, "double")
  expect_type(fit2$var.prior, "double")
  expect_type(fit2$s2.post, "double")
  expect_type(fit2$df.total, "double")
  expect_type(fit2$lods, "double")
  expect_type(fit2$F.p.value, "double")
  expect_type(fit2$stdev.unscaled, "double")
  expect_type(fit2$df.residual, "double")
  expect_type(fit2$pivot, "integer")
  expect_type(fit2$genes, "list")
  expect_type(fit2$method, "character")
  expect_type(fit2$contrasts, "double")
  expect_type(fit2$s2.prior, "double")
  expect_type(fit2$proportion, "double")
  expect_type(fit2$t, "double")
  expect_type(fit2$p.value, "double")
  expect_type(fit2$F, "double")

  # Print Top deferentially expressed genes
  tT <- calculateTopDifferentiallyExpressedGenes(fit2, adjustment)
  expect_type(tT, "list")
  expect_type(tT$ID, "character")
  expect_type(tT$t, "double")
  expect_type(tT$Gene.symbol, "character")
  expect_type(tT$adj.P.Val, "double")
  expect_type(tT$B, "double")
  expect_type(tT$Gene.title, "character")
  expect_type(tT$P.Value, "double")
  expect_type(tT$logFC, "double")
  expect_type(tT$Gene.ID, "character")

  # Non-Interactive Histogram
  fig <- nonInteractiveHistogramPlot(fit2, adjustment)
  expect_type(fig, "list")
  expect_type(fig$breaks, "double")
  expect_type(fig$counts, "integer")
  expect_type(fig$density, "double")
  expect_type(fig$mids, "double")
  expect_type(fig$xname, "character")
  expect_type(fig$equidist, "logical")

  # Interactive Histogram
  fig <- interactiveHistogramPlot(fit2, adjustment)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Summarize test results as "up", "down" or "not expressed"
  dT <- calculateDifferentialGeneExpressionSummary(fit2, adjustment, significanceLevelCutOff)
  expect_type(dT, 'double')
  expect_equal(ncol(dT),1)
  expect_equal(nrow(dT),897)

  # Non-Interactive Venn diagram
  fig <- nonInteractiveVennDiagramPlot(dT)

  # Non-Interactive Q-Q plot
  fig <- nonInteractiveQQPlot(fit2)
  expect_type(fig, 'list')
  expect_type(fig$y,"double")
  expect_type(fig$x,"double")

  # Interactive Q-Q plot
  ct <- 1
  fig <- interactiveQQPlot(fit2, dT, ct)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # Non-Interactive volcano plot (log P-value vs log fold change)
  fig <- nonInteractiveVolcanoPlot(fit2, dT, ct)

  # Interactive volcano plot (log P-value vs log fold change)
  fig <- interactiveVolcanoPlot(fit2, dT, ct)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')

  # MD plot (log fold change vs mean log expression)
  fig <- noninteractiveMeanDifferencePlot(fit2, dT, ct)

  # Plot Interactive Mean Difference of fit 2 data
  fig <- interactiveMeanDifferencePlot(fit2, dT, ct)
  fig
  expect_type(fig, 'list')
  expect_type(fig$elementId, 'NULL')
  expect_type(fig$height, 'NULL')
  expect_type(fig$width, 'NULL')
  expect_type(fig$x, 'list')
  expect_type(fig$sizingPolicy, 'list')
  expect_type(fig$dependencies, 'list')
  expect_type(fig$preRenderHook, 'closure')
  expect_type(fig$jsHooks, 'list')
})
