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
  geoAccessionCode <- "GSE18380"
  allGset <- getGeoObject(geoAccessionCode)
  ed <- experimentData(allGset[[1]])
  expect_equal(pubMedIds(ed), "19943900")
  ei <- expinfo(ed)
  expect_equal(ei[1], "Alan,D.,Grossman", ignore_attr = TRUE)
  expect_equal(ei[2], "", ignore_attr = TRUE) #lab
  expect_equal(ei[3], "clee2@mit.edu, adg@mit.edu", ignore_attr = TRUE)
  expect_equal(ei[4], "Autonomous plasmid-like replication of a conjugative transposon", ignore_attr = TRUE)
  expect_equal(ei[5], "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18380", ignore_attr = TRUE)

  # Extract platforms
  platforms <- extractPlatforms(allGset)
  platform <- platforms[1]
  expect_type(platforms, 'character')
  expect_type(platform, 'character')
  expect_equal(platform,"GPL4694")

  # Extract the GEO2R data from the specified platform
  gsetData <- extractPlatformGset(allGset, platform)
  expect_type(gsetData, 'S4')
  expect_s4_class(gsetData, 'ExpressionSet')
  expect_equal(nrow(pData(gsetData)),12)
  expect_equal(nrow(fData(gsetData)),4624)

  # Extract the experiment information
  experimentInformation <- extractExperimentInformation(gsetData)
  expect_type(experimentInformation, 'S4')
  expect_s4_class(experimentInformation, 'MIAME')
  expect_equal(experimentInformation@name,"Alan,D.,Grossman")
  expect_equal(experimentInformation@lab, "")
  expect_equal(experimentInformation@contact, "clee2@mit.edu, adg@mit.edu")
  expect_equal(experimentInformation@title, "Autonomous plasmid-like replication of a conjugative transposon")
  expect_equal(nchar(experimentInformation@title), 63)
  expect_equal(experimentInformation@url, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18380")
  expect_equal(experimentInformation@pubMedIds, "19943900")

  # Extract Sample Information
  sampleInfo <- extractSampleInformation(gsetData)
  expect_type(sampleInfo, 'list')
  expect_equal(nrow(sampleInfo),12)
  expect_equal(ncol(sampleInfo),45)

  # Extract expression data
  expressionData <- extractExpressionData(gsetData)
  expect_type(expressionData, 'double')
  expect_equal(ncol(expressionData),12)
  expect_equal(nrow(expressionData),4624)

  # Get column Details
  columnInfo <- extractSampleDetails(gsetData)
  expect_type(columnInfo, 'list')
  expect_equal(ncol(columnInfo),4)
  expect_equal(nrow(columnInfo),12)

  # Is log transformation auto applied
  autoLogInformation <- calculateAutoLogTransformApplication(expressionData)
  expect_type(autoLogInformation, 'character')
  expect_equal(autoLogInformation,"The auto-detect option did not apply log transformation.")

  # Get a list of all the columns
  columns <- extractSampleNames(expressionData)
  expect_type(columns, 'character')
  expect_equal(columns[1],"GSM455528")

  # Apply log transformation to expression data if necessary
  dataInput <- calculateLogTransformation(expressionData, logTransformation)
  expect_type(dataInput, 'double')
  expect_equal(ncol(dataInput),12)
  expect_equal(nrow(dataInput),4624)
  expect_equal(dataInput[1,1],-0.2039365)

  # Perform KNN transformation on log expression data if necessary
  knnDataInput <- calculateKnnImpute(dataInput, "Yes")
  expect_type(knnDataInput, 'double')
  expect_equal(ncol(knnDataInput),12)
  expect_equal(nrow(knnDataInput),4072)
  expect_equal(knnDataInput[1,1],-0.2039365)

  # Get a list of all the columns in the KNN output
  knnColumns <- extractSampleNames(knnDataInput)

  # Get knn output column Details
  knnColumnInfo <- extractSampleDetails(gsetData)
  knnColumnInfo <- knnColumnInfo[knnColumns,]

  # Remove all incomplete rows
  naOmitInput <- calculateNaOmit(knnDataInput)
  expect_type(naOmitInput, 'double')
  expect_equal(ncol(naOmitInput),12)
  expect_equal(nrow(naOmitInput),4072)
  expect_equal(naOmitInput[1,1],-0.2039365)

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
  expect_equal(extractedExperimentInformation[1], "<b> Paper Title: </b> <p> Autonomous plasmid-like replication of a conjugative transposon </p><b> Author's Name: </b> <p> Alan,D.,Grossman </p><b> Laboratory: </b> <p>  </p><b> Contact Information: </b> <p> clee2@mit.edu, adg@mit.edu </p><b> Paper URL: </b> <p> https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18380 </p><b> PubMed ID: </b> <p> 19943900 </p><b> Abstract: </b> <p> Integrative and conjugative elements (ICEs), a.k.a., conjugative transposons, are mobile genetic elements involved in many biological processes, including the spread of antibiotic resistance. Unlike conjugative plasmids that are extra-chromosomal and replicate autonomously, ICEs are integrated in the chromosome and replicate passively during chromosomal replication. It is generally thought that ICEs do not replicate autonomously. We found that when induced, Bacillus subtilis ICEBs1 replicates as a plasmid. The ICEBs1 origin of transfer (oriT) served as the origin of replication and the conjugal DNA relaxase served as the replication initiation protein. Autonomous replication of ICEBs1 conferred genetic stability to the excised element, but was not required for mating. The B. subtilis helicase PcrA that mediates unwinding and replication of Gram-positive rolling circle replicating plasmids was required for ICEBs1 replication and mating. Nicking of oriT by the relaxase and unwinding by PcrA likely directs transfer of a single-strand of ICEBs1 into recipient cells.\n\nThis SuperSeries is composed of the SubSeries listed below. </p>")

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
  expect_equal(column2[1], "GSM455783")
  expect_equal(column2[2], "GSM455784")
  expect_equal(column2[3], "GSM455785")
  expect_equal(column2[4], "GSM455786")
  expect_equal(column2[5], "GSM455787")
  expect_equal(column2[6], "NA")
  expect_equal(length(column2), 5)

  # Calculate gsms
  gsms <- calculateEachGroupsSamples(columnNames,group1, group2)
  expect_type(gsms, "character")
  expect_equal(gsms, "000000011111")
  expect_equal(nchar(gsms), 12)

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
  expect_type(tT$SPOT_ID, "character")
  expect_type(tT$GB_ACC, "character")
  expect_type(tT$adj.P.Val, "double")
  expect_type(tT$B, "double")
  expect_type(tT$RANGE_GB, "character")
  expect_type(tT$P.Value, "double")
  expect_type(tT$logFC, "double")
  expect_type(tT$RANGE_START, "double")

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
  expect_equal(nrow(dT),3234)

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