library(GEOexplorer)
context("Without Missing Values")

test_that("Microarray GSE without missing values is handled correctly by all
          functions",
          {
            input <- NULL
            all <- NULL
            input$logTransformation <- "Auto-Detect"
            input$knnTransformation <- "Yes"
            input$knn <- 2
            input$pValueAdjustment <- 
              "Benjamini & Hochberg (False discovery rate)"
            input$limmaPrecisionWeights <- "Yes"
            input$forceNormalization <- "Yes"
            input$platformAnnotation <- "NCBI generated"
            input$significanceLevelCutOff <- 0.05
            input$dataSource <- "GEO"
            all$typeOfData <- "Microarray"
            input$dataSetType <- "Single"

            # Get the GEO data for all platforms
            input$geoAccessionCode <- "GSE18388"
            allGset <- getGeoObject(input$geoAccessionCode)
            ed <- experimentData(allGset[[1]])
            expect_equal(pubMedIds(ed), "20213684")
            ei <- expinfo(ed)
            expect_equal(ei[1], "Ty,W,Lebsack", ignore_attr = TRUE)
            expect_equal(ei[2], "", ignore_attr = TRUE) #lab
            expect_equal(ei[3], "lebsack@email.arizona.edu",
                         ignore_attr = TRUE)
            expect_equal(
              ei[4],
              "Microarray Analysis of Space-flown Murine Thymus Tissue",
              ignore_attr = TRUE)
            expect_equal(
              ei[5],
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18388",
              ignore_attr = TRUE
            )

            # Extract platforms
            platforms <- extractPlatforms(allGset)
            platform <- platforms[1]
            expect_type(platforms, 'character')
            expect_type(platform, 'character')
            expect_equal(platform, "GPL6246")

            # Extract the GEO2R data from the specified platform
            all$gsetData <- extractPlatformGset(allGset, platform)
            expect_type(all$gsetData, 'S4')
            expect_s4_class(all$gsetData, 'ExpressionSet')
            expect_equal(nrow(pData(all$gsetData)), 8)
            expect_equal(nrow(fData(all$gsetData)), 35557)

            # Extract the experiment information
            experimentInformation <- extractExperimentInformation(all$gsetData)
            expect_type(experimentInformation, 'S4')
            expect_s4_class(experimentInformation, 'MIAME')
            expect_equal(experimentInformation@name, "Ty,W,Lebsack")
            expect_equal(experimentInformation@lab, "")
            expect_equal(experimentInformation@contact,
                         "lebsack@email.arizona.edu")
            expect_equal(
              experimentInformation@title,
              "Microarray Analysis of Space-flown Murine Thymus Tissue"
            )
            expect_equal(nchar(experimentInformation@title), 55)
            expect_equal(
              experimentInformation@url,
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18388"
            )
            expect_equal(experimentInformation@pubMedIds, "20213684")

            # Extract Sample Information
            sampleInfo <- extractSampleInformation(all$gsetData)
            expect_type(sampleInfo, 'list')
            expect_equal(nrow(sampleInfo), 8)
            expect_equal(ncol(sampleInfo), 35)

            # Extract expression data
            expressionData <- extractExpressionData(all$gsetData)
            expect_type(expressionData, 'double')
            expect_equal(ncol(expressionData), 8)
            expect_equal(nrow(expressionData), 35557)

            # Get column Details
            columnInfo <- extractSampleDetails(all$gsetData)
            expect_type(columnInfo, 'list')
            expect_equal(ncol(columnInfo), 5)
            expect_equal(nrow(columnInfo), 8)

            # Is log transformation auto applied
            autoLogInformation <-
              calculateAutoLogTransformApplication(expressionData)
            expect_type(autoLogInformation, 'character')
            expect_equal(
              autoLogInformation,
              "The auto-detect option did not apply log transformation.")

            # Get a list of all the columns
            columns <- extractSampleNames(expressionData)
            expect_type(columns, 'character')
            expect_equal(columns[1], "GSM458594")

            # Apply log transformation to expression data if necessary
            dataInput <-
              calculateLogTransformation(expressionData, 
                                         input$logTransformation)
            expect_type(dataInput, 'double')
            expect_equal(ncol(dataInput), 8)
            expect_equal(nrow(dataInput), 35557)
            expect_equal(dataInput[1, 1], 11.711505)

            # Perform input$knn transformation on log expression data if 
            # necessary
            all$knnDataInput <- calculateKnnImpute(dataInput, "Yes")
            expect_type(all$knnDataInput, 'double')
            expect_equal(ncol(all$knnDataInput), 8)
            expect_equal(nrow(all$knnDataInput), 35557)
            expect_equal(all$knnDataInput[1, 1], 11.711505)

            # Get a list of all the columns in the input$knn output
            knnColumns <- extractSampleNames(all$knnDataInput)

            # Get input$knn output column Details
            knnColumnInfo <- extractSampleDetails(all$gsetData)
            knnColumnInfo <- knnColumnInfo[knnColumns, ]

            # Remove all incomplete rows
            naOmitInput <- calculateNaOmit(all$knnDataInput)
            expect_type(naOmitInput, 'double')
            expect_equal(ncol(naOmitInput), 8)
            expect_equal(nrow(naOmitInput), 35557)
            expect_equal(naOmitInput[1, 1], 11.711505)

            # Perform Princomp PCA analysis on input$knn transformation
            # expression data
            pcaPrincompDataInput <- calculatePrincompPca(naOmitInput)
            expect_type(pcaPrincompDataInput, 'list')
            expect_s3_class(pcaPrincompDataInput, 'princomp')

            # Perform Prcomp PCA analysis on input$knn transformation 
            # expression data
            pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
            expect_type(pcaPrcompDataInput, 'list')
            expect_s3_class(pcaPrcompDataInput, 'prcomp')

            # Extract Experiment Information
            extractedExperimentInformation <- convertExperimentInformation(
              experimentInformation)
            expect_type(extractedExperimentInformation, 'character')

            # Non-Interactive Box-and-Whisker Plot
            fig <- nonInteractiveBoxAndWhiskerPlot(
              ex = all$knnDataInput)

            expect_type(fig, 'list')
            expect_type(fig$stats, 'double')
            expect_type(fig$n, 'double')
            expect_type(fig$conf, 'double')
            expect_type(fig$out, 'double')
            expect_type(fig$group, 'double')
            expect_type(fig$names, 'character')

            # Interactive Box-and-Whisker Plot
            fig <-
              interactiveBoxAndWhiskerPlot(
                all$knnDataInput)
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
            fig <- nonInteractiveDensityPlot(
              ex = naOmitInput)

            expect_type(fig, 'list')
            expect_type(fig$X, 'double')
            expect_type(fig$Y, 'double')

            # Interactive Density Plot
            fig <-
              interactiveDensityPlot(naOmitInput)
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
            fig <-
              interactiveThreeDDensityPlot(naOmitInput)
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
            fig <- interactiveUmapPlot(naOmitInput, input$knn)
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
            fig <-
              interactiveMeanVariancePlot(naOmitInput, all$gsetData)
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
            fig <- interactivePrincompPcaScreePlot(pcaPrincompDataInput)
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
            fig <- interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput,
                                                         all$gsetData)
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
            fig <- interactivePrincompPcaVariablesPlot(pcaPrincompDataInput)
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
            fig <-
              interactivePrcompPcaScreePlot(pcaPrcompDataInput)
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
            fig <- interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
                                                       all$gsetData)
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
            fig <- interactivePrcompPcaVariablesPlot(pcaPrcompDataInput)
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
            fig <- nonInteractiveUmapPlot(naOmitInput, input$knn)
            expect_type(fig, 'list')
            expect_type(fig$x, 'double')
            expect_type(fig$y, 'double')

            # Non-Interactive Mean Variance Plot
            fig <-
              nonInteractiveMeanVariancePlot(naOmitInput)

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

            # Non-Interactive Princomp PCA Individual and Variables Bilot
            fig <- nonInteractivePcaBiplotPlot(pcaPrincompDataInput)
            fig

            # Differential gene expression analysis functions
            # Get column names
            columnNames <- extractSampleNames(expressionData)

            # Define Groups
            numberOfColumns <- length(columnNames)
            numberOfColumns <- numberOfColumns + 1
            halfNumberOfColumns <- ceiling(numberOfColumns / 2)
            i <- 0

            group1 <- c()
            group2 <- c()

            for (name in columnNames) {
              if (i < halfNumberOfColumns) {
                group1 <- c(group1, name)
                i <- i + 1
              } else {
                group2 <- c(group2, name)
                i <- i + 1
              }
            }

            # Select columns in group2
            column2 <- calculateExclusiveColumns(columnNames, group1)
            expect_type(column2, "character")
            expect_equal(column2[1], "GSM458599")
            expect_equal(column2[2], "GSM458600")
            expect_equal(column2[3], "GSM458601")
            expect_equal(column2[4], "NA")
            expect_equal(length(column2), 3)

            # Calculate all$gsms
            all$gsms <- calculateEachGroupsSamples(columnNames, group1, group2)
            expect_type(all$gsms, "character")
            expect_equal(all$gsms, "00000111")
            expect_equal(nchar(all$gsms), 8)
            
            all$gsms <- "00001111"

            # Convert P value adjustment
            input$pValueAdjustment <- 
              "Benjamini & Hochberg (False discovery rate)"
            adjustment <- convertAdjustment(input$pValueAdjustment)
            expect_type(adjustment, "character")
            expect_equal(adjustment, "fdr")

            input$pValueAdjustment <- "Benjamini & Yekutieli"
            adjustment <- convertAdjustment(input$pValueAdjustment)
            expect_type(adjustment, "character")
            expect_equal(adjustment, "BY")

            input$pValueAdjustment <- "Bonferroni"
            adjustment <- convertAdjustment(input$pValueAdjustment)
            expect_type(adjustment, "character")
            expect_equal(adjustment, "bonferroni")

            input$pValueAdjustment <- "Holm"
            adjustment <- convertAdjustment(input$pValueAdjustment)
            expect_type(adjustment, "character")
            expect_equal(adjustment, "holm")

            input$pValueAdjustment <- "None"
            adjustment <- convertAdjustment(input$pValueAdjustment)
            expect_type(adjustment, "character")
            expect_equal(adjustment, "none")

            adjustment <- convertAdjustment(input$pValueAdjustment)

            # Get fit 2
            results <- calculateDifferentialGeneExpression(all$gsms,input,
                                                           all)
            expect_type(results$fit2, "list")
            expect_type(results$fit2$coefficients, "double")
            expect_type(results$fit2$sigma, "double")
            expect_type(results$fit2$cov.coefficients, "double")
            expect_type(results$fit2$rank, "integer")
            expect_type(results$fit2$Amean, "double")
            expect_type(results$fit2$design, "double")
            expect_type(results$fit2$df.prior, "double")
            expect_type(results$fit2$var.prior, "double")
            expect_type(results$fit2$s2.post, "double")
            expect_type(results$fit2$df.total, "double")
            expect_type(results$fit2$lods, "double")
            expect_type(results$fit2$F.p.value, "double")
            expect_type(results$fit2$stdev.unscaled, "double")
            expect_type(results$fit2$df.residual, "double")
            expect_type(results$fit2$pivot, "integer")
            expect_type(results$fit2$genes, "list")
            expect_type(results$fit2$method, "character")
            expect_type(results$fit2$contrasts, "double")
            expect_type(results$fit2$s2.prior, "double")
            expect_type(results$fit2$proportion, "double")
            expect_type(results$fit2$t, "double")
            expect_type(results$fit2$p.value, "double")
            expect_type(results$fit2$F, "double")

            # Print Top deferentially expressed genes
            all$tT <- calculateTopDifferentiallyExpressedGenes(results$fit2,
                                                           adjustment)
            expect_type(all$tT, "list")
            expect_type(all$tT$ID, "character")
            expect_type(all$tT$t, "double")
            expect_type(all$tT$Gene.symbol, "character")
            expect_type(all$tT$adj.P.Val, "double")
            expect_type(all$tT$B, "double")
            expect_type(all$tT$Gene.title, "character")
            expect_type(all$tT$P.Value, "double")
            expect_type(all$tT$logFC, "double")
            expect_type(all$tT$Gene.ID, "character")

            # Non-Interactive Histogram
            fig <- nonInteractiveHistogramPlot(results$fit2, adjustment)
            expect_type(fig, "list")
            expect_type(fig$breaks, "double")
            expect_type(fig$counts, "integer")
            expect_type(fig$density, "double")
            expect_type(fig$mids, "double")
            expect_type(fig$xname, "character")
            expect_type(fig$equidist, "logical")

            # Interactive Histogram
            fig <- interactiveHistogramPlot(results$fit2, adjustment)
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
            all$dT <- calculateDifferentialGeneExpressionSummary(
              results$fit2, adjustment, input$significanceLevelCutOff)
            expect_type(all$dT, 'double')
            expect_equal(ncol(all$dT), 1)
            expect_equal(nrow(all$dT), 35557)

            # Non-Interactive Venn diagram
            fig <- nonInteractiveVennDiagramPlot(all$dT)

            # Non-Interactive Q-Q plot
            fig <- nonInteractiveQQPlot(results$fit2)
            expect_type(fig, 'list')
            expect_type(fig$y, "double")
            expect_type(fig$x, "double")

            # Interactive Q-Q plot
            ct <- 1
            fig <- interactiveQQPlot(results$fit2, all$dT, ct)
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
            fig <- nonInteractiveVolcanoPlot(results$fit2, all$dT, ct)

            # Interactive volcano plot (log P-value vs log fold change)
            fig <- interactiveVolcanoPlot(results$fit2, all$dT, ct)
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
            fig <- noninteractiveMeanDifferencePlot(results$fit2, all$dT, ct)

            # Plot Interactive Mean Difference of fit 2 data
            fig <- interactiveMeanDifferencePlot(results$fit2, all$dT, ct)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')

            # Plot Interactive Heatmap Plot
            numberOfGenes <- 20
            fig <- interactiveDGEHeatMapPlot(results$ex,
                                             input$limmaPrecisionWeights,
                                             numberOfGenes, all$tT)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            output <- NULL
            session <- NULL
            errorChecks <- NULL
            
            differentiallyExressedGeneAnnotation <- 
              createGeneAnnotationTable(input, output, session, errorChecks, 
                                        all)
            expect_type(differentiallyExressedGeneAnnotation, 'list')
            
            
            x <- differentiallyExressedGeneAnnotation[,c(4,ncol(
              differentiallyExressedGeneAnnotation))]
            expect_type(x, 'list')
            
            colnames(x) <- c("Gene.symbol", "Group1-Group2")
            expect_type(x, 'list')
            
            # Extract differentially expressed gene symbols
            differemtiallyExpressedGeneSymbols <-
              extractGeneSymbols(x, "Gene.symbol")
            expect_type(differemtiallyExpressedGeneSymbols, 'character')
            
            # enrich Differentially Expressed Genes
            enrichedDifferentiallyExpressedGenes <-
              enrichGenes(differemtiallyExpressedGeneSymbols, 
                          "GO_Biological_Process_2015")
            expect_type(enrichedDifferentiallyExpressedGenes, 'list')
            
            
            enrichedDifferentiallyExpressedGenes <-calculateLogPValue(
              enrichedDifferentiallyExpressedGenes)
            expect_type(enrichedDifferentiallyExpressedGenes, 'list')
            
            
            enrichedDifferentiallyExpressedGenes <- calculateOverlapFractions(
              enrichedDifferentiallyExpressedGenes)
            expect_type(enrichedDifferentiallyExpressedGenes, 'list')
            
            
            columnToSort <- "P.value"
            recordsToDisplay <- 20
            sortDecreasingly <- TRUE
            
            fig <- interactiveGeneEnrichmentVolcanoPlot(
              enrichedDifferentiallyExpressedGenes)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            fig <- interactiveGeneEnrichmentManhattanPlot(
              enrichedDifferentiallyExpressedGenes, columnToSort)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            sortedEnrichedDifferentiallyExpressedGenes <-
              sortGeneEnrichmentTable(enrichedDifferentiallyExpressedGenes, 
                                      columnToSort,
                                      sortDecreasingly)
            expect_type(sortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            
            topSortedEnrichedDifferentiallyExpressedGenes <-
              selectTopGeneEnrichmentRecords(
                enrichedDifferentiallyExpressedGenes, recordsToDisplay)
            expect_type(topSortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            fig <- interactiveGeneEnrichmentBarPlot(
              topSortedEnrichedDifferentiallyExpressedGenes, columnToSort)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            # Extract differentially expressed genes
            differentiallyExpressedGenes <- 
              extractDifferenitallyExpressedGenes(all$tT, all$dT)
            expect_type(differentiallyExpressedGenes, 'list')
            
            
            # Extract differentially expressed gene symbols
            differemtiallyExpressedGeneSymbols <-
              extractGeneSymbols(differentiallyExpressedGenes, "Gene.symbol")
            expect_type(differemtiallyExpressedGeneSymbols, 'character')
            
            
            # enrich Differentially Expressed Genes
            enrichedDifferentiallyExpressedGenes <-
              enrichGenes(differemtiallyExpressedGeneSymbols, 
                          "GO_Biological_Process_2015")
            expect_type(enrichedDifferentiallyExpressedGenes, 'list')
            
            # Extract Upregulated genes
            upregulatedGenes <- extractUpregulatedGenes(
              differentiallyExpressedGenes)
            expect_type(upregulatedGenes, 'list')
            
            # Extract upregulated gene symbols
            upregulatedGenesGeneSymbols <-
              extractGeneSymbols(upregulatedGenes, "Gene.symbol")
            expect_type(upregulatedGenesGeneSymbols, 'character')
            
            # enrich upregulated Genes
            enrichedUpregulatedGenes <-
              enrichGenes(upregulatedGenesGeneSymbols, 
                          "GO_Biological_Process_2015")
            expect_type(enrichedUpregulatedGenes, 'list')
            
            enrichedUpregulatedGenes <-calculateLogPValue(
              enrichedUpregulatedGenes)
            expect_type(enrichedUpregulatedGenes, 'list')
            
            enrichedUpregulatedGenes <- calculateOverlapFractions(
              enrichedUpregulatedGenes)
            expect_type(enrichedUpregulatedGenes, 'list')
            
            # enrich downregulated Genes
            # Extract downregulated genes
            downregulatedGenes <- extractdowregulatedGenes(
              differentiallyExpressedGenes)
            expect_type(downregulatedGenes, 'list')
            
            
            # Extract downregulated gene symbols
            downregulatedGenesGeneSymbols <-
              extractGeneSymbols(downregulatedGenes, "Gene.symbol")
            expect_type(downregulatedGenesGeneSymbols, 'character')
            
            enrichedDownregulatedGenes <-
              enrichGenes(downregulatedGenesGeneSymbols, 
                          "GO_Biological_Process_2015")
            expect_type(enrichedDownregulatedGenes, 'list')
            
            enrichedDownregulatedGenes <-calculateLogPValue(
              enrichedDownregulatedGenes)
            expect_type(enrichedDownregulatedGenes, 'list')
            
            
            enrichedDownregulatedGenes <- calculateOverlapFractions(
              enrichedDownregulatedGenes)
            expect_type(enrichedDownregulatedGenes, 'list')
            
            fig <- interactiveGeneEnrichmentVolcanoPlot(
              enrichedUpregulatedGenes)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            fig <- interactiveGeneEnrichmentManhattanPlot(
              enrichedUpregulatedGenes, columnToSort)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            sortedEnrichedDifferentiallyExpressedGenes <-
              sortGeneEnrichmentTable(enrichedUpregulatedGenes, columnToSort,
                                      sortDecreasingly)
            expect_type(sortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            topSortedEnrichedDifferentiallyExpressedGenes <-
              selectTopGeneEnrichmentRecords(
                sortedEnrichedDifferentiallyExpressedGenes,recordsToDisplay)
            expect_type(topSortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            fig <- interactiveGeneEnrichmentVolcanoPlot(
              enrichedDownregulatedGenes)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            fig <- interactiveGeneEnrichmentManhattanPlot(
              enrichedDownregulatedGenes, columnToSort)
            expect_type(fig, 'list')
            expect_type(fig$elementId, 'NULL')
            expect_type(fig$height, 'NULL')
            expect_type(fig$width, 'NULL')
            expect_type(fig$x, 'list')
            expect_type(fig$sizingPolicy, 'list')
            expect_type(fig$dependencies, 'list')
            expect_type(fig$preRenderHook, 'closure')
            expect_type(fig$jsHooks, 'list')
            
            sortedEnrichedDifferentiallyExpressedGenes <-
              sortGeneEnrichmentTable(enrichedDownregulatedGenes, columnToSort,
                                      sortDecreasingly)
            expect_type(sortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            topSortedEnrichedDifferentiallyExpressedGenes <-
              selectTopGeneEnrichmentRecords(enrichedDownregulatedGenes,
                                             recordsToDisplay)
            expect_type(topSortedEnrichedDifferentiallyExpressedGenes, 'list')
            
            
            fig <- interactiveGeneEnrichmentBarPlot(enrichedDownregulatedGenes, 
                                                    columnToSort)
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
