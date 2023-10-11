library(GEOexplorer)
context("Non-Log Values")

test_that("Microarray GSE with non-log values is handled correctly
          by all functions",
          {
            # Input Values
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

            # Get the GEO data for all all$platforms
            input$geoAccessionCode <- "GSE2"
            all$allGset <- getGeoObject(input$geoAccessionCode)
            all$ed <- experimentData(all$allGset[[1]])
            expect_equal(pubMedIds(all$ed), "")
            all$ei <- expinfo(all$ed)
            expect_equal(all$ei[1], "Yoshihiro,,Kagami", ignore_attr = TRUE)
            expect_equal(all$ei[2], "", ignore_attr = TRUE) #lab
            expect_equal(all$ei[3], "ykagami@brain.riken.go.jp",
                         ignore_attr = TRUE)
            expect_equal(all$ei[4], "Cerebellar development", 
                         ignore_attr = TRUE)
            expect_equal(
              all$ei[5],
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2",
              ignore_attr = TRUE)

            # Extract all$platforms
            all$platforms <- extractPlatforms(all$allGset)
            all$platform <- all$platforms[1]
            expect_type(all$platforms, 'character')
            expect_type(all$platform, 'character')
            expect_equal(all$platform, "GPL8")

            # Extract the GEO2R data from the specified all$platform
            all$gsetData <- extractPlatformGset(all$allGset, all$platform)
            expect_type(all$gsetData, 'S4')
            expect_s4_class(all$gsetData, 'ExpressionSet')
            expect_equal(nrow(pData(all$gsetData)), 5)
            expect_equal(nrow(fData(all$gsetData)), 897)

            # Extract the experiment information
            all$experimentInformation <-
              extractExperimentInformation(all$gsetData)
            expect_type(all$experimentInformation, 'S4')
            expect_s4_class(all$experimentInformation, 'MIAME')
            expect_equal(all$experimentInformation@name, "Yoshihiro,,Kagami")
            expect_equal(all$experimentInformation@lab, "")
            expect_equal(all$experimentInformation@contact,
                         "ykagami@brain.riken.go.jp")
            expect_equal(all$experimentInformation@title, 
                         "Cerebellar development")
            expect_equal(nchar(all$experimentInformation@title), 22)
            expect_equal(
              all$experimentInformation@url,
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2"
            )
            expect_equal(all$experimentInformation@pubMedIds, "")

            # Extract Sample Information
            all$sampleInfo <- extractSampleInformation(all$gsetData)
            expect_type(all$sampleInfo, 'list')
            expect_equal(nrow(all$sampleInfo), 5)
            expect_equal(ncol(all$sampleInfo), 30)

            # Extract expression data
            all$expressionData <- extractExpressionData(all$gsetData)
            expect_type(all$expressionData, 'double')
            expect_equal(ncol(all$expressionData), 5)
            expect_equal(nrow(all$expressionData), 897)

            # Get column Details
            all$columnInfo <- extractSampleDetails(all$gsetData)
            expect_type(all$columnInfo, 'list')
            expect_equal(ncol(all$columnInfo), 5)
            expect_equal(nrow(all$columnInfo), 5)

            # Is log transformation auto applied
            all$autoLogInformation <-
              calculateAutoLogTransformApplication(all$expressionData)
            expect_type(all$autoLogInformation, 'character')
            expect_equal(all$autoLogInformation,
                         "The auto-detect option applied log transformation.")

            # Get a list of all the all$columns
            all$columns <- extractSampleNames(all$expressionData)
            expect_type(all$columns, 'character')
            expect_equal(all$columns[1], "GSM50")

            # Apply log transformation to expression data if necessary
            all$dataInput <-
              calculateLogTransformation(all$expressionData, 
                                         input$logTransformation)
            expect_type(all$dataInput, 'double')
            expect_equal(ncol(all$dataInput), 5)
            expect_equal(nrow(all$dataInput), 897)
            expect_equal(all$dataInput[1, 1], 6.5833085)

            # Perform input$knn transformation on log expression data 
            # if necessary
            all$knnDataInput <- calculateKnnImpute(all$dataInput, "Yes")
            expect_type(all$knnDataInput, 'double')
            expect_equal(ncol(all$knnDataInput), 5)
            expect_equal(nrow(all$knnDataInput), 897)
            expect_equal(all$knnDataInput[1, 1], 6.5833085)

            # Get a list of all the all$columns in the input$knn output
            all$knnColumns <- extractSampleNames(all$knnDataInput)

            # Get input$knn output column Details
            all$knnColumnInfo <- extractSampleDetails(all$gsetData)
            all$knnColumnInfo <- all$knnColumnInfo[all$knnColumns,]

            # Remove all incomplete rows
            all$naOmitInput <- calculateNaOmit(all$knnDataInput)
            expect_type(all$naOmitInput, 'double')
            expect_equal(ncol(all$naOmitInput), 5)
            expect_equal(nrow(all$naOmitInput), 897)
            expect_equal(all$naOmitInput[1, 1], 6.5833085)

            # Perform Princomp PCA analysis on input$knn transformation
            # expression data
            all$pcaPrincompDataInput <-
              calculatePrincompPca(all$naOmitInput)
            expect_type(all$pcaPrincompDataInput, 'list')
            expect_s3_class(all$pcaPrincompDataInput, 'princomp')

            # Perform Prcomp PCA analysis on input$knn transformation 
            # expression data
            all$pcaPrcompDataInput <- calculatePrcompPca(all$naOmitInput)
            expect_type(all$pcaPrcompDataInput, 'list')
            expect_s3_class(all$pcaPrcompDataInput, 'prcomp')

            # Extract Experiment Information
            extractedExperimentInformation <-
              convertExperimentInformation(all$experimentInformation)
            expect_type(extractedExperimentInformation, 'character')
            expect_equal(nchar(extractedExperimentInformation[1]), 1216)

            # Non-Interactive Box-and-Whisker Plot
            fig <-
              nonInteractiveBoxAndWhiskerPlot(ex = all$knnDataInput)
            expect_type(fig, 'list')
            expect_type(fig$stats, 'double')
            expect_type(fig$n, 'double')
            expect_type(fig$conf, 'double')
            expect_type(fig$out, 'double')
            expect_type(fig$group, 'double')
            expect_type(fig$names, 'character')

            # Interactive Box-and-Whisker Plot
            fig <-
              interactiveBoxAndWhiskerPlot(all$knnDataInput)
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
            fig <-
              nonInteractiveDensityPlot(ex = all$naOmitInput)
            expect_type(fig, 'list')
            expect_type(fig$X, 'double')
            expect_type(fig$Y, 'double')

            # Interactive Density Plot
            fig <-
              interactiveDensityPlot(all$naOmitInput)
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
              interactiveThreeDDensityPlot(all$naOmitInput)
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
            fig <-
              interactiveUmapPlot(all$naOmitInput, input$knn)
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
              interactiveMeanVariancePlot(all$naOmitInput, all$gsetData)
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
            fig <-
              interactivePrincompPcaScreePlot(all$pcaPrincompDataInput)
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
            fig <-
              interactivePrincompPcaIndividualsPlot(all$pcaPrincompDataInput,
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
            fig <-
              interactivePrincompPcaVariablesPlot(all$pcaPrincompDataInput)
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
              interactivePrcompPcaScreePlot(all$pcaPrcompDataInput)
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
            fig <-
              interactivePrcompPcaIndividualsPlot(all$pcaPrcompDataInput, 
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
            fig <-
              interactivePrcompPcaVariablesPlot(all$pcaPrcompDataInput)
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
            fig <- interactiveHeatMapPlot(all$naOmitInput)
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
            fig <-
              nonInteractiveUmapPlot(all$naOmitInput, input$knn)
            expect_type(fig, 'list')
            expect_type(fig$x, 'double')
            expect_type(fig$y, 'double')

            # Non-Interactive Mean Variance Plot
            fig <-
              nonInteractiveMeanVariancePlot(all$naOmitInput)

            # Non-Interactive Princomp PCA Scree Plot
            fig <- nonInteractivePcaScreePlot(all$pcaPrincompDataInput)
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
            fig <-
              nonInteractivePcaIndividualsPlot(all$pcaPrincompDataInput)
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
            fig <-
              nonInteractivePcaVariablesPlot(all$pcaPrincompDataInput)
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
            fig <- nonInteractivePcaBiplotPlot(all$pcaPrincompDataInput)
            fig

            # Differential gene expression analysis functions
            # Get column names
            columnNames <- extractSampleNames(all$expressionData)

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

            # Select all$columns in group2
            column2 <-
              calculateExclusiveColumns(columnNames, group1)
            expect_type(column2, "character")
            expect_equal(column2[1], "GSM53")
            expect_equal(column2[2], "GSM54")
            # expect_equal(column2[3], "NA")
            expect_equal(length(column2), 2)

            # Calculate all$gsms
            all$gsms <-
              calculateEachGroupsSamples(columnNames, group1, group2)
            expect_type(all$gsms, "character")
            expect_equal(all$gsms, "00011")
            expect_equal(nchar(all$gsms), 5)
            all$gsms <- "11X00"

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
            results <-
              calculateDifferentialGeneExpression(all$gsms,
                                                  input,
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
            expect_equal(nrow(all$dT), 897)

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
            
            # Extract differential expressed gene symbols
            differemtiallyExpressedGeneSymbols <-
              extractGeneSymbols(x, "Gene.symbol")
            expect_type(differemtiallyExpressedGeneSymbols, 'character')
            
            # enrich Differential Expressed Genes
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
            
            head(x)
            
            x[x["Group1-Group2"] == 1, ] 
            
            # Extract Upregulated genes
            upregulatedGenes <- extractUpregulatedGenes(
              x)
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
              x)
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
