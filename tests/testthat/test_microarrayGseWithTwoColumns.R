library(GEOexplorer)
context("With Two Columns")

test_that("Microarray GSE with two columns is handled correctly by
          all functions",
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
            input$geoAccessionCode <- "GSE25743"
            all$allGset <- getGeoObject(input$geoAccessionCode)
            ed <- experimentData(all$allGset[[1]])
            expect_equal(pubMedIds(ed), "21281469")
            ei <- expinfo(ed)
            expect_equal(ei[1], "Gilbert ,,Spizzo", ignore_attr = TRUE)
            expect_equal(ei[2], "", ignore_attr = TRUE) #lab
            expect_equal(ei[3], "", ignore_attr = TRUE)
            expect_equal(ei[4], "EpCAM expressing Hs578T cells",
                         ignore_attr = TRUE)
            expect_equal(
              ei[5],
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25743",
              ignore_attr = TRUE
            )

            # Extract platforms
            platforms <- extractPlatforms(all$allGset)
            platform <- platforms[1]
            expect_type(platforms, 'character')
            expect_type(platform, 'character')
            expect_equal(platform, "GPL570")

            # Extract the GEO2R data from the specified platform
            all$gsetData <- extractPlatformGset(all$allGset, platform)
            expect_type(all$gsetData, 'S4')
            expect_s4_class(all$gsetData, 'ExpressionSet')
            expect_equal(nrow(pData(all$gsetData)), 2)
            expect_equal(nrow(fData(all$gsetData)), 54675)

            # Extract the experiment information
            experimentInformation <-
              extractExperimentInformation(all$gsetData)
            expect_type(experimentInformation, 'S4')
            expect_s4_class(experimentInformation, 'MIAME')
            expect_equal(experimentInformation@name, "Gilbert ,,Spizzo")
            expect_equal(experimentInformation@lab, "")
            expect_equal(experimentInformation@contact, "")
            expect_equal(experimentInformation@title,
            "EpCAM expressing Hs578T cells")
            expect_equal(nchar(experimentInformation@title), 29)
            expect_equal(
              experimentInformation@url,
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25743"
            )
            expect_equal(experimentInformation@pubMedIds, "21281469")

            # Extract Sample Information
            sampleInfo <- extractSampleInformation(all$gsetData)
            expect_type(sampleInfo, 'list')
            expect_equal(nrow(sampleInfo), 2)
            expect_equal(ncol(sampleInfo), 37)

            # Extract expression data
            expressionData <- extractExpressionData(all$gsetData)
            expect_type(expressionData, 'double')
            expect_equal(ncol(expressionData), 2)
            expect_equal(nrow(expressionData), 54675)
            expressionData <- as.matrix(expressionData)
            
            # Get column Details
            columnInfo <- extractSampleDetails(all$gsetData)
            expect_type(columnInfo, 'list')
            expect_equal(ncol(columnInfo), 5)
            expect_equal(nrow(columnInfo), 2)

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
            expect_equal(columns[1], "GSM632607")

            # Apply log transformation to expression data if necessary
            dataInput <-
              calculateLogTransformation(expressionData, 
                                         input$logTransformation)
            expect_type(dataInput, 'double')
            expect_equal(ncol(dataInput), 2)
            expect_equal(nrow(dataInput), 54675)
            expect_equal(dataInput[1, 1], 7.7585824)

            # Perform input$knn transformation on log expression data if 
            # necessary
            expect_error(all$knnDataInput <-
                           calculateKnnImpute(dataInput, "Yes"))
            all$knnDataInput <- calculateKnnImpute(dataInput, "No")
            expect_type(all$knnDataInput, 'double')
            expect_equal(ncol(all$knnDataInput), 2)
            expect_equal(nrow(all$knnDataInput), 54675)
            expect_equal(all$knnDataInput[1, 1], 7.7585824)

            # Get a list of all the columns in the input$knn output
            knnColumns <- extractSampleNames(all$knnDataInput)

            # Get input$knn output column Details
            knnColumnInfo <- extractSampleDetails(all$gsetData)
            knnColumnInfo <- knnColumnInfo[knnColumns,]

            # Remove all incomplete rows
            naOmitInput <- calculateNaOmit(all$knnDataInput)
            expect_type(naOmitInput, 'double')
            expect_equal(ncol(naOmitInput), 2)
            expect_equal(nrow(naOmitInput), 54675)
            expect_equal(naOmitInput[1, 1], 7.7585824)

            # Perform Princomp PCA analysis on input$knn transformation
            # expression data
            pcaPrincompDataInput <-
              calculatePrincompPca(naOmitInput)
            expect_type(pcaPrincompDataInput, 'list')
            expect_s3_class(pcaPrincompDataInput, 'princomp')

            # Perform Prcomp PCA analysis on input$knn transformation 
            # expression data
            pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)
            expect_type(pcaPrcompDataInput, 'list')
            expect_s3_class(pcaPrcompDataInput, 'prcomp')

            # Extract Experiment Information
            extractedExperimentInformation <-
              convertExperimentInformation(experimentInformation)
            expect_type(extractedExperimentInformation, 'character')

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
              nonInteractiveDensityPlot(ex = naOmitInput)
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
            fig <-
              interactiveUmapPlot(naOmitInput, input$knn)
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
            fig <-
              interactivePrincompPcaScreePlot(pcaPrincompDataInput)
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
              interactivePrincompPcaIndividualsPlot(pcaPrincompDataInput,
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
              interactivePrincompPcaVariablesPlot(pcaPrincompDataInput)
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
            fig <-
              interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput, 
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
              interactivePrcompPcaVariablesPlot(pcaPrcompDataInput)
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
            fig <-
              nonInteractiveUmapPlot(naOmitInput, input$knn)
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
            fig <-
              nonInteractivePcaIndividualsPlot(pcaPrincompDataInput)
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
              nonInteractivePcaVariablesPlot(pcaPrincompDataInput)
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
            column2 <-
              calculateExclusiveColumns(columnNames, group1)
            expect_type(column2, "NULL")
            expect_equal(length(column2), 0)

            # Calculate all$gsms
            all$gsms <-
              calculateEachGroupsSamples(columnNames, group1, group2)
            expect_type(all$gsms, "character")
            expect_equal(all$gsms, "00")
            expect_equal(nchar(all$gsms), 2)

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
            expect_error(
              results <- calculateDifferentialGeneExpression(all$gsms,input,
                                                             all)
            )
          })
