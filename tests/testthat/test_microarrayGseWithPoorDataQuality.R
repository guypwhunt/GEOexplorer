library(GEOexplorer)
context("With Blank Column")

test_that("Microarray GSE with blank column
          is handled correctly by all functions",
          {
            # Input Values
            logTransformation <- "Auto-Detect"
            knnTransformation <- "Yes"
            knn <- 2
            pValueAdjustment <-
              "Benjamini & Hochberg (False discovery rate)"
            limmaPrecisionWeights <- "Yes"
            forceNormalization <- "Yes"
            platformAnnotation <- "NCBI generated"
            significanceLevelCutOff <- 0.05

            # Get the GEO data for all platforms
            geoAccessionCode <- "GSE18459"
            allGset <- getGeoObject(geoAccessionCode)
            ed <- experimentData(allGset[[1]])
            expect_equal(pubMedIds(ed), "20856801")
            ei <- expinfo(ed)
            expect_equal(ei[1], "Amy,,Pointon", ignore_attr = TRUE)
            expect_equal(ei[2], "", ignore_attr = TRUE) #lab
            expect_equal(ei[3], "", ignore_attr = TRUE)
            expect_equal(
              ei[4],
              "Mouse cardiac tissues treated with doxorubicin and DMNQ",
              ignore_attr = TRUE)
            expect_equal(
              ei[5],
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18459",
              ignore_attr = TRUE
            )

            # Extract platforms
            platforms <- extractPlatforms(allGset)
            platform <- platforms[1]
            expect_type(platforms, 'character')
            expect_type(platform, 'character')
            expect_equal(platform, "GPL9106")

            # Extract the GEO2R data from the specified platform
            gsetData <- extractPlatformGset(allGset, platform)
            expect_type(gsetData, 'S4')
            expect_s4_class(gsetData, 'ExpressionSet')
            expect_equal(nrow(pData(gsetData)), 139)
            expect_equal(nrow(fData(gsetData)), 38468)

            # Extract the experiment information
            experimentInformation <-
              extractExperimentInformation(gsetData)
            expect_type(experimentInformation, 'S4')
            expect_s4_class(experimentInformation, 'MIAME')
            expect_equal(experimentInformation@name, "Amy,,Pointon")
            expect_equal(experimentInformation@lab, "")
            expect_equal(experimentInformation@contact, "")
            expect_equal(
              experimentInformation@title,
              "Mouse cardiac tissues treated with doxorubicin and DMNQ"
            )
            expect_equal(nchar(experimentInformation@title), 55)
            expect_equal(
              experimentInformation@url,
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18459"
            )
            expect_equal(experimentInformation@pubMedIds, "20856801")

            # Extract Sample Information
            sampleInfo <- extractSampleInformation(gsetData)
            expect_type(sampleInfo, 'list')
            expect_equal(nrow(sampleInfo), 139)
            expect_equal(ncol(sampleInfo), 60)

            # Extract expression data
            expressionData <- extractExpressionData(gsetData)
            expect_type(expressionData, 'double')
            expect_equal(ncol(expressionData), 139)
            expect_equal(nrow(expressionData), 38468)

            # Get column Details
            columnInfo <- extractSampleDetails(gsetData)
            expect_type(columnInfo, 'list')
            expect_equal(ncol(columnInfo), 5)
            expect_equal(nrow(columnInfo), 139)

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
            expect_equal(columns[1], "GSM444703")

            # Apply log transformation to expression data if necessary
            dataInput <-
              calculateLogTransformation(expressionData, logTransformation)
            expect_type(dataInput, 'double')
            expect_equal(ncol(dataInput), 139)
            expect_equal(nrow(dataInput), 38468)
            expect_equal(dataInput[1, 1], as.double(NA))

            # Perform KNN transformation on log expression data if necessary
            expect_error(knnDataInput <-
                           calculateKnnImpute(dataInput, "Yes"))
            knnDataInput <- calculateKnnImpute(dataInput, "No")
            expect_type(knnDataInput, 'double')
            expect_equal(ncol(knnDataInput), 139)
            expect_equal(nrow(knnDataInput), 38468)
            expect_equal(knnDataInput[1, 1], as.double(NA))

            # Get a list of all the columns in the KNN output
            knnColumns <- extractSampleNames(knnDataInput)

            # Get knn output column Details
            knnColumnInfo <- extractSampleDetails(gsetData)
            knnColumnInfo <- knnColumnInfo[knnColumns, ]

            # Remove all incomplete rows
            naOmitInput <- calculateNaOmit(knnDataInput)
            expect_type(naOmitInput, 'double')
            expect_equal(ncol(naOmitInput), 139)
            expect_equal(nrow(naOmitInput), 0)

            # Perform Princomp PCA analysis on KNN transformation
            # expression data
            expect_error(pcaPrincompDataInput <-
                           calculatePrincompPca(naOmitInput))

            # Perform Prcomp PCA analysis on KNN transformation expression data
            expect_error(pcaPrcompDataInput <-
                           calculatePrcompPca(naOmitInput))

            # Extract Experiment Information
            extractedExperimentInformation <-
              convertExperimentInformation(experimentInformation)
            expect_type(extractedExperimentInformation, 'character')

            # Non-Interactive Box-and-Whisker Plot
            fig <-
              nonInteractiveBoxAndWhiskerPlot(
                ex = knnDataInput,
                geoAccessionCode = geoAccessionCode,
                platform = platform)
            expect_type(fig, 'list')
            expect_type(fig$stats, 'double')
            expect_type(fig$n, 'double')
            expect_type(fig$conf, 'double')
            expect_type(fig$out, 'double')
            expect_type(fig$group, 'double')
            expect_type(fig$names, 'character')

            # Interactive Box-and-Whisker Plot
            fig <-
              interactiveBoxAndWhiskerPlot(knnDataInput,
                                           geoAccessionCode, platform)
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
            expect_error(
              fig <-
                nonInteractiveDensityPlot(
                  ex = naOmitInput,
                  geoAccessionCode = geoAccessionCode,
                  platform = platform
                )
            )

            # Interactive Density Plot
            expect_error(fig <-
                           interactiveDensityPlot(naOmitInput,
                                                  geoAccessionCode, platform))

            # 3D Interactive Density Plot
            expect_error(fig <-
                           interactiveThreeDDensityPlot(naOmitInput,
                                                        geoAccessionCode,
                                                        platform))

            # Interactive UMAP
            expect_error(fig <-
                           interactiveUmapPlot(naOmitInput, knn,
                                               geoAccessionCode))

            # Interactive Mean Variance Plot
            expect_error(fig <-
                           interactiveMeanVariancePlot(naOmitInput,
                                                       geoAccessionCode,
                                                       gsetData))

            # Correlation Matrix of samples
            expect_error(fig <- interactiveHeatMapPlot(naOmitInput))

            # Non-Interactive UMAP
            expect_error(fig <-
                           nonInteractiveUmapPlot(naOmitInput, knn,
                                                  geoAccessionCode))

            # Non-Interactive Mean Variance Plot
            expect_error(fig <-
                           nonInteractiveMeanVariancePlot(naOmitInput,
                                                          geoAccessionCode))


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
            expect_type(column2, "character")
            expect_equal(column2[1], "GSM444815")
            expect_equal(column2[2], "GSM444816")
            expect_equal(column2[3], "GSM444817")
            expect_equal(length(column2), 69)

            # Calculate gsms
            gsms <-
              calculateEachGroupsSamples(columnNames, group1, group2)
            expect_type(gsms, "character")
            expect_equal(nchar(gsms), 139)

            # Convert P value adjustment
            pValueAdjustment <-
              "Benjamini & Hochberg (False discovery rate)"
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
            expect_error(
              fit2 <-
                calculateDifferentialGeneExpression(
                  gsms,
                  limmaPrecisionWeights,
                  forceNormalization,
                  gsetData,
                  expressionData
                )
            )
          })
