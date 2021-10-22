#' A Function to Return the Server Component
#'
#' A Function to Return the Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceServer()
#' @importFrom DT renderDataTable JS
#' @importFrom shinyBS addTooltip
#' @importFrom utils write.csv
#' @author Guy Hunt
#' @noRd
sourceServer <- function(input, output, session) {
  datasetInformationServer <- ({
    # Add platform tool tip
    addTooltip(
      session,
      id = "platform",
      title = "Each platform relates to a different microarray experiment
      performed in the study.",
      placement = "top",
      trigger = "hover"
    )

    # Add Log Tool Tips
    addTooltip(
      session,
      id = "logTransformation",
      title = "The GEO database accepts a variety of data value types,
              including logged and unlogged data.
              Limma expects data values to be in log space.
              To address this, an auto-detect feature that checks
              the values of selected samples
              and automatically performs a log2 transformation on
              values determined not to be in log space.
              Alternatively, the user can select Yes
              to force log2 transformation,
              or No to override the auto-detect feature.
              The auto-detect feature only considers Sample values that
              have been assigned to a group, and applies the transformation in
              an all-or-none fashion",
      placement = "top",
      trigger = "hover"
    )

    # Add KNN Tool Tips
    addTooltip(
      session,
      id = "knnTransformation",
      title = "Rows with over 50% missing values are imputed using the overall
              mean per sample. Columns with over
              80% will cause an error in the KNN
              computation.",
      placement = "top",
      trigger = "hover"
    )

    # Add P-value adjustment Tool Tips
    addTooltip(
      session,
      id = "pValueAdjustment",
      title =
        "P-value adjustments can be applied to reduce the likelihood of
        a false positive occurring. The P-value adjustment 'None'
        indicates
        no P-value adjustment will be applied and is the least
        conservative
        P-value adjustment. The Benjamini & Hochberg
        (False discovery rate)
        and Benjamini & Yekutieli methods are slightly more
        conservative and aim to control the false discovery rate.
        The Bonferroni
        and Holm methods are the most conservative as they aim to
        control the family-wise error rate.",
      placement = "top",
      trigger = "hover"
    )

    # Add Limma Precision Weights adjustment Tool Tips
    addTooltip(
      session,
      id = "limmaPrecisionWeights",
      title =
        "Limma precision weights should be applied if there is
                  a strong
                  mean-variance trend as can be identified from the
                  'Mean-Variance Plot' tab. By applying limma precision
                  weights the
                  limma vooma function is used to estimate the mean-variance
                  relationship and uses this to compute appropriate
                  observational-level weights.",
      placement = "top",
      trigger = "hover"
    )

    # Add Limma Precision Weights adjustment Tool Tips
    addTooltip(
      session,
      id = "forceNormalization",
      title =
        "Force normalisation should be selected if the gene
                  expression
                  dataset is not normally distributed, as can be identified
                  from the
                  'Box-and-Whisper Plot', the 'Expression Density Plot'
                  and the
                  '3D Expression Density Plot'. By selecting force
                  normalisation
                  quantile normalisation  is applied to the expression dataset
                  making all selected samples have identical value
                distributions.",
      placement = "top",
      trigger = "hover"
    )

    #Add Significance Level Cutoff Tool Tips
    addTooltip(
      session,
      id = "significanceLevelCutOff",
      title =
        "The significance level cut-off is used to identify genes
                  that are
            differentially expressed between the two groups. Genes with
            adjusted P-values less than the significance level cut-off are
            determined to be differentially expressed.",
      placement = "top",
      trigger = "hover"
    )


    # Data Extraction Functions
    # Get the GEO2R data for all platforms
    allGset <- reactive({
      tryCatch({
        # Notify the user the GEO accession code
        # is not a GEO series accession code
        if (substr(input$geoAccessionCode, 1, 3) != "GSE")
        {
          showNotification("Please input a GEO series
                           accession code with the format GSEXXX",
                           type = "warning")
        }

        # Error handling preventing errors
        # caused by non GEO series accession codes
        validate(need(
          substr(input$geoAccessionCode, 1, 3) == "GSE",
          "Please input a GEO series accession code (GSEXXXX)"
        ))

        getGeoObject(input$geoAccessionCode)
      }, error = function(err) {
        # Return null if there is a error in the getGeoObject function
        return(NULL)
      })
    })

    # Error handling to prevent invalid GEO series accession codes being used
    errorCheck <- reactive({
      is.null(allGset())
    })

    # Get a list of all the platforms
    platforms <- reactive({
      extractPlatforms(allGset())
    })

    # Select the top platform
    platform <- reactive({
      platforms()[1]
    })

    # Update Platform Options
    platformObserve <- observe({
      updateSelectInput(session,
                        "platform",
                        choices = platforms(),
                        selected = platform())
    })

    # Exploratory data analysis visualisation
    observeEvent(input$exploratoryDataAnalysisButton, {
      # Clear unused memory
      gc()

      # Set all outputs to blank, this resets
      # all the visualizations to blank after clicking analyse
      output$experimentInfo <- renderUI({

      })
      output$columnTable <- renderDataTable({

      })
      output$table <- renderDataTable({

      })
      output$logTransformationText <- renderUI({

      })
      output$experimentInfo <- renderUI({

      })
      output$knnColumnTable <- renderDataTable({

      })
      output$interactiveBoxAndWhiskerPlot <- renderPlotly({

      })
      output$interactiveDensityPlot <- renderPlotly({

      })
      output$interactiveThreeDDensityPlot <- renderPlotly({

      })
      output$interactiveUmapPlot <- renderPlotly({

      })
      output$interactiveHeatMapPlot <- renderPlotly({

      })
      output$interactiveMeanVariancePlot <- renderPlotly({

      })
      output$interactivePcaScreePlot <- renderPlotly({

      })
      output$interactivePcaIndividualsPlot <- renderPlotly({

      })
      output$interactivePcaVariablesPlot <- renderPlotly({

      })
      output$dETable <- renderDataTable({

      })
      output$iDEHistogram <- renderPlotly({

      })
      output$dEVennDiagram <- renderPlot({

      })
      output$iDEQQ <- renderPlotly({

      })
      output$iDEVolcano <- renderPlotly({

      })
      output$iDEMd <- renderPlotly({

      })

      # Make Differential Gene Expression Action
      # Button Appear, this prevents users
      # trying to perform differential gene expression analysis
      # prior to exploratory data analysis
      output$differentialExpressionButton <- renderUI({
        actionButton("differentialExpressionButton", "Analyse")
      })

      # Error handling to display a notification if an
      # invalid GEO accession code is used.
      if (errorCheck() == TRUE) {
        showNotification(
          paste0(
            paste0("The GEO accession code ",
                   input$geoAccessionCode),
            " is not a valid microarray accession code.
                                Please enter a valid microarray accession code"
          ),
          type = "error"
        )
      } else {
        # Extract the GEO2R data from the specified platform
        gsetData <- tryCatch({
          extractPlatformGset(allGset(), input$platform)
        }, error = function(err) {
          # Return null if there is a error in the getGeoObject function
          return(NULL)
        })

        # Error handling to prevent users
        # trying to run exploratory data analysis
        # without selecting a platform
        if (is.null(gsetData) == TRUE) {
          showNotification("Please select a platform.",
                           type = "error")
        } else {
          # Extract expression data
          expressionData <- extractExpressionData(gsetData)

          # Error handling to prevent issues
          # due to expression data with no samples
          if (length(expressionData) == 0) {
            showNotification(
              "The expression data is empty
            and therefore can not be analysed.
            This may indicate the GEO accession
            code relates to an RNA sequence experiment
                             rather than a microarray experiment.",
              type = "error"
            )
          } else {
            # Extract the experiment information
            experimentInformation <-
              extractExperimentInformation(gsetData)

            # Extract Column Information
            columnInfo <- extractSampleDetails(gsetData)

            # Get a list of all the columns
            columns <- extractSampleNames(expressionData)

            # Error handling to prevent non-microarray GEO
            # accession codes from being used
            if (is.double(expressionData) == FALSE) {
              showNotification(
                paste0(
                  paste0(
                    "It appears that the GEO accession code ",
                    input$geoAccessionCode
                  ),
                  " is not a valid microarray gene
                                    expression GEO accession code.
                                    Please enter a valid microarray gene
                                    expression GEO accession code."
                ),
                type = "error"
              )
              # Experimental Information Display
              output$experimentInfo <- renderUI({
                convertExperimentInformation(experimentInformation)
              })

              # Column Set Plot
              output$columnTable <- renderDataTable({
                columnInfo
              })

              # Expression dataset table
              output$table <- renderDataTable({
                expressionData
              })

            } else {
              # Data Transformation Functions
              # Apply log transformation to expression
              #data if necessary
              dataInput <- tryCatch({
                calculateLogTransformation(expressionData,
                                           input$logTransformation)
              }, error = function(cond) {
                return(NULL)
              })

              # Error handling to display a notification if
              # there was an error in log transformation
              if (is.null(dataInput) == TRUE) {
                showNotification(
                  "There was an error applying log
              transformation to the expression data.
                               Therefore, the original dataset will be used.",
                  type = "warning"
                )
                dataInput <- expressionData
              }

              # Is log transformation auto applied
              autoLogInformation <- tryCatch({
                calculateAutoLogTransformApplication(expressionData)
              }, error = function(cond) {
                return(
                  "There was an error calculating if log transformation
                       would automatically be applied."
                )
              })

              # Perform KNN transformation on log
              # expression data if necessary
              knnDataInput <- tryCatch({
                calculateKnnImpute(dataInput,
                                   input$knnTransformation)
              }, error = function(cond) {
                return(NULL)
              })

              # Error handling to display a notification if
              # there was an error in KNN imputation
              if (is.null(knnDataInput) == TRUE) {
                showNotification(
                  "There was an error applying KNN imputation to the
                               expression data. Therefore, the log-transformed/
                               original dataset
                               will be used instead.",
                  type = "warning"
                )
                knnDataInput <- dataInput
              }
              # Remove all incomplete rows
              naOmitInput <- calculateNaOmit(knnDataInput)

              # Perform PCA analysis on KNN transformation
              # expression data using princomp
              pcaPrcompDataInput <- tryCatch({
                calculatePrcompPca(naOmitInput)
              }, error = function(cond) {
                return(NULL)
              })

              # Data Visualisation Functions
              # Update if log transformation took place
              output$logTransformationText <- renderUI({
                helpText(autoLogInformation)
              })

              # Experimental Information Display
              output$experimentInfo <- renderUI({
                convertExperimentInformation(experimentInformation)
              })

              # Column Set Plot
              output$columnTable <- renderDataTable({
                columnInfo
              })

              # KNN Column Set Plot
              knnColumns <- extractSampleNames(knnDataInput)
              knnColumnInfo <- extractSampleDetails(gsetData)

              # Could turn the below into a function
              knnColumnInfo <- knnColumnInfo[knnColumns, ]

              for (i in seq_len(nrow(knnColumnInfo))) {
                knnColumnInfo$group[i] <- as.character(selectInput(
                  paste0("sel", i),
                  "",
                  choices = unique(c(
                    "N/A", "Group 1", "Group 2"
                  )),
                  width = "100px"
                ))
              }

              output$knnColumnTable <- renderDataTable(
                knnColumnInfo,
                escape = FALSE,
                selection = 'none',
                server = FALSE,
                options = list(
                  dom = 't',
                  paging = FALSE,
                  ordering = FALSE
                ),
                callback =
                  JS(
                    "table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());"
                  )
              )

              # Expression dataset table
              output$table <- renderDataTable({
                knnDataInput
              })

              # Interactive Box-and-Whisker Plot
              output$interactiveBoxAndWhiskerPlot <- renderPlotly({
                interactiveBoxAndWhiskerPlot(knnDataInput,
                                             input$geoAccessionCode,
                                             input$platform)
              })

              # Interactive Density Plot
              output$interactiveDensityPlot <- renderPlotly({
                interactiveDensityPlot(naOmitInput,
                                       input$geoAccessionCode,
                                       input$platform)
              })

              # 3D Interactive Density Plot
              output$interactiveThreeDDensityPlot <- renderPlotly({
                interactiveThreeDDensityPlot(naOmitInput,
                                             input$geoAccessionCode,
                                             input$platform)
              })

              # Error handling to prevent errors caused by
              # expression datasets with only one column
              if (ncol(expressionData) > 1) {
                # Update UMAP KNN max
                updateNumericInput(session, inputId = "knn",
                                   value =2, max = ncol(expressionData))

                # Interactive UMAP Plot
                output$interactiveUmapPlot <- renderPlotly({
                  interactiveUmapPlot(naOmitInput,
                                      input$knn,
                                      input$geoAccessionCode)
                })

                # Heatmap Plot
                output$interactiveHeatMapPlot <- renderPlotly({
                  interactiveHeatMapPlot(naOmitInput)
                })

                # Interactive Mean Variance Plot
                output$interactiveMeanVariancePlot <- renderPlotly({
                  interactiveMeanVariancePlot(naOmitInput,
                                              input$geoAccessionCode,
                                              gsetData)
                })

                # Error handling to display a notification
                # if there was an error in PCA
                if (is.null(pcaPrcompDataInput) == TRUE) {
                  showNotification(
                    "There was an error performing principal component
                                 analysis on the expression data.
                                 Therefore the PCA visualisations
                                 will not be displayed.",
                    type = "warning"
                  )
                } else {
                  # Interactive PCA Scree Plot
                  output$interactivePcaScreePlot <- renderPlotly({
                    interactivePrcompPcaScreePlot(pcaPrcompDataInput,
                                                  input$geoAccessionCode)
                  })

                  # Interactive PCA Individual Plot
                  output$interactivePcaIndividualsPlot <-
                    renderPlotly({
                      interactivePrcompPcaIndividualsPlot(
                        pcaPrcompDataInput,
                        input$geoAccessionCode,
                        gsetData)
                    })

                  # Interactive PCA Variables Plot
                  output$interactivePcaVariablesPlot <-
                    renderPlotly({
                      interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                        input$geoAccessionCode)
                    })
                  showNotification("Exploratory data analysis complete!",
                                   type = "message")
                }
              } else{
                # A notification to the user that only
                # certain data visulisations will be created
                showNotification(
                  "As the expression dataset had only one
                               column only the Box-and-Whisper Plot
                               and Expression Density
                               Plots will be produced.",
                  type = "warning"
                )
              }
            }
          }
        }
      }
    })

    # Differential Gene Expression Functions
    observeEvent(input$differentialExpressionButton, {
      # Clear unused memory
      gc()

      # Set all differential gene expression
      # analysis outputs to blank, this resets
      # all the visualizations to blank after
      # clicking analyse
      output$dETable <- renderDataTable({

      })
      output$iDEHistogram <- renderPlotly({

      })
      output$dEVennDiagram <- renderPlot({

      })
      output$iDEQQ <- renderPlotly({

      })
      output$iDEVolcano <- renderPlotly({

      })
      output$iDEMd <- renderPlotly({

      })

      # Error handling to display a notification if an
      # invalid GEO accession code is used.
      if (errorCheck() == TRUE) {
        showNotification(
          paste0(
            paste0("The GEO accession code ",
                   input$geoAccessionCode),
            " is not a valid microarray accession code.
                                Please enter a valid microarray accession code"
          ),
          type = "error"
        )
      } else {
        # Extract the GEO2R data from the specified platform
        gsetData <- tryCatch({
          extractPlatformGset(allGset(), input$platform)
        }, error = function(err) {
          # Return null if there is a error in the getGeoObject function
          return(NULL)
        })

        # Error handling to prevent users
        # trying to run exploratory data analysis
        # without selecting a platform
        if (is.null(gsetData) == TRUE) {
          showNotification("Please select a platform.",
                           type = "error")
        } else {
          # Extract expression data
          expressionData <- extractExpressionData(gsetData)

          # Error handling to prevent issues
          # due to expression data with no samples
          if (length(expressionData) == 0) {
            showNotification(
              "The expression data is empty
            and therefore can not be analysed.
            This may indicate the GEO accession
            code relates to an RNA sequence experiment
                             rather than a microarray experiment.",
              type = "error"
            )
          } else {
            # Error handling to prevent non-microarray GEO
            # accession codes from being used
            if (is.double(expressionData) == FALSE) {
              showNotification(
                paste0(
                  paste0(
                    "It appears that the GEO accession code ",
                    input$geoAccessionCode
                  ),
                  " is not a valid microarray gene
                                    expression GEO accession code.
                                    Please enter a valid microarray gene
                                    expression GEO accession code."
                ),
                type = "error"
              )
            } else {
              # Data Transformation Functions
              # Apply log transformation to expression
              #data if necessary
              dataInput <- tryCatch({
                calculateLogTransformation(expressionData,
                                           input$logTransformation)
              }, error = function(cond) {
                return(NULL)
              })

              # Error handling to display a notification if
              # there was an error in log transformation
              if (is.null(dataInput) == TRUE) {
                dataInput <- expressionData
              }

              # Perform KNN transformation on log
              # expression data if necessary
              knnDataInput <- tryCatch({
                calculateKnnImpute(dataInput,
                                   input$knnTransformation)
              }, error = function(cond) {
                return(NULL)
              })

              # Error handling to display a notification if
              # there was an error in KNN imputation
              if (is.null(knnDataInput) == TRUE) {
                knnDataInput <- dataInput
              }

              # KNN Column Set Plot
              knnColumns <- extractSampleNames(knnDataInput)
              knnColumnInfo <- extractSampleDetails(gsetData)

              # Could turn the below into a function
              knnColumnInfo <- knnColumnInfo[knnColumns,]

              # Differential gene expression analysis
              gsms <- tryCatch({
                calculateEachGroupsSamplesFromDataFrame(
                  as.data.frame(
                    sapply(
                      seq_len(
                        nrow(knnColumnInfo)
                      ),
                      function(a)
                        input[[paste0("sel", a)]])))

              }, error = function(cond) {
                return(NULL)
              })

              # Error handling to prevent differential gene expression
              # analysis being performed before exploratory data analysis
              if (is.null(gsms)) {
                showNotification(
                  "There was an error running differential gene expression
                  analysis. Please ensure you have performed exploratory data
                  analysis first.",
                  type = "error"
                )
              } else {
                # Error handling to ensure at least one
                # group has two samples and the other group
                # has at least one sample
                if ((lengths(regmatches(
                  gsms, gregexpr("0", gsms)
                )) > 0 &
                lengths(regmatches(
                  gsms, gregexpr("1", gsms)
                )) > 1) |
                (lengths(regmatches(
                  gsms, gregexpr("0", gsms)
                )) > 1 &
                lengths(regmatches(
                  gsms, gregexpr("1", gsms)
                )) > 0)) {
                  fit2 <- tryCatch({
                    calculateDifferentialGeneExpression(
                      gsms,
                      input$limmaPrecisionWeights,
                      input$forceNormalization,
                      gsetData,
                      knnDataInput
                    )
                  }
                  , error = function(cond) {
                    return(NULL)
                  })

                  # Error handling to ensure Differential Gene
                  # Expression Analysis worked
                  if (is.null(fit2)) {
                    showNotification(
                      "There was an error calculating the
                             differential gene expression analysis!",
                      type = "error"
                    )
                  } else {
                    adjustment <- convertAdjustment(input$pValueAdjustment)
                    tT <-
                      calculateTopDifferentiallyExpressedGenes(fit2,
                                                               adjustment)
                    dT <-
                      calculateDifferentialGeneExpressionSummary(
                        fit2,
                        adjustment,
                        input$significanceLevelCutOff)

                    ct <- 1

                    # Differential gene expression table
                    output$dETable <- renderDataTable({
                      as.data.frame(tT)
                    })

                    # Interactive Histogram Plot
                    output$iDEHistogram <- renderPlotly({
                      interactiveHistogramPlot(fit2, adjustment)
                    })

                    # Venn Diagram Plot
                    output$dEVennDiagram <- renderPlot({
                      nonInteractiveVennDiagramPlot(dT)
                    })

                    # Interactive QQ Plot
                    output$iDEQQ <- renderPlotly({
                      interactiveQQPlot(fit2, dT, ct)
                    })

                    # Interactive Volcano Plot
                    output$iDEVolcano <- renderPlotly({
                      interactiveVolcanoPlot(fit2, dT, ct)
                    })

                    # Interactive Mean Difference Plot
                    output$iDEMd <- renderPlotly({
                      interactiveMeanDifferencePlot(fit2, dT, ct)
                    })

                    # Download Top Differentially Expressed Genes Table
                    output$downloadData <- downloadHandler(
                      filename = function() {
                        "top_differentially_expressed_genes.csv"
                      },
                      content = function(file) {
                        write.csv(tT, file, row.names = FALSE)
                      }
                    )
                    showNotification("Differential gene
                           expression analysis complete!",
                                     type = "message")
                  }
                } else {
                  showNotification(
                    "One group needs at
          least 2 samples and the other
                           group needs at least 1 sample",
                    type = "error"
                  )
                }
              }
            }
          }
        }
      }
    })
  })
  return(datasetInformationServer)
}
