#' A Function to Return the Server Component
#'
#' A Function to Return the Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceServer()
#' @importFrom DT renderDataTable JS
#' @author Guy Hunt
#' @noRd
sourceServer <- function(input, output, session) {
  datasetInformationServer <- ({
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
                           type = "error")
        }

        # Error handling preventing errors
        # caused by non GEO series accession codes
        validate(need(
          substr(input$geoAccessionCode, 1, 3) == "GSE",
          "Please input a GEO series accession code (GSEXXXX)"
        ))

        getGeoObject(input$geoAccessionCode)
      }, error = function(err)
        # Return null if there is a error in the getGeoObject function
        return(NULL))
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
        gsetData <<- tryCatch({
          extractPlatformGset(allGset(), input$platform)
        }, error = function(err)
          # Return null if there is a error in the getGeoObject function
          return(NULL))

        # Error handling to prevent users
        # trying to run exploratory data analysis
        # without selecting a platform
        if (is.null(gsetData) == TRUE) {
          showNotification("Please select a platform.",
                           type = "error")
        } else {
          # Extract expression data
          expressionData <<- extractExpressionData(gsetData)

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
              dataInput <<- tryCatch({
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
                dataInput <<- expressionData
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
              knnDataInput <<- tryCatch({
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
                knnDataInput <<- dataInput
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
              knnColumnInfo <<- knnColumnInfo[knnColumns, ]

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
                      interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
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


      # Error handling to prevent non-microarray
      # datasets being used
      if (is.double(expressionData) == FALSE) {
        showNotification(
          paste0(
            paste0(
              "It appears that the GEO accession code ",
              input$geoAccessionCode
            ),
            " is not a valid microarray gene expression GEO accession code.
            Or the expression data is not in a numerical format."
          ),
          type = "error"
        )
      } else {
        # Differential gene expression analysis
        gsms <-
          calculateEachGroupsSamplesFromDataFrame(as.data.frame(sapply(seq_len(
            nrow(knnColumnInfo)
          ),
          function(a)
            input[[paste0("sel", a)]])))

        # Error handling to ensure at least one
        # group has two samples and the other group
        # has at least one sample
        if ((lengths(regmatches(gsms, gregexpr("0", gsms))) > 0 &
             lengths(regmatches(gsms, gregexpr("1", gsms))) > 1) |
            (lengths(regmatches(gsms, gregexpr("0", gsms))) > 1 &
             lengths(regmatches(gsms, gregexpr("1", gsms))) > 0)) {
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
            tT <- calculateTopDifferentiallyExpressedGenes(fit2,
                                                           adjustment)
            dT <- calculateDifferentialGeneExpressionSummary(fit2,
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
    })

  })
  return(datasetInformationServer)
}
