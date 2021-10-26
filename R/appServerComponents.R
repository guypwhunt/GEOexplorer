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
    # Define variables
    all <- reactiveValues()
    ct <- 1

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

    # Add Force normalisation adjustment Tool Tips
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

    observeEvent(input$dataSource, {
      # Update UI side bar with GEO widgets
      if (input$dataSource == "GEO") {
        output$output1 <- renderUI({
          helpText(
            "Input a GEO series accession code (GSEXXXX format)
      to examine the gene expression data.
      This can be obtained from https://www.ncbi.nlm.nih.gov/gds."
          )
        })
        output$output2 <- renderUI({
          textInput("geoAccessionCode", "GEO accession code", "")
        })
        output$output3 <- renderUI({
          helpText("Select the platform of interest.")
        })
        output$output4 <- renderUI({
          selectInput("platform", "Platform", c())
        })

        # Add platform tool tip
        addTooltip(
          session,
          id = "platform",
          title = "Each platform relates to a different microarray experiment
      performed in the study.",
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
        observeEvent(errorCheck(), {
          if (errorCheck() == FALSE) {
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
          }
        })

        # Exploratory data analysis visualisation
        observeEvent(input$exploratoryDataAnalysisButton, {
          # Clear unused memory
          gc()

          # Set all outputs to blank, this resets
          # all the visualizations to blank after clicking analyse
          output$experimentInfo <- renderUI({

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

          output$interactivePcaIndividualsPlot <- renderPlotly({

          })
          output$interactivePcaVariablesPlot <- renderPlotly({

          })
          output$interactive3DPcaVariablesPlot <- renderPlotly({

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

          # Update Experimental informatin
          output$experimentInfo <-  renderUI({"Experimental Information is
          not available when processing user-uploaded files!"})

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
            all$gsetData <- tryCatch({
              extractPlatformGset(allGset(), input$platform)
            }, error = function(err) {
              # Return null if there is a error in the getGeoObject function
              return(NULL)
            })

            # Error handling to prevent users
            # trying to run exploratory data analysis
            # without selecting a platform
            if (is.null(all$gsetData) == TRUE) {
              showNotification("Please select a platform.",
                               type = "error")
            } else {
              # Extract expression data
              all$expressionData <- extractExpressionData(all$gsetData)

              # Error handling to prevent issues
              # due to expression data with no samples
              if (length(all$expressionData) == 0) {
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
                  extractExperimentInformation(all$gsetData)

                # Extract Column Information
                columnInfo <- extractSampleDetails(all$gsetData)

                # Get a list of all the columns
                columns <- extractSampleNames(all$expressionData)

                # Error handling to prevent non-microarray GEO
                # accession codes from being used
                if (is.double(all$expressionData) == FALSE) {
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
                  output$experimentInfo <- tryCatch({
                    renderUI({
                      convertExperimentInformation(experimentInformation)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })


                  # Column Set Plot
                  output$columnTable <-
                    tryCatch({
                      renderDataTable({
                        columnInfo
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })

                  # Expression dataset table
                  output$table <-
                    tryCatch({
                      renderDataTable({
                        all$expressionData
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                } else {
                  # Data Transformation Functions
                  # Apply log transformation to expression
                  #data if necessary
                  all$dataInput <- tryCatch({
                    calculateLogTransformation(all$expressionData,
                                               input$logTransformation)
                  }, error = function(cond) {
                    return(NULL)
                  })

                  # Error handling to display a notification if
                  # there was an error in log transformation
                  if (is.null(all$dataInput) == TRUE) {
                    showNotification(
                      "There was an error applying log
              transformation to the expression data.
                               Therefore, the original dataset will be used.",
                      type = "warning"
                    )
                    all$dataInput <- all$expressionData
                  }

                  # Is log transformation auto applied
                  autoLogInformation <- tryCatch({
                    calculateAutoLogTransformApplication(all$expressionData)
                  }, error = function(cond) {
                    return(
                      "There was an error calculating if log transformation
                       would automatically be applied."
                    )
                  })

                  # Perform KNN transformation on log
                  # expression data if necessary
                  all$knnDataInput <- tryCatch({
                    calculateKnnImpute(all$dataInput,
                                       input$knnTransformation)
                  }, error = function(cond) {
                    return(NULL)
                  })

                  # Error handling to display a notification if
                  # there was an error in KNN imputation
                  if (is.null(all$knnDataInput) == TRUE) {
                    showNotification(
                      "There was an error applying KNN imputation to the
                               expression data. Therefore, the log-transformed/
                               original dataset
                               will be used instead.",
                      type = "warning"
                    )
                    all$knnDataInput <- all$dataInput
                  }
                  # Remove all incomplete rows
                  naOmitInput <- calculateNaOmit(all$knnDataInput)

                  # Perform PCA analysis on KNN transformation
                  # expression data using princomp
                  pcaPrcompDataInput <- tryCatch({
                    calculatePrcompPca(naOmitInput)
                  }, error = function(cond) {
                    return(NULL)
                  })

                  # Data Visualisation Functions
                  # Update if log transformation took place
                  output$logTransformationText <-
                    tryCatch({
                      renderUI({
                        helpText(autoLogInformation)
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                  # Experimental Information Display
                  output$experimentInfo <-
                    tryCatch({
                      renderUI({
                        convertExperimentInformation(experimentInformation)
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })

                  # Column Set Plot
                  output$columnTable <-
                    tryCatch({
                      renderDataTable({
                        columnInfo
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })

                  # KNN Column Set Plot
                  all$knnColumns <- extractSampleNames(all$knnDataInput)
                  all$knnColumnInfo <-
                    extractSampleDetails(all$gsetData)

                  # Could turn the below into a function
                  all$knnColumnInfo <-
                    all$knnColumnInfo[all$knnColumns,]

                  for (i in seq_len(nrow(all$knnColumnInfo))) {
                    all$knnColumnInfo$group[i] <- as.character(selectInput(
                      paste0("sel", i),
                      "",
                      choices = unique(c(
                        "N/A", "Group 1", "Group 2"
                      )),
                      width = "100px"
                    ))
                  }

                  output$knnColumnTable <- tryCatch({
                    renderDataTable(
                      all$knnColumnInfo,
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
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })

                  # Expression dataset table
                  output$table <-
                    tryCatch({
                      renderDataTable({
                        all$knnDataInput
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                  # Interactive Box-and-Whisker Plot
                  output$interactiveBoxAndWhiskerPlot <-
                    tryCatch({
                      renderPlotly({
                        interactiveBoxAndWhiskerPlot(all$knnDataInput,
                                                     input$geoAccessionCode,
                                                     input$platform)
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                  # Interactive Density Plot
                  output$interactiveDensityPlot <-
                    tryCatch({
                      renderPlotly({
                        interactiveDensityPlot(naOmitInput,
                                               input$geoAccessionCode,
                                               input$platform)
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                  # 3D Interactive Density Plot
                  output$interactiveThreeDDensityPlot <-
                    tryCatch({
                      renderPlotly({
                        interactiveThreeDDensityPlot(naOmitInput,
                                                     input$geoAccessionCode,
                                                     input$platform)
                      })
                    },
                    error = function(e) {
                      # return a safeError if a parsing error occurs
                      stop(safeError(e))
                    })


                  # Error handling to prevent errors caused by
                  # expression datasets with only one column
                  if (ncol(all$expressionData) > 1) {
                    # Update UMAP KNN max
                    updateNumericInput(
                      session,
                      inputId = "knn",
                      value = 2,
                      max = ncol(all$expressionData)
                    )

                    # Interactive UMAP Plot
                    output$interactiveUmapPlot <-
                      tryCatch({
                        renderPlotly({
                          interactiveUmapPlot(naOmitInput,
                                              input$knn,
                                              input$geoAccessionCode)
                        })
                      },
                      error = function(e) {
                        # return a safeError if a parsing error occurs
                        stop(safeError(e))
                      })


                    # Heatmap Plot
                    output$interactiveHeatMapPlot <-
                      tryCatch({
                        renderPlotly({
                          interactiveHeatMapPlot(naOmitInput)
                        })
                      },
                      error = function(e) {
                        # return a safeError if a parsing error occurs
                        stop(safeError(e))
                      })

                    # Interactive Mean Variance Plot
                    output$interactiveMeanVariancePlot <-
                      tryCatch({
                        renderPlotly({
                          interactiveMeanVariancePlot(naOmitInput,
                                                      input$geoAccessionCode,
                                                      all$gsetData)
                        })
                      },
                      error = function(e) {
                        # return a safeError if a parsing error occurs
                        stop(safeError(e))
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
                      # Interactive PCA Individual Plot
                      output$interactivePcaIndividualsPlot <-
                        tryCatch({
                          renderPlotly({
                            interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
                                                                input$geoAccessionCode,
                                                                all$gsetData)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })


                      # Interactive PCA Variables Plot
                      output$interactivePcaVariablesPlot <-
                        tryCatch({
                          renderPlotly({
                            interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                              input$geoAccessionCode)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })

                      # Interactive 3D PCA Variables Plot
                      output$interactive3DPcaVariablesPlot <-
                        tryCatch({
                          renderPlotly({
                            interactive3DPrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                              input$geoAccessionCode)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
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
            # Error handling to prevent users
            # trying to run exploratory data analysis
            # without selecting a platform
            if (is.null(all$gsetData) == TRUE) {
              showNotification("Please select a platform.",
                               type = "error")
            } else {
              # Error handling to prevent issues
              # due to expression data with no samples
              if (length(all$expressionData) == 0) {
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
                if (is.double(all$expressionData) == FALSE) {
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
                  # Error handling to display a notification if
                  # there was an error in log transformation
                  if (is.null(all$dataInput) == TRUE) {
                    all$dataInput <- all$expressionData
                  }

                  # Error handling to display a notification if
                  # there was an error in KNN imputation
                  if (is.null(all$knnDataInput) == TRUE) {
                    all$knnDataInput <- all$dataInput
                  }

                  # Differential gene expression analysis
                  gsms <- tryCatch({
                    calculateEachGroupsSamplesFromDataFrame(as.data.frame(sapply(seq_len(nrow(all$knnColumnInfo)),
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
                          all$gsetData,
                          all$knnDataInput
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
                          calculateDifferentialGeneExpressionSummary(fit2,
                                                                     adjustment,
                                                                     input$significanceLevelCutOff)
                        # Differential gene expression table
                        output$dETable <- tryCatch({
                          renderDataTable({
                            as.data.frame(tT)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })

                        # Interactive Histogram Plot
                        output$iDEHistogram <- tryCatch({
                          renderPlotly({
                            interactiveHistogramPlot(fit2, adjustment)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })


                        # Venn Diagram Plot
                        output$dEVennDiagram <- tryCatch({
                          renderPlot({
                            nonInteractiveVennDiagramPlot(dT)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })


                        # Interactive QQ Plot
                        output$iDEQQ <-
                          tryCatch({
                            renderPlotly({
                              interactiveQQPlot(fit2, dT, ct)
                            })
                          },
                          error = function(e) {
                            # return a safeError if a parsing error occurs
                            stop(safeError(e))
                          })


                        # Interactive Volcano Plot
                        output$iDEVolcano <- tryCatch({
                          renderPlotly({
                            interactiveVolcanoPlot(fit2, dT, ct)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
                        })


                        # Interactive Mean Difference Plot
                        output$iDEMd <- tryCatch({
                          renderPlotly({
                            interactiveMeanDifferencePlot(fit2, dT, ct)
                          })
                        },
                        error = function(e) {
                          # return a safeError if a parsing error occurs
                          stop(safeError(e))
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

      } else if
      (input$dataSource == "Upload") {
        # Update UI side bar with User Upload widgets
        output$output1 <- renderUI({
          fileInput("file2", "Upload CSV Experimental Conditions File",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv"))
        })
        output$output2 <- renderUI({
          fileInput("file1", "Upload CSV Gene Expression Count File",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv"))
        })
        output$output <- renderUI({
          radioButtons(
            "typeOfData",
            label = "Microarray or RNA Sequencing Data?",
            choices = list("Microarray", "RNA Sequencing"),
            selected = "Microarray"
          )
        })
        output$output4 <- renderUI({
          radioButtons(
            "cpmTransformation",
            label = "Convert to count per million:",
            choices = list("Yes", "No"),
            selected = "No"
          )
        })
        # Add CPM tool tip
        bsTooltip(
          id = "cpmTransformation",
          title = "This is recommended for raw RNA sequence data.",
          placement = "top",
          trigger = "hover"
        )

        # Define Variables
        gsetData <- NULL
        geoAccessionCode <- ""

        # Extract Expression Data
        expressionData <- reactive({
          # Ensure a file has been uploaded
          req(input$file1)

          # Extract CSV
          tryCatch({
            expressionDataDf <- readCsvFile(input$file1$datapath)
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })

          try({
            # Preprocess the data
            expressionDataDf <-
              preProcessGeneExpressionData(expressionDataDf)
          })

          return(expressionDataDf)

        })


        # Expression dataset table
        output$table <-
          tryCatch({
            renderDataTable(expressionData())
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })


        # Extract Column Information
        columnInfo <- reactive({
          # Ensure a file has been uploaded
          req(input$file2)

          # Try reading the CSV
          tryCatch({
            columnInfoDf <- readCsvFile(input$file2$datapath)
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })

          return(columnInfoDf)
        })


        # Experimental conditions table
        output$columnTable <-
          tryCatch({
            renderDataTable({
              columnInfo()
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })


        # Exploratory Data Analysis steps
        observeEvent(input$exploratoryDataAnalysisButton, {
          # Clear unused memory
          gc()

          # Set all outputs to blank, this resets
          # all the visualizations to blank after clicking analyse
          output$logTransformationText <- renderUI({

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
          output$interactivePcaIndividualsPlot <- renderPlotly({

          })
          output$interactivePcaVariablesPlot <- renderPlotly({

          })
          output$interactive3DPcaVariablesPlot <- renderPlotly({

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

          # Get a list of all the columns
          columns <- extractSampleNames(expressionData())

          if (input$typeOfData == "RNA Sequencing")
          {
            # Raw counts are converted to counts-per-million (CPM)
            all$cpm <-
              calculateCountsPerMillion(expressionData(), input$cpmTransformation)
          } else if (input$typeOfData == "Microarray")
          {
            all$cpm <- expressionData()
          }

          autoLogInformation <-
            calculateAutoLogTransformApplication(all$cpm)

          # Apply log transformation to expression data if necessary
          all$dataInput <-
            calculateLogTransformation(all$cpm, input$logTransformation)

          # Perform KNN transformation on log expression data if necessary
          all$knnDataInput <-
            calculateKnnImpute(all$dataInput, input$knnTransformation)

          # Get a list of all the columns in the KNN output
          knnColumns <- extractSampleNames(all$knnDataInput)

          # Get knn output column Details
          all$knnColumnInfo <- columnInfo()
          row.names(all$knnColumnInfo) <- all$knnColumnInfo$column
          all$knnColumnInfo <- all$knnColumnInfo[knnColumns,]

          # Remove all incomplete rows
          naOmitInput <- calculateNaOmit(all$knnDataInput)

          # Perform Prcomp PCA analysis on KNN transformation expression data
          pcaPrcompDataInput <- calculatePrcompPca(naOmitInput)

          # Data Visualisation Functions
          # Update if log transformation took place
          output$logTransformationText <-
            tryCatch({
              renderUI({
                helpText(autoLogInformation)
              })
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })


          # Generate Differential Gene Expression Table
          for (i in seq_len(nrow(all$knnColumnInfo))) {
            all$knnColumnInfo$group[i] <- as.character(selectInput(
              paste0("sel", i),
              "",
              choices = unique(c("N/A", "Group 1", "Group 2")),
              width = "100px"
            ))
          }

          output$knnColumnTable <-
            tryCatch({
              renderDataTable(
                all$knnColumnInfo,
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
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })



          # Expression dataset table
          output$table <- tryCatch({
            renderDataTable({
              all$knnDataInput
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })

          # Interactive Box-and-Whisker Plot
          output$interactiveBoxAndWhiskerPlot <-
            tryCatch({
              renderPlotly({
                interactiveBoxAndWhiskerPlot(all$knnDataInput,
                                             geoAccessionCode,
                                             input$platform)
              })
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })



          # Interactive Density Plot
          output$interactiveDensityPlot <-
            tryCatch({
              renderPlotly({
                interactiveDensityPlot(naOmitInput,
                                       geoAccessionCode,
                                       input$platform)
              })
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })


          # 3D Interactive Density Plot
          output$interactiveThreeDDensityPlot <-
            tryCatch({
              renderPlotly({
                interactiveThreeDDensityPlot(naOmitInput,
                                             geoAccessionCode,
                                             input$platform)
              })
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })



          # Error handling to prevent errors caused by
          # expression datasets with only one column
          if (ncol(naOmitInput) > 1) {
            # Update UMAP KNN max
            updateNumericInput(
              session,
              inputId = "knn",
              value = 2,
              max = ncol(naOmitInput)
            )

            # Interactive UMAP Plot
            output$interactiveUmapPlot <-
              tryCatch({
                renderPlotly({
                  interactiveUmapPlot(naOmitInput,
                                      input$knn,
                                      geoAccessionCode)
                })
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
              })




            # Heatmap Plot
            output$interactiveHeatMapPlot <-
              tryCatch({
                renderPlotly({
                  interactiveHeatMapPlot(naOmitInput)
                })
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
              })


            # Interactive Mean Variance Plot
            output$interactiveMeanVariancePlot <-
              tryCatch({
                renderPlotly({
                  interactiveMeanVariancePlot(naOmitInput,
                                              geoAccessionCode,
                                              gsetData)
                })
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
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
              # Interactive PCA Individual Plot
              output$interactivePcaIndividualsPlot <-
                tryCatch({
                  renderPlotly({
                    interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
                                                        geoAccessionCode,
                                                        gsetData)
                  })
                },
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError(e))
                })

              # Interactive PCA Variables Plot
              output$interactivePcaVariablesPlot <-
                tryCatch({
                  renderPlotly({
                    interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                      geoAccessionCode)
                  })
                },
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError(e))
                })

              # Interactive 3D PCA Variables Plot
              output$interactive3DPcaVariablesPlot <-
                tryCatch({
                  renderPlotly({
                    interactive3DPrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                      geoAccessionCode)
                  })
                },
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError(e))
                })

              # Show notification that Exploratory data analysis finished
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
        })

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

          # Differential gene expression analysis
          gsms <- tryCatch({
            calculateEachGroupsSamplesFromDataFrame(as.data.frame(sapply(seq_len(
              nrow(all$knnColumnInfo)
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
            if ((lengths(regmatches(gsms, gregexpr("0", gsms))) > 0 &
                 lengths(regmatches(gsms, gregexpr("1", gsms))) > 1) |
                (lengths(regmatches(gsms, gregexpr("0", gsms))) > 1 &
                 lengths(regmatches(gsms, gregexpr("1", gsms))) > 0)) {
              if (input$typeOfData == "RNA Sequencing")
              {
                fit2 <-
                  tryCatch({
                    calculateDifferentialGeneExpressionRnaSeq(
                      all$knnDataInput,
                      gsms,
                      input$limmaPrecisionWeights,
                      input$forceNormalization
                    )
                  }
                  , error = function(cond) {
                    return(NULL)
                  })
                if (is.null(fit2)) {
                  # Try again with non-log data
                  knnDataInput <-
                    calculateKnnImpute(all$cpm, input$knnTransformation)
                  fit22 <- tryCatch({
                    calculateDifferentialGeneExpressionRnaSeq(
                      knnDataInput,
                      gsms,
                      input$limmaPrecisionWeights,
                      input$forceNormalization
                    )
                  }
                  , error = function(cond) {
                    return(NULL)
                  })

                  if (is.null(fit22)) {
                    showNotification(
                      "There was an error calculating the
                             differential gene expression analysis!",
                      type = "error"
                    )
                  } else
                  {
                    # Update fit2
                    fit2 <- fit22

                    # Show warning that non-log data was used
                    showNotification(
                      "There was an error calculating the
                             differential gene expression analysis
              using the log data. So the non-log data was used instead!",
                      type = "warning"
                    )
                  }
                }
              } else if (input$typeOfData == "Microarray") {
                fit2 <- tryCatch({
                  calculateDifferentialGeneExpressionMicroarray(
                    all$knnDataInput,
                    gsms,
                    input$limmaPrecisionWeights,
                    input$forceNormalization
                  )
                }
                , error = function(cond) {
                  return(NULL)
                })
              }
              if (is.null(fit2) == FALSE) {
                # Convert the UI adjustment into the value needed for the backend
                adjustment <- convertAdjustment(input$pValueAdjustment)

                # Calculate the top differentially expressed genes
                tT <-
                  calculateTopDifferentiallyExpressedGenes(fit2,
                                                           adjustment)

                # Calculate genes that are upregulated and downregulated
                dT <-
                  calculateDifferentialGeneExpressionSummary(fit2,
                                                             adjustment,
                                                             input$significanceLevelCutOff)

                # Differential gene expression table
                output$dETable <-
                  tryCatch({
                    renderDataTable({
                      as.data.frame(tT)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })

                # Interactive Histogram Plot
                output$iDEHistogram <-
                  tryCatch({
                    renderPlotly({
                      interactiveHistogramPlot(fit2, adjustment)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })


                # Venn Diagram Plot
                output$dEVennDiagram <-
                  tryCatch({
                    renderPlot({
                      nonInteractiveVennDiagramPlot(dT)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })



                # Interactive QQ Plot
                output$iDEQQ <-
                  tryCatch({
                    renderPlotly({
                      interactiveQQPlot(fit2, dT, ct)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })



                # Interactive Volcano Plot
                output$iDEVolcano <-
                  tryCatch({
                    renderPlotly({
                      interactiveVolcanoPlot(fit2, dT, ct)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
                  })

                # Interactive Mean Difference Plot
                output$iDEMd <-
                  tryCatch({
                    renderPlotly({
                      interactiveMeanDifferencePlot(fit2, dT, ct)
                    })
                  },
                  error = function(e) {
                    # return a safeError if a parsing error occurs
                    stop(safeError(e))
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
              } else {
                showNotification("There was an error in differential gene expression
              analysis.",
                                 type = "error")
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
      }
    })

  })
  return(datasetInformationServer)
}
