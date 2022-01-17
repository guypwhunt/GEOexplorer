#' A Function to Return the Server Component
#'
#' A Function to Return the Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceServer()
#' @importFrom DT renderDataTable
#' @importFrom utils write.csv object.size
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @importFrom stringr str_trim
#' @import markdown
#' @importFrom knitr knit
#' @author Guy Hunt
#' @noRd
sourceServer <- function(input, output, session) {
  serverComponents <- ({
    # Update the max number of expressions
    try(options(expressions = 500000))

    # Update the timeout
    try(options(timeout = 300))

    # Update default max file upload
    try(options(shiny.maxRequestSize = 15*1024^2))

    # Common steps
    # Define variables
    all <- reactiveValues()
    errorChecks <- reactiveValues()
    ct <- 1
    userUploadExperimentInformation <- HTML(
      "<p>Experimental
          Information is not available when processing
          user-uploaded files!</p><br>"
    )

    # Update enrichment database list
    databaseNames <- tryCatch({extractDatabaseNamesFromEnrichR()},
                              error = function(e) {
                                # return a safeError if a parsing error occurs
                                return(NULL)
                              })
    if (!is.null(databaseNames)) {
      updateSelectInput(session, "enrichDatabases", choices = databaseNames,
                        selected = "GO_Biological_Process_2021")
    }

    # Define error checks
    errorChecks <- resetErrorChecks(errorChecks)

    # Search GEO server actions
    observeEvent(input$searchGeo, {
      # Define Variables
      firstResultNumber <- "0"
      lastResultNumber <- "50"

      # Update Variables
      all$geoSearchTerm <- input$geoSearchTerm
      all$geoSearchResults <- tryCatch({
        searchGeo(input$geoSearchTerm, firstResultNumber, lastResultNumber)},
        error = function(e) {
          # return a safeError if a parsing error occurs
          return(NULL)
        })


      # Load search results table
      output$geoSearchResults <- tryCatch({
        renderDataTable(
          all$geoSearchResults,
          server = FALSE,
          escape = FALSE,
          selection = 'none'
        )
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })
    })

    # Load GEO Accession Code from GEO Search
    observeEvent(input$loadGeoSearchAsFirstDataset, {
      selectedRow <- as.numeric(
        strsplit(input$loadGeoSearchAsFirstDataset, "_")[[1]][2])

      # Change to Home Tab
      updateTabsetPanel(session, "geoexplorerNavBar",
                        selected = "Home")

      # Update the Radio button to enable the dataset to be processed
      updateRadioButtons(session, inputId = "dataSource", selected = "GEO")

      # Update UI
      loadDataSetUiComponents(input, output, session, errorChecks, all,
                              userUploadExperimentInformation)

      # Add GEO accession input
      output$output5 <- renderUI({
        textInput("geoAccessionCode",
                  "GEO accession code",
                  all$geoSearchResults[selectedRow, 1])
      })
    })

    # Load the example dataset, configurations and perform exploratory data
    # analysis and differential gene expression analysis
    observeEvent(input$loadExampleData, {
      # Update the two Radio buttons to enable the dataset to be processed
      updateRadioButtons(session, inputId = "dataSetType", selected = "Single")
      updateRadioButtons(session, inputId = "dataSource", selected = "GEO")

      # Change to Home Tab
      updateTabsetPanel(session, "geoexplorerNavBar",
                        selected = "Home")

      # Add GEO accession input
      output$output5 <- renderUI({
        textInput("geoAccessionCode", "GEO accession code", "GSE18388")
      })

      # Trigger load GEO dataset when the GEO accession code updated
      observeEvent(input$geoAccessionCode, {
        # Load GEO dataset
        try(loadGeoDataset(input, output, session, errorChecks, all))

        # Perform Exploratory data analysis
        try(performExploratoryDataAnalysis(input,
                                       output,
                                       session,
                                       errorChecks,
                                       all,
                                       userUploadExperimentInformation))

        # Perform Differential Gene Expression Analysis
        try(performDifferentialGeneExpressionAnalysis(input, output, session,
                                                  errorChecks, all, ct,
                                                  exampleDataSet = TRUE))

        try(performGeneEnrichmentAnalysis(input,
                                          output,
                                          session,
                                          errorChecks,
                                          all,
                                          databaseNames,
                                          4))
        })
    })

    # Download gene expression template
    geneExpressionTemplate <- tryCatch({
      as.matrix(geneExpressionTemplate)
    }, error = function(e) {
      return(NULL)
    })
    output$downloadGeneExpressionFileTemplate <-
      try(
        dowloadFile("gene_expression_template.csv", geneExpressionTemplate))

    # Download microarray example dataset
    microarrayExampleDataset <- tryCatch({
      as.matrix(microarrayExampleGeneExpressionCsv)
    }, error = function(e) {
      return(NULL)
    })
    output$downloadMicroarrayExample <-
      try(
        dowloadFile(
          "microarray_example_dataset.csv", microarrayExampleDataset))

    # Download RNAseq example dataset
    rnaSeqExampleDataset <-
      tryCatch({
        as.matrix(rnaSeqExampleGeneExpressionCsv)
      }, error = function(e) {
        return(NULL)
      })
    output$downloadRnaSeqExample <-
      try(
        dowloadFile("rna_seq_example_dataset.csv", rnaSeqExampleDataset)
    )

    # Load logic to update UI
    observeEvent(input$dataSource, loadDataSetUiComponents(
      input, output, session, errorChecks, all,
      userUploadExperimentInformation))
    observeEvent(input$dataSetType,
                 loadDataSetCombinationUiComponents(input, output,
                                                    session, errorChecks, all))

    # Exploratory data analysis visualisation
    observeEvent(input$exploratoryDataAnalysisButton,
                 performExploratoryDataAnalysis(input, output, session,
                                                errorChecks, all,
                                                userUploadExperimentInformation
                                                ))

    # Differential Gene Expression Functions
    observeEvent(input$differentialExpressionButton,
                 performDifferentialGeneExpressionAnalysis(input,
                                                            output,
                                                            session,
                                                            errorChecks,
                                                            all,
                                                            ct))

    # Gene Enrichment Functions
    observeEvent(input$enrichmentAnalysisButton,
                 performGeneEnrichmentAnalysis(input,
                                               output,
                                               session,
                                               errorChecks,
                                               all,
                                               databaseNames))

  })
  return(serverComponents)
}



#' A Function to Return the Exploratory Data Analysis Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT renderDataTable
#' @importFrom utils write.csv object.size
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @importFrom stringr str_trim
#' @import markdown
#' @importFrom knitr knit
#' @author Guy Hunt
#' @noRd
performExploratoryDataAnalysis <- function(input,
                                           output,
                                           session,
                                           errorChecks,
                                           all,
                                           userUploadExperimentInformation)
{
  exploratoryDataAnalysisServerComponents <- {
    # Clear unused memory
    gc()

    # Set all outputs to blank, this resets
    # all the visualizations to blank after clicking analyse
    output$table <- renderDataTable({})
    output$logTransformationText <- renderUI({})
    output$experimentInfo <- renderUI({})
    output$knnColumnTable <- renderDataTable({})
    output$boxAndWhiskerPlot <- renderUI({})
    output$interactiveDensityPlot <- renderPlotly({})
    output$interactiveThreeDDensityPlot <- renderPlotly({})
    output$interactiveUmapPlot <- renderPlotly({})
    output$interactiveHeatMapPlot <- renderPlotly({})
    output$interactiveMeanVariancePlot <- renderPlotly({})
    output$interactivePcaScreePlot <- renderPlotly({})
    output$interactivePcaIndividualsPlot <- renderPlotly({})
    output$interactivePcaVariablesPlot <- renderPlotly({})
    output$interactive3DPcaVariablesPlot <- renderPlotly({})
    output$dETable <- renderDataTable({})
    output$iDEHistogram <- renderPlotly({})
    output$dEVennDiagram <- renderPlot({})
    output$iDEQQ <- renderPlotly({})
    output$iDEVolcano <- renderPlotly({})
    output$iDEMd <- renderPlotly({})
    output$iHeatmap <- renderPlotly({})

    # Extract information from GSET including expression data
    if (errorChecks$continueWorkflow) {
      if (input$dataSource == "GEO") {
        # Extract the GEO data from the specified platform
        all$gsetData <- tryCatch({
          extractPlatformGset(all$allGset(), input$platform)
        }, error = function(err) {
          # Return null if there is a error in the getGeoObject function
          return(NULL)
        })

        # Error handling to prevent users
        # trying to run exploratory data analysis
        # without selecting a platform
        if (is.null(all$gsetData)) {
          # Update Error Checks
          errorChecks$geoPlatform <- FALSE
          errorChecks$continueWorkflow <- FALSE

          # Show error
          showNotification("Please select a platform.",
                           type = "error")
        } else
        {
          errorChecks$geoPlatform <- TRUE
          errorChecks$continueWorkflow <- TRUE

          # Extract expression data
          all$expressionData <-
            extractExpressionData(all$gsetData)

          # Extract the experiment information
          all$experimentInformation <-
            extractExperimentInformation(all$gsetData)

          # Convert experiment information to HTML
          all$convertedExperimentInformation <-
            convertExperimentInformation(all$experimentInformation)

          # Extract Column Information
          all$columnInfo <- extractSampleDetails(all$gsetData)
        }

      } else
      {
        # Error handling to prevent no file being uploaded
        if (is.null(input$file1)) {
          showNotification("Please ensure you have uploaded a file before
                             clicking analyse.",
                           type = "error")
        } else
        {
          # Error handling to prevent non-csvs being uploaded
          if (file_ext(input$file1$name) %in% c('text/csv',
                                                'text/comma-separated-values',
                                                'text/plain',
                                                'csv')) {
            # Update error checks
            errorChecks$uploadFile <- TRUE
            errorChecks$continueWorkflow <- TRUE
            # Ensure a file has been uploaded
            req(input$file1$datapath)
            # Extract Expression Data from CSV
            all$expressionData <- tryCatch({
              readCsvFile(input$file1$datapath)
            },
            error = function(e) {
              # return null if there is an error
              return(NULL)
            })

            # Define experimental information
            all$convertedExperimentInformation <-
              userUploadExperimentInformation

            # Preprocess the data
            all$expressionData <-
              tryCatch({
                preProcessGeneExpressionData(all$expressionData)
              },
              error = function(e) {
                # return null if there is an error
                return(NULL)
              })

            # Expression Error Check
            if (is.null(all$expressionData)) {
              # Update error checks
              errorChecks$expressionData <- FALSE
              errorChecks$continueWorkflow <- FALSE
            } else
            {
              # Update error checks
              errorChecks$expressionData <- TRUE
              errorChecks$continueWorkflow <- TRUE
              # Extract Column Information
              all$columnInfo <-
                convertExpressionDataToExperimentInformation(
                  all$expressionData)
            }
          } else
          {
            # Update error checks
            errorChecks$uploadFile <- FALSE
            errorChecks$continueWorkflow <- FALSE
            # Show notification
            showNotification(
              "The gene expression file does not have the correct
              file extension. Please upload a CSV.",
              type = "error"
            )
          }
        }
      }
    }

    # Combining datasets workflow
    if (input$dataSetType == "Combine") {
      if (input$dataSource2 == "GEO") {
        if (input$dataSource == "GEO") {
          if (input$platform != input$platform2) {
            showNotification(
              "The two GEO series platforms are not the same.
                               This might cause an error if the datasets do not
                               have the same row names."
              ,
              type = "warning"
            )
          }
        }
        if (errorChecks$continueWorkflow &
            errorChecks$continueWorkflow2) {
          # Extract the GEO data from the specified platform
          all$gsetData2 <- tryCatch({
            extractPlatformGset(all$allGset2(), input$platform2)
          }, error = function(err) {
            # Return null if there is a error in the getGeoObject function
            return(NULL)
          })

          # Error handling to prevent users
          # trying to run exploratory data analysis
          # without selecting a platform
          if (is.null(all$gsetData2)) {
            # Update Error Checks
            errorChecks$geoPlatform2 <- FALSE
            errorChecks$continueWorkflow2 <- FALSE

            # Show error
            showNotification("Please select a platform.",
                             type = "error")
          } else
          {
            errorChecks$geoPlatform2 <- TRUE
            errorChecks$continueWorkflow2 <- TRUE

            # Extract expression data
            all$expressionData2 <-
              extractExpressionData(all$gsetData2)

            # Extract the experiment information
            all$experimentInformation2 <-
              extractExperimentInformation(all$gsetData2)

            # Extract column information
            all$columnInfo2 <- extractSampleDetails(all$gsetData2)

            # Extract experiment information
            all$convertedExperimentInformation2 <-
              convertExperimentInformation(all$experimentInformation2)

            # Combine experiment information
            all$convertedExperimentInformation <-
              convertTwoExperimentInformation(
                all$convertedExperimentInformation,
                all$convertedExperimentInformation2
              )
          }
        }
      } else {
        # Error handling to prevent non-csvs being uploaded
        if (file_ext(input$file2$name) %in% c('text/csv',
                                              'text/comma-separated-values',
                                              'text/plain',
                                              'csv')) {
          # Update error checks
          errorChecks$uploadFile2 <- TRUE
          errorChecks$continueWorkflow2 <- TRUE
          # Ensure a file has been uploaded
          req(input$file2)
          # Extract Expression Data from CSV
          all$expressionData2 <- tryCatch({
            readCsvFile(input$file2$datapath)
          },
          error = function(e) {
            # return null if there is an error
            return(NULL)
          })

          # Preprocess the data
          all$expressionData2 <-
            tryCatch({
              preProcessGeneExpressionData(all$expressionData2)
            },
            error = function(e) {
              # return null if there is an error
              return(NULL)
            })

          # Expression Error Check
          if (is.null(all$expressionData2)) {
            # Update error checks
            errorChecks$expressionData2 <- FALSE
            errorChecks$continueWorkflow2 <- FALSE
          } else
          {
            # Update error checks
            errorChecks$expressionData2 <- TRUE
            errorChecks$continueWorkflow2 <- TRUE
            # Extract Column Information
            all$columnInfo2 <-
              convertExpressionDataToExperimentInformation(all$expressionData2)

            # Define Experimental Information
            all$convertedExperimentInformation2 <- HTML(
              "<b>Experimental
          Information is not available when processing
          user-uploaded files!</b>"
            )
          }
        } else {
          # Update error checks
          errorChecks$uploadFile2 <- FALSE
          errorChecks$continueWorkflow2 <- FALSE
          # Show notification
          showNotification(
            "The gene expression file does not have the correct
              file extension. Please upload a CSV.",
            type = "error"
          )
        }
      }

      # Combine the expression datasets
      combinedExpressionData <- tryCatch({
        # Combine the two dataframes
        combineExpressionData(all$expressionData, all$expressionData2)
      }, error = function(err) {
        # Return null if there is a error in the getGeoObject function
        return(NULL)
      })

      if (is.null(combinedExpressionData)) {
        # Show error
        showNotification(
          "The two gene expression datasets
                                 could not be merged. Please make sure
                                 they have the same platform.  Only the first
                                 gene expression datasets was
                                 processed as a result.",
          type = "warning"
        )
      } else
      {
        # Perform batch correction
        combinedExpressionDataBatchRemoved <- tryCatch({
          calculateBatchCorrection(
            all$expressionData,
            all$expressionData2,
            combinedExpressionData,
            input$batchCorrection
          )
        }, error = function(err) {
          # Return null if there is a error in the getGeoObject function
          return(NULL)
        })

        if (is.null(combinedExpressionDataBatchRemoved)) {
          # Show error
          showNotification(
            "There was an error performing batch
                                   correction. Therefore the non-batch
                                   corrected data was used.",
            type = "warning"
          )

          # Update expression data with non-batch corrected data
          all$expressionData <- combinedExpressionData

        } else
        {
          # Update expression data
          all$expressionData <-
            combinedExpressionDataBatchRemoved
        }
        # Combine experimental conditions
        all$columnInfo <-
          rbind(all$columnInfo , all$columnInfo2)
      }
    }

    # Process Expression Data
    if (errorChecks$continueWorkflow) {
      # Error handling to detect wrong format expression data
      if (!isNumeric(all$expressionData)) {
        errorChecks$expressionData <- FALSE
        errorChecks$continueWorkflow <- FALSE
        # Display error message
        showNotification(
          "The gene expression data has non-numerical values.
              Please ensure the gene expression data has only numerical values.
              ",
          type = "error"
        )
      } else if (length(all$expressionData) == 0)
      {
        errorChecks$expressionData <- FALSE
        errorChecks$continueWorkflow <- FALSE
        # Error handling to prevent issues
        # due to expression data with no samples
        showNotification(
          "The expression data is empty
            and therefore can not be analysed.
            This may indicate the GEO accession
            code relates to an RNA sequence experiment
                             rather than a microarray experiment.",
          type = "error"
        )
      }
      else if ((isNumeric(all$expressionData)) &
               ((length(all$expressionData) == 0) == FALSE)) {
        # Error handling to prevent errors caused by
        # expression datasets with only one column
        if (ncol(all$expressionData) <= 1) {
          # Update error check
          errorChecks$expressionDataOverOneColumns <- FALSE
          errorChecks$expressionDataOverTwoColumns <- FALSE
          # Display notification
          showNotification(
            "As the expression dataset had only one column only the
                Box-and-Whisper Plot and Expression Density Plots will be
                produced.",
            type = "warning"
          )
        } else if (ncol(all$expressionData) <= 2)
        {
          # Update error check
          errorChecks$expressionDataOverTwoColumns <- FALSE
          # Display notification
          showNotification(
            "As the gene expression data has less than 3 columns, the
                3D PCA Variables Plot will not be produced.",
            type = "warning"
          )
        }

        if (input$dataSource == "Upload") {
          if (input$typeOfData == "RNA Sequencing") {
            # Raw counts are converted to counts-per-million (CPM)
            all$cpm <- tryCatch({
              calculateCountsPerMillion(all$expressionData,
                                        input$cpmTransformation)
            },
            error = function(e) {
              # return null if there is an error
              return(NULL)
            })

            if (is.null(all$cpm)) {
              # Update cpm
              all$cpm <- all$expressionData

              showNotification(
                "There was an error calculating CPM. Therefore, the
                    original expression data will be used.",
                type = "warning"
              )
            }
          } else
          {
            all$cpm <- all$expressionData
          }
        } else
        {
          all$cpm <- all$expressionData
        }

        # Data Transformation Functions
        # Apply log transformation to expression
        #data if necessary
        all$dataInput <- tryCatch({
          calculateLogTransformation(all$cpm,
                                     input$logTransformation)
        }, error = function(cond) {
          return(NULL)
        })
        # Error handling to display a notification if
        # there was an error in log transformation
        if (is.null(all$dataInput)) {
          # Update error check
          errorChecks$dataInput <- FALSE

          # Display error notification
          showNotification(
            "There was an error applying log transformation to the
                expression data. Therefore, the original expression data
                will be used.",
            type = "warning"
          )
          all$dataInput <- all$cpm
        }
        # Is log transformation auto applied
        autoLogInformation <- tryCatch({
          calculateAutoLogTransformApplication(all$cpm)
        }, error = function(cond) {
          return(
            "There was an error calculating if log transformation
                       would automatically be applied."
          )
        })

        if (input$dataSource == "Upload") {
          if (input$typeOfData == "RNA Sequencing") {
            # Perform KNN transformation on log
            # expression data if necessary
            all$knnDataInput <- tryCatch({
              calculateKnnImpute(all$dataInput,
                                 "No")
            }, error = function(cond) {
              return(NULL)
            })
          } else
          {
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
          }
        } else {
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
        }

        # Error handling to display a notification if
        # there was an error in KNN imputation
        if (is.null(all$knnDataInput)) {
          # Update error check
          errorChecks$knnDataInput <- FALSE
          # Display notification
          showNotification(
            "There was an error applying KNN imputation to the
                expression data. Therefore, the log transformed/original
                expression data will be used.",
            type = "warning"
          )
          all$knnDataInput <- all$dataInput
        }


        # Download Gene Expression Dataset
        output$downloadGeneExpression <- try(downloadHandler(
          filename = "transformed_gene_expression_dataset.csv",
          content = function(file) {
            write.csv(all$knnDataInput,
                      file,
                      row.names = TRUE)
          }
        ))

        # KNN Column Set Plot
        all$knnColumns <-
          extractSampleNames(all$knnDataInput)

        # Update col info
        all$columnInfo <-
          all$columnInfo[all$knnColumns,]

        # Remove all incomplete rows
        naOmitInput <- calculateNaOmit(all$knnDataInput)

        # Perform PCA analysis on KNN transformation
        # expression data using princomp
        pcaPrcompDataInput <- tryCatch({
          calculatePrcompPca(naOmitInput)
        }, error = function(cond) {
          return(NULL)
        })
        # Error handling to display a notification
        # if there was an error in PCA
        if (is.null(pcaPrcompDataInput)) {
          # Update error check
          errorChecks$pcaPrcompDataInput <- FALSE
          # Display notification
          showNotification(
            "There was an error performing principal component analysis
                (PCA) on the expression data. Therefore, the PCA
                visualisations will not be displayed.",
            type = "warning"
          )
        }



      }
    }

    # Process Data Visualisations
    if (errorChecks$continueWorkflow) {
      # Experimental Information Display
      output$experimentInfo <- tryCatch({
        renderUI({
          all$convertedExperimentInformation
        })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })

      # Column Set Plot
      output$columnTable <-
        tryCatch({
          renderDataTable(all$columnInfo, selection = 'none')
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      # Expression dataset table
      output$table <-
        tryCatch({
          renderDataTable(all$knnDataInput, selection = 'none')
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
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

      output$knnColumnTableOne <- tryCatch({
        renderDataTable(all$columnInfo
                        ,
                        selection = 'multiple'
                        ,
                        server = FALSE)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })

      observeEvent(input$knnColumnTableOne_rows_selected, {
        all$knnColumnTableTwo <- all$columnInfo[
          -input$knnColumnTableOne_rows_selected,]

        output$knnColumnTableTwo <- tryCatch({
          renderDataTable(all$knnColumnTableTwo
                          ,
                          selection = 'multiple'
                          ,
                          server = FALSE)
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      })



      # Expression dataset table
      output$table <-
        tryCatch({
          renderDataTable(all$knnDataInput, selection = 'none')
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })

      if (object.size(all$knnDataInput) < 10000000) {
        output$boxAndWhiskerPlot <- renderUI({
          plotlyOutput('interactiveBoxAndWhiskerPlot')
          })

        # Interactive Box-and-Whisker Plot
        output$interactiveBoxAndWhiskerPlot <-
          tryCatch({
            renderPlotly({
              interactiveBoxAndWhiskerPlot(all$knnDataInput)
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })
      } else if (object.size(naOmitInput) < 10000000) {
        showNotification(
          "Due to the size of the dataset, the
                                     Box-and-Whisker plot was created after
                                     rows containing missing values were
                                     removed.",
          type = "warning"
        )
        output$boxAndWhiskerPlot <- renderUI({
          plotlyOutput('interactiveBoxAndWhiskerPlot')
        })

        # Interactive Box-and-Whisker Plot
        output$interactiveBoxAndWhiskerPlot <-
          tryCatch({
            renderPlotly({
              interactiveBoxAndWhiskerPlot(naOmitInput)
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })
      } else {
        showNotification(
          "Due to the size of the dataset, a static
                           Box-and-Whisker plot was created instead
                           of an interactive one.",
          type = "warning"
        )

        output$boxAndWhiskerPlot <- renderUI({
          plotOutput('nonInteractiveBoxAndWhiskerPlot')
        })

        output$nonInteractiveBoxAndWhiskerPlot <- tryCatch({
          renderPlot({
            nonInteractiveBoxAndWhiskerPlot(all$knnDataInput)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      }



      # Interactive Density Plot
      output$interactiveDensityPlot <-
        tryCatch({
          renderPlotly({
            interactiveDensityPlot(naOmitInput)
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
            interactiveThreeDDensityPlot(naOmitInput)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      # Error handling to display a notification
      # if there was an error in PCA
      if (errorChecks$pcaPrcompDataInput) {
        # Interactive PCA Scree Plot
        output$interactivePcaScreePlot <- tryCatch({
          renderPlotly({
            interactivePrcompPcaScreePlot(pcaPrcompDataInput)
          })
        }, error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })
      }
      # Error handling to prevent errors caused by
      # expression datasets with only one column
      if (errorChecks$expressionDataOverOneColumns) {
        # Update UMAP KNN max
        updateNumericInput(
          session,
          inputId = "knn",
          value = 2,
          max = ncol(all$cpm)
        )

        # Interactive UMAP Plot
        output$interactiveUmapPlot <-
          tryCatch({
            renderPlotly({
              interactiveUmapPlot(naOmitInput,
                                  input$knn)
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
                                          all$gsetData)
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })
        if (errorChecks$pcaPrcompDataInput) {
          # Interactive PCA Individual Plot
          output$interactivePcaIndividualsPlot <-
            tryCatch({
              renderPlotly({
                interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
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
                interactivePrcompPcaVariablesPlot(pcaPrcompDataInput)
              })
            },
            error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
            })

          # Only Display 3D PCA Variables Plot if there are more
          # than two experimental samples
          if (errorChecks$expressionDataOverTwoColumns) {
            # Interactive 3D PCA Variables Plot
            output$interactive3DPcaVariablesPlot <-
              tryCatch({
                renderPlotly({
                  interactive3DPrcompPcaVariablesPlot(pcaPrcompDataInput)
                })
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
              })
          }
        }
      }
    }

    if (errorChecks$continueWorkflow) {
      showNotification("Exploratory data analysis complete!",
                       type = "message")

      # Make Differential Gene Expression Action
      # Button Appear, this prevents users
      # trying to perform differential gene expression analysis
      # prior to exploratory data analysis
      output$output100 <- renderUI({
        actionButton("differentialExpressionButton", "Analyse")
      })
    }
  }
  return(exploratoryDataAnalysisServerComponents)
}

#' A Function to Return the GEO sourcing Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT renderDataTable
#' @importFrom utils write.csv object.size
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @importFrom stringr str_trim
#' @import markdown
#' @importFrom knitr knit
#' @author Guy Hunt
#' @noRd
loadGeoDataset <- function (input,
                            output,
                            session,
                            errorChecks,
                            all) {
  loadGeoDatasetServerComponents <- {
    # Get the GEO data for all platforms
    all$allGset <- reactive({
      tryCatch({
        # Error handling to ensure geoAccessionCode is populated
        req(input$geoAccessionCode)
        # Notify the user the GEO accession code
        # is not a GEO series accession code
        if (substr(str_trim(input$geoAccessionCode), 1, 3) != "GSE")
        {
          showNotification("Please input a GEO series accession code
                                 with the format GSEXXX",
                           type = "warning")
          return(NULL)
        } else {
          return(getGeoObject(input$geoAccessionCode))
        }
      }, error = function(err) {
        # Return null if there is a error in the getGeoObject function
        return(NULL)
      })
    })

    # Update error check
    if (is.null(all$allGset())) {
      # Update error check
      errorChecks$geoAccessionCode <- FALSE
      errorChecks$continueWorkflow <- FALSE
      if (input$geoAccessionCode != "") {
        # Display notification
        showNotification(
          "There was an error obtaining the GEO dataset.
                           Please ensure you entered the correct GEO Accession
                           Code.",
          type = "warning"
        )
      }
    } else {
      # Update error checks
      errorChecks$geoAccessionCode <- TRUE
      errorChecks$continueWorkflow <- TRUE
    }

    if (errorChecks$continueWorkflow) {
      # Get a list of all the platforms
      platforms <- reactive({
        extractPlatforms(all$allGset())
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
  }
  return(loadGeoDatasetServerComponents)
}


#' A Function to Return the Differential Gene Expression
#' Analysis Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT renderDataTable
#' @importFrom utils write.csv object.size
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @importFrom stringr str_trim
#' @import markdown
#' @importFrom knitr knit
#' @author Guy Hunt
#' @noRd
performDifferentialGeneExpressionAnalysis <- function (input,
                            output,
                            session,
                            errorChecks,
                            all,
                            ct,
                            exampleDataSet = FALSE) {
  differentialGeneExpressionAnalysisServerComponents <- {
    # Clear unused memory
    gc()

    # Set all differential gene expression
    # analysis outputs to blank, this resets
    # all the visualizations to blank after
    # clicking analyse
    output$dETable <- renderDataTable({})
    output$iDEHistogram <- renderPlotly({})
    output$dEVennDiagram <- renderPlot({})
    output$iDEQQ <- renderPlotly({})
    output$iDEVolcano <- renderPlotly({})
    output$iDEMd <- renderPlotly({})
    output$iHeatmap <- renderPlotly({})
    output$geneAnnotationTable <- renderDataTable({})

    if (errorChecks$continueWorkflow)
    { if (!exampleDataSet) {
      gsms <- tryCatch({
        calculateEachGroupsSamplesGsms(
          all$columnInfo,
          row.names(all$columnInfo[input$knnColumnTableOne_rows_selected, ]),
          row.names(
            all$knnColumnTableTwo[input$knnColumnTableTwo_rows_selected, ])
        )

      }, error = function(cond) {
        return(NULL)
      })
    } else {
      gsms <- "11110000"
    }
      # Error handling to prevent differential gene expression
      # analysis being performed before exploratory data analysis
      if (is.null(gsms)) {
        # Update error check
        errorChecks$continueWorkflow <- FALSE
        showNotification(
          "There was an error running differential gene expression
                  analysis. Please ensure you have performed exploratory data
                  analysis first.",
          type = "error"
        )
      } else
      {
        # Update error check
        errorChecks$continueWorkflow <- TRUE
      }

      if (errorChecks$continueWorkflow) {
        # Error handling to ensure at least one
        # group has two samples and the other group
        # has at least one sample
        if ((lengths(regmatches(gsms, gregexpr("0", gsms))) > 0 &
             lengths(regmatches(gsms, gregexpr("1", gsms))) > 1) |
            (lengths(regmatches(gsms, gregexpr("0", gsms))) > 1 &
             lengths(regmatches(gsms, gregexpr("1", gsms))) > 0)) {
          if (input$dataSource == "GEO") {
            if (input$dataSetType == "Single") {
              all$results <- tryCatch({
                calculateDifferentialGeneExpression(
                  gsms,
                  input$limmaPrecisionWeights,
                  input$forceNormalization,
                  all$gsetData,
                  all$knnDataInput,
                  input$dataSource
                )
              }
              , error = function(cond) {
                return(NULL)
              })
            } else
            {
              all$results <- tryCatch({
                calculateDifferentialGeneExpression(
                  gsms,
                  input$limmaPrecisionWeights,
                  input$forceNormalization,
                  all$gsetData,
                  all$knnDataInput,
                  input$dataSource,
                  NULL,
                  input$dataSetType
                )
              }
              , error = function(cond) {
                return(NULL)
              })
            }
          } else
          {
            all$results <- tryCatch({
              calculateDifferentialGeneExpression(
                gsms,
                input$limmaPrecisionWeights,
                input$forceNormalization,
                all$gsetData,
                all$knnDataInput,
                input$dataSource,
                input$typeOfData
              )
            }
            , error = function(cond) {
              return(NULL)
            })
            if (is.null(all$results)) {
              if (input$typeOfData == "RNA Sequencing") {
                # Try again with non-log data
                knnDataInput <-
                  calculateKnnImpute(all$cpm, input$knnTransformation)
                all$results <- tryCatch({
                  calculateDifferentialGeneExpression(
                    gsms,
                    input$limmaPrecisionWeights,
                    input$forceNormalization,
                    all$gsetData,
                    knnDataInput,
                    input$dataSource,
                    input$typeOfData
                  )
                }
                , error = function(cond) {
                  return(NULL)
                })
                # Show warning that non-log data was used
                showNotification(
                  "There was an error calculating the
                             differential gene expression analysis
              using the log data. So the non-log data was used instead!",
                  type = "warning"
                )
              }
            }

          }


          # Error handling to ensure Differential Gene
          # Expression Analysis worked
          if (is.null(all$results)) {
            # Update Error Check
            errorChecks$continueWorkflow <- FALSE
            # Display notification
            showNotification(
              "There was an error calculating the
                 differential gene expression analysis!",
              type = "error"
            )
          } else
          {
            # Update error check
            errorChecks$continueWorkflow <- TRUE
          }
        } else {
          # Update error check
          errorChecks$continueWorkflow <- FALSE
          # Display notification
          showNotification(
            "One group needs at least 2 samples and the other
                              group needs at least 1 sample",
            type = "error"
          )

        }

      }
    }

    if (errorChecks$continueWorkflow) {
      all$adjustment <- convertAdjustment(input$pValueAdjustment)
      all$tT <-
        calculateTopDifferentiallyExpressedGenes(all$results$fit2,
                                                 all$adjustment,
                                                 nrow(all$knnDataInput))

      all$significanceLevelCutOff <- input$significanceLevelCutOff

      all$dT <-
        calculateDifferentialGeneExpressionSummary(all$results$fit2,
                                                   all$adjustment,
                                                   all$significanceLevelCutOff)
      # Differential gene expression table
      output$dETable <- tryCatch({
        renderDataTable(as.data.frame(all$tT), selection = 'none')
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })

      # Interactive Histogram Plot
      output$iDEHistogram <- tryCatch({
        renderPlotly({
          interactiveHistogramPlot(all$results$fit2, all$adjustment)
        })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })


      # Venn Diagram Plot
      output$dEVennDiagram <- tryCatch({
        renderPlot({
          nonInteractiveVennDiagramPlot(all$dT)
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
            interactiveQQPlot(all$results$fit2, all$dT, ct)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })


      # Interactive Volcano Plot
      output$iDEVolcano <- tryCatch({
        renderPlotly({
          interactiveVolcanoPlot(all$results$fit2, all$dT, ct)
        })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })


      # Interactive Mean Difference Plot
      output$iDEMd <- tryCatch({
        renderPlotly({
          interactiveMeanDifferencePlot(all$results$fit2, all$dT, ct)
        })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      })

      # Update Interactive Heatmap Plot Gene number max
      updateNumericInput(session, "numberOfGenes", max = nrow(all$tT))

      # Interactive Heatmap Plot
      output$iHeatmap <-
        tryCatch({
          renderPlotly({
            interactiveDGEHeatMapPlot(all$results$ex,
                                      input$limmaPrecisionWeights,
                                      input$numberOfGenes,
                                      all$tT)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })

      # Download Top Differentially Expressed Genes Table
      output$downloadData <- try(downloadHandler(
        filename = "top_differentially_expressed_genes.csv",
        content = function(file) {
          write.csv(all$tT,
                    file,
                    row.names = TRUE)
        }
      ))

      showNotification("Differential gene
                           expression analysis complete!",
                       type = "message")

      # Make Differential Gene Expression Action
      # Button Appear, this prevents users
      # trying to perform differential gene expression analysis
      # prior to exploratory data analysis
      output$output101 <- renderUI({
        actionButton("enrichmentAnalysisButton", "Analyse")
      })

      # Define Gene Annotation Table
      all$geneAnnotationTable <- tryCatch({
        differentiallyExpressedGenes <-
          all$dT[!(all$dT[, "Group1-Group2"] == 0),]

        differentiallyExpressedGenes <- as.data.frame(
          differentiallyExpressedGenes)

        try({differentiallyExpressedGenes[, "Gene.symbol"] <- NA
        differentiallyExpressedGenes <-
          differentiallyExpressedGenes[, c(2,1)]})

        differentiallyExpressedGenes
      }, error = function(e) {
        # return a safeError if a parsing error occurs
        return(NULL)})

      if(input$dataSource == "GEO") {
        all$geneAnnotationTable <- try({
          createGeneAnnotationTable(
          input, output, session, errorChecks, all)
        })

      } else if (input$dataSetType == "Combine") {
        if (input$dataSource2 == "GEO") {
          all$geneAnnotationTable <- try({
            createGeneAnnotationTable(
              input, output, session, errorChecks, all)
          })
        }}

      output$geneAnnotationTable <-
        tryCatch({
          renderDataTable(
            all$geneAnnotationTable,
            server = FALSE,
            escape = FALSE,
            options = list(pageLength = 3),
            editable = TRUE,
            selection = list(target = 'column', mode = "single")
          )},
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })

      all$updatedGeneAnnotationTable <- all$geneAnnotationTable

      # Update gene annotation data with user input values
      observeEvent(input$geneAnnotationTable_cell_edit, {
        info <- input$geneAnnotationTable_cell_edit
        try(all$updatedGeneAnnotationTable[info$row,info$col] <- info$value)
        })

    }
    # Reset error check
    errorChecks$continueWorkflow <- TRUE
  }
  return(differentialGeneExpressionAnalysisServerComponents)
}

#' A Function to Return the Gene Enrichment Analysis Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @importFrom DT renderDataTable
#' @importFrom utils write.csv object.size
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @importFrom stringr str_trim
#' @import markdown
#' @importFrom knitr knit
#' @author Guy Hunt
#' @noRd
performGeneEnrichmentAnalysis <- function (input,
                                           output,
                                           session,
                                           errorChecks,
                                           all,
                                           databaseNames,
                                           columnNumber = NULL
                                           ) {
  geneEnrichmentAnalysisServerComponents <- {
    # Reset all Visualisations
    output$differentiallyExpressedGenesEnrichmentTable <- renderDataTable({})
    output$differentiallyExpressedGenesEnrichmentPlot <- renderPlot({})
    output$upregulatedGenesEnrichmentTable <- renderDataTable({})
    output$upregulatedGenesEnrichmentPlot <- renderPlot({})
    output$downregulatedGenesEnrichmentTable <- renderDataTable({})
    output$downregulatedGenesEnrichmentPlot <- renderPlot({})

    if(errorChecks$continueWorkflow) {
      if (!is.null(databaseNames)) {
        if (is.null(columnNumber)) {
          if (is.null(input$geneAnnotationTable_columns_selected)) {
            showNotification("Please select the gene symbol column.",
                             type = "error")
          } else {
            columnNumber <- try(input$geneAnnotationTable_columns_selected)
          }
        }

        if (!is.null(columnNumber)) {
          differentiallyExpressedGenes <- tryCatch({
            differentiallyExpressedGenes <- all$updatedGeneAnnotationTable[
              ,c(columnNumber,ncol(all$updatedGeneAnnotationTable))]
            colnames(differentiallyExpressedGenes) <- c("Gene.symbol",
                                                        "Group1-Group2")
            differentiallyExpressedGenes
          }, error = function(e) {
            # return a safeError if a parsing error occurs
            return(NULL)
          })

          # Error handling if there are no differentially expressed genes
          if (is.null(differentiallyExpressedGenes) |
              nrow(differentiallyExpressedGenes)==0) {
            showNotification("There are no differentially expressed genes.
                           Therefore, enrichment analysis will not be
                           performed.", type = "error")
          } else
          {
            # Analyse Differential Expressed Genes
            # Extract differentially expressed gene symbols
            differemtiallyExpressedGeneSymbols <- tryCatch({
              extractGeneSymbols(differentiallyExpressedGenes, "Gene.symbol")
            }, error = function(e) {
              # return a safeError if a parsing error occurs
              return(NULL)
            })

            # enrich Differentially Expressed Genes
            enrichedDifferentiallyExpressedGenes <- tryCatch({
              enrichGenes(differemtiallyExpressedGeneSymbols,
                          input$enrichDatabases)
            }, error = function(e) {
              # return a safeError if a parsing error occurs
              return(NULL)
            })

            enrichedDifferentiallyExpressedGenes <- tryCatch({
              calculateLogPValue(enrichedDifferentiallyExpressedGenes)},
              error = function(e) {
                # return a safeError if a parsing error occurs
                return(enrichedDifferentiallyExpressedGenes)
              })

            enrichedDifferentiallyExpressedGenes <- tryCatch({
              calculateOverlapFractions(
                enrichedDifferentiallyExpressedGenes)}, error = function(e) {
                # return a safeError if a parsing error occurs
                return(enrichedDifferentiallyExpressedGenes)
              })

            # Extract Upregulated genes
            upregulatedGenes <- tryCatch({
              extractUpregulatedGenes(differentiallyExpressedGenes)
            }, error = function(e) {
              # return a safeError if a parsing error occurs
              return(NULL)
            })

            if (is.null(upregulatedGenes) | nrow(upregulatedGenes) == 0) {
              showNotification("There are no upregulated genes.
                           Therefore, upregulated gene enrichment analysis will
                           not be performed.", type = "warning")

              enrichedUpregulatedGenes <- NULL
            } else {
              # Extract upregulated gene symbols
              upregulatedGenesGeneSymbols <- tryCatch({
                extractGeneSymbols(upregulatedGenes, "Gene.symbol")
              }, error = function(e) {
                # return a safeError if a parsing error occurs
                return(NULL)
              })

              # enrich upregulated Genes
              enrichedUpregulatedGenes <- tryCatch({
                enrichGenes(upregulatedGenesGeneSymbols, input$enrichDatabases)
              }, error = function(e) {
                # return a safeError if a parsing error occurs
                return(NULL)
              })

              enrichedUpregulatedGenes <- tryCatch({
                calculateLogPValue(enrichedUpregulatedGenes)},
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  return(enrichedUpregulatedGenes)
                })

              enrichedUpregulatedGenes <- tryCatch({
                calculateOverlapFractions(
                  enrichedUpregulatedGenes)}, error = function(e) {
                    # return a safeError if a parsing error occurs
                    return(enrichedUpregulatedGenes)
                  })
            }


            # Extract downregulated genes
            downregulatedGenes <- tryCatch({
              extractdowregulatedGenes(differentiallyExpressedGenes)
            }, error = function(e) {
              # return a safeError if a parsing error occurs
              return(NULL)
            })

            # Error handling for no downregulated genes
            if (is.null(downregulatedGenes) | nrow(downregulatedGenes) == 0) {
              showNotification("There are no downregulated genes.
                           Therefore, downregulated gene
                           enrichment analysis will
                           not be performed.", type = "warning")
              enrichedDownregulatedGenes <- NULL
            } else {
              # Extract downregulated gene symbols
              downregulatedGenesGeneSymbols <- tryCatch({
                extractGeneSymbols(downregulatedGenes, "Gene.symbol")
              }, error = function(e) {
                # return a safeError if a parsing error occurs
                return(NULL)
              })

              # enrich downregulated Genes
              enrichedDownregulatedGenes <- tryCatch({
                enrichGenes(downregulatedGenesGeneSymbols,
                            input$enrichDatabases)
              }, error = function(e) {
                # return a safeError if a parsing error occurs
                return(NULL)
              })

              enrichedDownregulatedGenes <- tryCatch({
                calculateLogPValue(enrichedDownregulatedGenes)},
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  return(enrichedDownregulatedGenes)
                })

              enrichedDownregulatedGenes <- tryCatch({
                calculateOverlapFractions(
                  enrichedDownregulatedGenes)}, error = function(e) {
                    # return a safeError if a parsing error occurs
                    return(enrichedDownregulatedGenes)
                  })
            }

            observeEvent(input$geneEnrichmentDataBarchartPlot,{
              updateSelectInput(session, "geneEnrichmentDataManhattanPlot",
                                selected =
                                  input$geneEnrichmentDataBarchartPlot)
              updateSelectInput(session, "geneEnrichmentDataVolcanoPlot",
                                selected =
                                  input$geneEnrichmentDataBarchartPlot)

              updateSelectInput(session, "geneEnrichmentDataTable",
                                selected =
                                  input$geneEnrichmentDataBarchartPlot)

              if (input$geneEnrichmentDataBarchartPlot == "All differentially
                                         expressed genes") {
                enrichedGenes <- enrichedDifferentiallyExpressedGenes
              } else if (input$geneEnrichmentDataBarchartPlot ==
              "Upregulated genes") {
                enrichedGenes <- enrichedUpregulatedGenes
              } else {
                enrichedGenes <- enrichedDownregulatedGenes
              }

              # Display table of differentially expressed genes
              output$differentiallyExpressedGenesEnrichmentTable <- tryCatch({
                renderDataTable(
                  enrichedGenes,
                  server = FALSE,
                  escape = FALSE,
                  selection = 'none'
                )},
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError(e))
                })

              # Download differentially expressed gene enrichment
              output$downloadDifferentiallyExpressedGenesEnrichmentTable <-
                try(downloadHandler(
                  filename = "gene_enrichment.csv",
                  content = function(file) {
                    write.csv(enrichedGenes,
                              file,
                              row.names = TRUE)
                  }
                ))

              # Plot Differentially Expressed Genes
              output$differentiallyExpressedGenesEnrichmentPlot <- tryCatch({
                renderPlot({
                  plotGeneEnrichmentinformation(enrichedGenes)})},
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError(e))
                })

              output$genesEnrichmentVolcanoPlot <- tryCatch({renderPlotly({
                interactiveGeneEnrichmentVolcanoPlot(
                  enrichedGenes)
              })},
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
              })

              output$genesEnrichmentManhattanPlot <- tryCatch({renderPlotly({
                interactiveGeneEnrichmentManhattanPlot(
                  enrichedGenes)
              })},
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
              })

              reloadBarChart <-
                reactive(
                  c(
                    input$sortDecreasingly,
                    input$columnToSort,
                    input$recordsToDisplay,
                    input$geneEnrichmentDataBarchartPlot
                  )
                )

              observeEvent(reloadBarChart(),
                           {sortDecreasingly <- convertUiSortingMethod(
                             input$sortDecreasingly)

                           sortedEnrichedDifferentiallyExpressedGenes <- try(
                             sortGeneEnrichmentTable(enrichedGenes,
                                                     input$columnToSort,
                                                     sortDecreasingly))

                           topSortedEnrichedDifferentiallyExpressedGenes <-
                             try(selectTopGeneEnrichmentRecords(
                               sortedEnrichedDifferentiallyExpressedGenes,
                               input$recordsToDisplay))

                           output$genesEnrichmentBarchartPlot <- tryCatch({
                             renderPlotly({
                             interactiveGeneEnrichmentBarPlot(
                               topSortedEnrichedDifferentiallyExpressedGenes,
                               input$columnToSort)
                           })},
                           error = function(e) {
                             # return a safeError if a parsing error occurs
                             stop(safeError(e))
                           })})

              updateSliderInput(
                session,
                "recordsToDisplay",
                max = nrow(enrichedGenes)
              )


            })

            observeEvent(input$geneEnrichmentDataManhattanPlot,{
              updateSelectInput(session, "geneEnrichmentDataBarchartPlot",
                                selected =
                                  input$geneEnrichmentDataManhattanPlot)
              updateSelectInput(session, "geneEnrichmentDataVolcanoPlot",
                                selected =
                                  input$geneEnrichmentDataManhattanPlot)
              updateSelectInput(session, "geneEnrichmentDataTable",
                                selected =
                                  input$geneEnrichmentDataManhattanPlot)
            })

            observeEvent(input$geneEnrichmentDataVolcanoPlot,{
              updateSelectInput(session, "geneEnrichmentDataBarchartPlot",
                                selected =
                                  input$geneEnrichmentDataVolcanoPlot)
              updateSelectInput(session, "geneEnrichmentDataManhattanPlot",
                                selected =
                                  input$geneEnrichmentDataVolcanoPlot)
              updateSelectInput(session, "geneEnrichmentDataTable",
                                selected =
                                  input$geneEnrichmentDataVolcanoPlot)
            })

            observeEvent(input$geneEnrichmentDataTable,{
              updateSelectInput(session, "geneEnrichmentDataBarchartPlot",
                                selected =
                                  input$geneEnrichmentDataTable)
              updateSelectInput(session, "geneEnrichmentDataManhattanPlot",
                                selected =
                                  input$geneEnrichmentDataTable)
              updateSelectInput(session, "geneEnrichmentDataVolcanoPlot",
                                selected =
                                  input$geneEnrichmentDataTable)
            })

            showNotification("Gene enrichment analysis completed!",
                             type = "message")
          }
        }

      } else
        {
        showNotification("The enrichment website enrichR
                         is unavailable. Therefore, enrichment
                         analysis can not be performed at this
                         time.",type = "error")
      }
    }
  }
  return(geneEnrichmentAnalysisServerComponents)
}

#' A Function to load the UI for different datasets
#'
#' A Function to load the UI for different datasets
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @author Guy Hunt
#' @noRd
loadDataSetUiComponents <- function(input,
                                    output,
                                    session,
                                    errorChecks,
                                    all,
                                    userUploadExperimentInformation) {
  dataSetUiComponents <- {
    # Refresh error checks
    errorChecks <- resetErrorChecks(errorChecks)
    # Update UI side bar with GEO widgets
    if (input$dataSource == "GEO") {
      # GEO help text
      output$output4 <- renderUI({
        helpText(
          "Input a GEO series accession code (GSEXXXX format)
      to examine the gene expression data.
      This can be obtained from https://www.ncbi.nlm.nih.gov/gds."
        )
      })
      # GEO accession input
      output$output5 <- renderUI({
        textInput("geoAccessionCode", "GEO accession code", "")
      })
      # Platform
      output$output6 <- renderUI({
        selectInput("platform", "Platform", c())
      })
      # KNN Imputation Radio Button
      output$output13 <- renderUI({
        radioButtons(
          "knnTransformation",
          label = "Apply k-nearest neighbors (KNN) algorithm to predict
      null data:",
          choices = list("Yes", "No"),
          selected = "No"
        )
      })

      resetErrorChecksVariables <-
        reactive(
          c(
            input$logTransformation,
            input$platform,
            input$knnTransformation,
            input$geoAccessionCode
          )
        )

      # Reset error checks when input variables are updated
      observeEvent(resetErrorChecksVariables(), {
        # Reset error checks
        errorChecks <- resetErrorChecks(errorChecks)
      })

      observeEvent(
        input$geoAccessionCode,
        loadGeoDataset(input,
                       output,
                       session,
                       errorChecks,
                       all)
      )
    } else
    {
      # Define variables
      all$gsetData <- NULL

      # Update variables if combining the dataset with a GEO
      # Dataset
      if (input$dataSetType == "Combine") {
        if (input$dataSource2 == "GEO") {
          all$gsetData <- all$gsetData2
        }
      } else {
        all$convertedExperimentInformation2 <-
          userUploadExperimentInformation
      }

      # Update UI side bar with User Upload widgets
      observeEvent(input$dataSetType, {
        # Reset error checks when data set type is changed
        errorChecks <- resetErrorChecks(errorChecks)

        # Microarray vs RNA Seq Widget
        if (input$dataSetType == "Combine") {
          reactiveDataSources <- reactive(c(input$dataSource,
                                            input$dataSource2))
          observeEvent(reactiveDataSources(), {
            if ((input$dataSource == "GEO") | (input$dataSource2 == "GEO"))
            {
              output$output4 <- renderUI({
                radioButtons(
                  "typeOfData",
                  label = "Is the data from Microarray or RNA Sequencing?",
                  choices = list("Microarray"),
                  selected = "Microarray"
                )
              })
            } else {
              output$output4 <- renderUI({
                radioButtons(
                  "typeOfData",
                  label = "Is the data from Microarray or RNA Sequencing?",
                  choices = list("Microarray", "RNA Sequencing"),
                  selected = "Microarray"
                )
              })
            }
          })
        } else
        {
          output$output4 <- renderUI({
            radioButtons(
              "typeOfData",
              label = "Is the data from Microarray or RNA Sequencing?",
              choices = list("Microarray", "RNA Sequencing"),
              selected = "Microarray"
            )
          })
        }
      })

      # File Upload Widget
      output$output5 <- renderUI({
        fileInput(
          "file1",
          "Upload CSV Gene Expression Count File",
          multiple = TRUE,
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        )
      })

      # Reset error checks when a new file is uploaded
      observeEvent(input$file1, {
        # Reset error checks
        errorChecks <- resetErrorChecks(errorChecks)
      })

      # Blank Widgets
      output$output6 <- renderUI({
      })
    }


    # Add or remove CPM radio button
    observeEvent(input$typeOfData, {
      # Define error checks
      errorChecks <- resetErrorChecks(errorChecks)

      if (input$typeOfData == "RNA Sequencing") {
        # Add CPM widget if the dataset is microarray
        output$output13 <- renderUI({
          radioButtons(
            "cpmTransformation",
            label = "Convert data to count per million:",
            choices = list("Yes", "No"),
            selected = "No"
          )
        })
        # Reset error checks when CPM transformation is updated
        observeEvent(input$cpmTransformation, {
          # Reset error checks
          errorChecks <- resetErrorChecks(errorChecks)
        })

      } else
      {
        # Add KNN Imputation if the dataset is microarray
        output$output13 <- renderUI({
          radioButtons(
            "knnTransformation",
            label = "Apply k-nearest neighbors (KNN) algorithm to predict
      null data:",
            choices = list("Yes", "No"),
            selected = "No"
          )
        })

        # Reset error checks when KNN transformation is updated
        observeEvent(input$knnTransformation, {
          # Reset error checks
          errorChecks <- resetErrorChecks(errorChecks)
        })
      }
    })
  }
  return(dataSetUiComponents)
}

#' A Function to load the UI when using different numbers of datasets
#'
#' A Function to load the UI when using different numbers of datasets
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @author Guy Hunt
#' @noRd
loadDataSetCombinationUiComponents <- function(input, output, session,
                                    errorChecks, all, geoAccessionCode = "") {
  dataSetCombinationUiComponents <- {
    # Refresh error checks
    errorChecks <- resetErrorChecks(errorChecks)

    if (input$dataSetType == "Combine"){
      # Define error checks
      errorChecks$continueWorkflow2 <- TRUE
      errorChecks$geoAccessionCode2 <- TRUE
      errorChecks$geoMicroarrayAccessionCode2 <- TRUE
      errorChecks$geoPlatform2 <- TRUE
      errorChecks$expressionData2 <- TRUE
      errorChecks$dataInput2 <- TRUE
      errorChecks$knnDataInput2 <- TRUE
      errorChecks$pcaPrcompDataInput2 <- TRUE
      errorChecks$expressionDataOverTwoColumns2 <- TRUE
      errorChecks$expressionDataOverOneColumns2 <- TRUE
      errorChecks$differentialGeneExpression2 <- TRUE
      errorChecks$differentialGeneExpressionGroup2 <- TRUE
      errorChecks$uploadFile2 <- TRUE
      errorChecks$uploadFileExtension2 <- TRUE
      errorChecks$uploadLogData2 <- TRUE

      # First Data Set Information Widget
      output$output2 <- renderUI({
        HTML(
          "<b>First Gene Expression Dataset Information</b><br></br>"
        )
      })
      # Second Data Set Information Widget
      output$output7 <- renderUI({
        HTML(
          "<b>Second Gene Expression Dataset Information</b><br></br>"
        )
      })
      # Second Data Source Widget
      output$output8 <- renderUI({
        radioButtons(
          "dataSource2",
          label = "Would you like to upload the gene expression data
      or source the data from GEO?",
          choices = list("GEO", "Upload"),
          selected = "GEO"
        )})

      # Second Data Source Widget
      output$output14 <- renderUI({
        radioButtons(
          "batchCorrection",
          label = "Batch correction method:",
          choices = list("Empirical Bayes", "Linear Model", "None"),
          selected = "None"
        )})

      inputErrorCheckVariables <- reactive(c(input$batchCorrection,
                                             input$dataSource2))

      # Reset error checks when input variables are updated
      observeEvent(inputErrorCheckVariables(), {
        # Reset error checks
        errorChecks <- resetErrorChecks(errorChecks)
      })

      # Load UI components
      observeEvent(input$dataSource2, loadDataSet2UiComponents(input,
                                                              output,
                                                              session,
                                                              errorChecks,
                                                              all,
                                                              geoAccessionCode)
                   )
    } else
    {
      # Set all UI widgets to blank
      output$output2 <- renderUI({})
      output$output7 <- renderUI({})
      output$output8 <- renderUI({})
      output$output9 <- renderUI({})
      output$output10 <- renderUI({})
      output$output11 <- renderUI({})
      output$output14 <- renderUI({})
    }
  }
  return(dataSetCombinationUiComponents)
}

#' A Function to load the UI for different datasets
#'
#' A Function to load the UI for different datasets
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @author Guy Hunt
#' @noRd
loadDataSet2UiComponents <- function(input, output, session, errorChecks, all,
                                     geoAccessionCode = "")
{
  dataSet2UiComponents <- {
    if (input$dataSource2 == "GEO") {
      # GEO Help Text Widget
      output$output9 <- renderUI({
        helpText(
          "Input a GEO series accession code (GSEXXXX format)
      to examine the gene expression data.
      This can be obtained from https://www.ncbi.nlm.nih.gov/gds."
        )
      })
      # GEO Accession Code Input Widget
      output$output10 <- renderUI({
        textInput("geoAccessionCode2", "GEO accession code", geoAccessionCode)
      })

      # Platform input text
      output$output11 <- renderUI({
        selectInput("platform2", "Platform", c())
      })
      # KNN Input
      output$output13 <- renderUI({
        radioButtons(
          "knnTransformation",
          label = "Apply k-nearest neighbors (KNN) algorithm to predict
      null data:",
          choices = list("Yes", "No"),
          selected = "No"
        )
      })

      resetErrorChecksVariables2 <-
        reactive(
          c(
            input$logTransformation,
            input$platform2,
            input$knnTransformation,
            input$geoAccessionCode2
          )
        )

      # Reset error checks when Platform is updated
      observeEvent(resetErrorChecksVariables2(), {
        # Reset error checks
        errorChecks <- resetErrorChecks(errorChecks)
      })

      # Process second GEO accession code
      observeEvent(input$geoAccessionCode2, {
        # Get the GEO data for all platforms
        all$allGset2 <- reactive({
          tryCatch({
            # Error handling to ensure geoAccessionCode is populated
            req(input$geoAccessionCode2)
            # Notify the user the GEO accession code
            # is not a GEO series accession code
            if (substr(str_trim(input$geoAccessionCode2), 1, 3) != "GSE")
            {
              showNotification("Please input a GEO series accession code
                                 with the format GSEXXX",
                               type = "warning")
              return(NULL)
            } else {
              return(getGeoObject(input$geoAccessionCode2))
            }
          }, error = function(err) {
            # Return null if there is a error in the
            # getGeoObject function
            return(NULL)
          })
        })

        # Update error check
        if (is.null(all$allGset2())) {
          # Update error check
          errorChecks$geoAccessionCode2 <- FALSE
          errorChecks$continueWorkflow2 <- FALSE
          if (input$geoAccessionCode2 != "") {
            # Display notification
            showNotification(
              "There was an error obtaining the GEO dataset.
                           Please ensure you entered the correct GEO Accession
                           Code.",
              type = "warning"
            )
          }
        } else {
          # Update error checks
          errorChecks$geoAccessionCode2 <- TRUE
          errorChecks$continueWorkflow2 <- TRUE
        }

        if (errorChecks$continueWorkflow2) {
          # Get a list of all the platforms
          platforms2 <- reactive({
            extractPlatforms(all$allGset2())
          })

          # Select the top platform
          platform2 <- reactive({
            platforms2()[1]
          })

          # Update Platform Options
          platformObserve2 <- observe({
            updateSelectInput(session,
                              "platform2",
                              choices = platforms2(),
                              selected = platform2())
          })
        }
      })
    } else
    {
      # File upload widget
      output$output9 <- renderUI({
        fileInput(
          "file2",
          "Upload CSV Gene Expression Count File",
          multiple = TRUE,
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        )
      })

      # Reset error checks when a file is uploaded
      observeEvent(input$file2, {
        # Reset error checks
        errorChecks <- resetErrorChecks(errorChecks)
      })

      # Blank widgets
      output$output10 <- renderUI({
      })
      output$output11 <- renderUI({
      })
    }
  }
  return(dataSet2UiComponents)
}

#' A Function to create the gene annotation table
#'
#' A Function to create the gene annotation table
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @author Guy Hunt
#' @noRd
createGeneAnnotationTable <- function(input, output, session, errorChecks, all)
  {
  # Output gene annotation table
  geneAnnotation <- tryCatch({fData(all$gsetData)},
                             error = function(e) {
                               # return a safeError if a parsing
                               # error occurs
                               return(NULL)})

  # Extract the differentially Expressed Gene Annotation
  differentiallyExressedGeneAnnotation <- tryCatch({
    extractDifferenitallyExpressedGenes(geneAnnotation, all$dT)},
    error = function(e) {
      # return a safeError if a parsing
      # error occurs
      return(NULL)})
}

