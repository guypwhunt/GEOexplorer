#' A Function to Return the Server Component
#'
#' A Function to Return the Server Component
#' @rawNamespace import(shiny, except = c(dataTableOutput, renderDataTable))
#' @examples sourceServer()
#' @importFrom DT renderDataTable JS
#' @importFrom shinyBS addTooltip
#' @importFrom utils write.csv
#' @importFrom htmltools HTML
#' @importFrom xfun file_ext
#' @author Guy Hunt
#' @noRd
sourceServer <- function(input, output, session) {
  datasetInformationServer <- ({
    ###############
    # Common steps
    # Define variables
    all <- reactiveValues()
    errorChecks <- reactiveValues()
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
              mean per sample. Columns with over 80% will cause an error in
              the KNN computation. This is only desinged to be used on
              microarray gene expression datasets.",
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

    # Define error checks
    errorChecks$continueWorkflow <- TRUE
    errorChecks$geoAccessionCode <- TRUE
    errorChecks$geoMicroarrayAccessionCode <- TRUE
    errorChecks$geoPlatform <- TRUE
    errorChecks$expressionData <- TRUE
    errorChecks$dataInput <- TRUE
    errorChecks$knnDataInput <- TRUE
    errorChecks$pcaPrcompDataInput <- TRUE
    errorChecks$expressionDataOverTwoColumns <- TRUE
    errorChecks$expressionDataOverOneColumns <- TRUE
    errorChecks$differentialGeneExpression <- TRUE
    errorChecks$differentialGeneExpressionGroup <- TRUE
    errorChecks$uploadFile <- TRUE
    errorChecks$uploadFileExtension <- TRUE
    errorChecks$uploadLogData <- TRUE

    ###############

    observeEvent(input$dataSource, {
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

        # Add platform tool tip
        addTooltip(
          session,
          id = "platform",
          title = "Each platform relates to a different microarray experiment
      performed in the study.",
          placement = "top",
          trigger = "hover"
        )

        # Add knn tool tip
        addTooltip(
          session,
          id = "knnTransformation",
          title = "Rows with over 50% missing values are
          imputed using the overall
              mean per sample. Columns with over
              80% will cause an error in the KNN
              computation.",
          placement = "top",
          trigger = "hover"
        )


        observeEvent(input$geoAccessionCode, {
          # Define error checks
          errorChecks$continueWorkflow <- TRUE
          errorChecks$geoAccessionCode <- TRUE
          errorChecks$geoMicroarrayAccessionCode <- TRUE
          errorChecks$geoPlatform <- TRUE
          errorChecks$expressionData <- TRUE
          errorChecks$dataInput <- TRUE
          errorChecks$knnDataInput <- TRUE
          errorChecks$pcaPrcompDataInput <- TRUE
          errorChecks$expressionDataOverTwoColumns <- TRUE
          errorChecks$expressionDataOverOneColumns <- TRUE
          errorChecks$differentialGeneExpression <- TRUE
          errorChecks$differentialGeneExpressionGroup <- TRUE
          errorChecks$uploadFile <- TRUE
          errorChecks$uploadFileExtension <- TRUE
          errorChecks$uploadLogData <- TRUE

          # Get the GEO data for all platforms
          all$allGset <- reactive({
            tryCatch({
              # Error handling to ensure geoAccessionCode is populated
              req(input$geoAccessionCode)
              # Notify the user the GEO accession code
              # is not a GEO series accession code
              if (substr(input$geoAccessionCode, 1, 3) != "GSE")
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

          if (errorChecks$continueWorkflow == TRUE) {
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
        })
      } else if
      (input$dataSource == "Upload") {
        # Define variables
        all$gsetData <- NULL

        # Update UI side bar with User Upload widgets
        observeEvent(input$dataSetType, {
          # Microarray vs RNA Seq Widget
          if (input$dataSetType == "Combine") {
            output$output4 <- renderUI({
              radioButtons(
                "typeOfData",
                label = "Is the data from Microarray or RNA Sequencing?",
                choices = list("Microarray"),
                selected = "Microarray"
              )
            })
          } else if (input$dataSetType == "Single") {
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
        # Blank Widgets
        output$output6 <- renderUI({})
      }


      # Add or remove CPM radio button
      observeEvent(input$typeOfData, {
        # Define error checks
        errorChecks$continueWorkflow <- TRUE
        errorChecks$geoAccessionCode <- TRUE
        errorChecks$geoMicroarrayAccessionCode <- TRUE
        errorChecks$geoPlatform <- TRUE
        errorChecks$expressionData <- TRUE
        errorChecks$dataInput <- TRUE
        errorChecks$knnDataInput <- TRUE
        errorChecks$pcaPrcompDataInput <- TRUE
        errorChecks$expressionDataOverTwoColumns <- TRUE
        errorChecks$expressionDataOverOneColumns <- TRUE
        errorChecks$differentialGeneExpression <- TRUE
        errorChecks$differentialGeneExpressionGroup <- TRUE
        errorChecks$uploadFile <- TRUE
        errorChecks$uploadFileExtension <- TRUE
        errorChecks$uploadLogData <- TRUE

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
          # Add CPM tool tip
          bsTooltip(
            id = "cpmTransformation",
            title = "This is recommended for raw RNA sequence data.",
            placement = "top",
            trigger = "hover"
          )
        } else if (input$typeOfData == "Microarray") {
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

          # Add knn tool tip
          addTooltip(
            session,
            id = "knnTransformation",
            title = "Rows with over 50% missing values are imputed
              using the overall
              mean per sample. Columns with over
              80% will cause an error in the KNN
              computation.",
            placement = "top",
            trigger = "hover"
          )
        }
      })
    })


    observeEvent(input$dataSetType,{
      if (input$dataSetType == "Combine"){
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
        observeEvent(input$dataSource2,{
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
              textInput("geoAccessionCode2", "GEO accession code", "")
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
            # Process second GEO accession code
            observeEvent(input$geoAccessionCode2, {
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

              # Get the GEO data for all platforms
              all$allGset2 <- reactive({
                tryCatch({
                  # Error handling to ensure geoAccessionCode is populated
                  req(input$geoAccessionCode2)
                  # Notify the user the GEO accession code
                  # is not a GEO series accession code
                  if (substr(input$geoAccessionCode2, 1, 3) != "GSE")
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

              if (errorChecks$continueWorkflow2 == TRUE) {
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
          } else if (input$dataSource2 == "Upload")
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
            # Blank widgets
            output$output10 <- renderUI({})
            output$output11 <- renderUI({})
          }
        })
      } else if (input$dataSetType == "Single"){
        # Set all UI widgets to blank
        output$output2 <- renderUI({})
        output$output7 <- renderUI({})
        output$output8 <- renderUI({})
        output$output9 <- renderUI({})
        output$output10 <- renderUI({})
        output$output11 <- renderUI({})
        output$output14 <- renderUI({})
      }
    })

    # Exploratory data analysis visualisation
    observeEvent(input$exploratoryDataAnalysisButton, {
      # Clear unused memory
      gc()

      # Set all outputs to blank, this resets
      # all the visualizations to blank after clicking analyse
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
      output$iHeatmap <- renderPlotly({

      })

      # Extract information from GSET including expression data
      if (errorChecks$continueWorkflow == TRUE) {
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
          if (is.null(all$gsetData) == TRUE) {
            # Update Error Checks
            errorChecks$geoPlatform <- FALSE
            errorChecks$continueWorkflow <- FALSE

            # Show error
            showNotification("Please select a platform.",
                             type = "error")
          } else if (is.null(all$gsetData) == FALSE) {
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

        } else if
        (input$dataSource == "Upload") {
          # Error handling to prevent non-csvs being uploaded
          if (file_ext(input$file1$name) %in% c('text/csv',
                                                'text/comma-separated-values',
                                                'text/plain',
                                                'csv')) {
            # Update error checks
            errorChecks$uploadFile <- TRUE
            errorChecks$continueWorkflow <- TRUE
            # Ensure a file has been uploaded
            req(input$file1)
            # Extract Expression Data from CSV
            all$expressionData <- tryCatch({
              readCsvFile(input$file1$datapath)
            },
            error = function(e) {
              # return null if there is an error
              return(NULL)
            })


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
            if (is.null(all$expressionData) == TRUE) {
              # Update error checks
              errorChecks$expressionData <- FALSE
              errorChecks$continueWorkflow <- FALSE
            } else if (is.null(all$expressionData) == FALSE) {
              # Update error checks
              errorChecks$expressionData <- TRUE
              errorChecks$continueWorkflow <- TRUE
              # Extract Column Information
              all$columnInfo <- convertExpressionDataToExperimentInformation(
                all$expressionData)
            }
          } else {
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

      # Combining datasets workflow
      if (input$dataSetType == "Combine") {
        if (input$dataSource2 == "GEO"){
          if (errorChecks$continueWorkflow == TRUE &
              errorChecks$continueWorkflow2 == TRUE) {

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
            if (is.null(all$gsetData2) == TRUE) {
              # Update Error Checks
              errorChecks$geoPlatform2 <- FALSE
              errorChecks$continueWorkflow2 <- FALSE

              # Show error
              showNotification("Please select a platform.",
                               type = "error")
            } else if (is.null(all$gsetData2) == FALSE) {
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
                  all$convertedExperimentInformation2)
            }
          }
        } else if
        (input$dataSource2 == "Upload") {
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
            if (is.null(all$expressionData2) == TRUE) {
              # Update error checks
              errorChecks$expressionData2 <- FALSE
              errorChecks$continueWorkflow2 <- FALSE
            } else if (is.null(all$expressionData2) == FALSE) {
              # Update error checks
              errorChecks$expressionData2 <- TRUE
              errorChecks$continueWorkflow2 <- TRUE
              # Extract Column Information
              all$columnInfo2 <- convertExpressionDataToExperimentInformation(
                all$expressionData2)
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

        if (is.null(combinedExpressionData) == TRUE) {
          # Show error
          showNotification("The two gene expression datasets
                                 could not be merged. Please make sure
                                 they have the same platform.  Only the first
                                 gene expression datasets was
                                 processed as a result.",
                           type = "warning")
        } else if (is.null(combinedExpressionData) == FALSE)
        {
          # Perform batch correction
          combinedExpressionDataBatchRemoved <- tryCatch({
            calculateBatchCorrection(
              all$expressionData,
              all$expressionData2,
              combinedExpressionData,
              input$batchCorrection
            )}, error = function(err) {
              # Return null if there is a error in the getGeoObject function
              return(NULL)
            })

          if (is.null(combinedExpressionDataBatchRemoved) == TRUE) {
            # Show error
            showNotification("There was an error performing batch
                                   correction. Therefore the non-batch
                                   corrected data was used.",
                             type = "warning")

            # Update expression data with non-batch corrected data
            all$expressionData <- combinedExpressionData

          } else if (is.null(combinedExpressionDataBatchRemoved)
                     == FALSE)
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
      if (errorChecks$continueWorkflow == TRUE) {
        # Error handling to detect wrong format expression data
        if (isNumeric(all$expressionData) == FALSE) {
          errorChecks$expressionData <- FALSE
          errorChecks$continueWorkflow <- FALSE
          # Display error message
          showNotification(
            "The gene expression data has non-numerical values.
              Please ensure the gene expression data has only numerical values.
              ",
            type = "error"
          )
        } else if (length(all$expressionData) == 0) {
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
                 ((length(all$expressionData) == 0) == FALSE) == TRUE) {
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
          } else if (ncol(all$expressionData) <= 2) {
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

              if (is.null(all$cpm) == TRUE) {
                # Update cpm
                all$cpm <- all$expressionData

                showNotification(
                  "There was an error calculating CPM. Therefore, the
                    original expression data will be used.",
                  type = "warning"
                )
              }
            } else if (input$typeOfData == "Microarray") {
              all$cpm <- all$expressionData
            }
          } else if (input$dataSource == "GEO") {
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
          if (is.null(all$dataInput) == TRUE) {
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
          if (is.null(all$knnDataInput) == TRUE) {
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

          # KNN Column Set Plot
          all$knnColumns <-
            extractSampleNames(all$knnDataInput)

          # Update col info
          all$columnInfo <-
            all$columnInfo[all$knnColumns, ]

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
          if (is.null(pcaPrcompDataInput) == TRUE) {
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
      if (errorChecks$continueWorkflow == TRUE) {
        if (input$dataSource == "GEO") {
          # Experimental Information Display
          output$experimentInfo <- tryCatch({
            renderUI({
              all$convertedExperimentInformation            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })
        } else if (input$dataSource == "Upload") {
          # Update Experimental information
          output$experimentInfo <-  renderUI({
            HTML(
              "<b>Experimental
          Information is not available when processing
          user-uploaded files!</b>"
            )
          })
        }

        # Column Set Plot
        output$columnTable <-
          tryCatch({
            renderDataTable({
              all$columnInfo
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
              all$knnDataInput
            })
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
          renderDataTable(
            all$columnInfo,
            selection = 'multiple',
            server = FALSE
          )
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })

        observeEvent(input$knnColumnTableOne_rows_selected, {
          all$knnColumnTableTwo <- all$columnInfo[
            -input$knnColumnTableOne_rows_selected, ]

          output$knnColumnTableTwo <- tryCatch({
            renderDataTable(
              all$knnColumnTableTwo,
              selection = 'multiple',
              server = FALSE
            )
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })
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
              interactiveBoxAndWhiskerPlot(all$knnDataInput)
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
        # Error handling to prevent errors caused by
        # expression datasets with only one column
        if (errorChecks$expressionDataOverOneColumns == TRUE) {
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
        }
        # Error handling to display a notification
        # if there was an error in PCA
        if (errorChecks$pcaPrcompDataInput  == TRUE) {
          # Interactive PCA Scree Plot
          output$interactivePcaScreePlot <- tryCatch({
            renderPlotly({
              interactivePrcompPcaScreePlot(pcaPrcompDataInput)
            })
          }, error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })

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
          if (errorChecks$expressionDataOverTwoColumns == TRUE) {
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

      if (errorChecks$continueWorkflow == TRUE) {
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
      output$iHeatmap <- renderPlotly({

      })

      if (errorChecks$continueWorkflow == TRUE)
      {
        # Differential gene expression analysis
        gsms <- tryCatch({
          calculateEachGroupsSamplesGsms(
            all$columnInfo,
            row.names(
              all$columnInfo[
                input$knnColumnTableOne_rows_selected,
              ]),
            row.names(
              all$knnColumnTableTwo[
                input$knnColumnTableTwo_rows_selected,
              ])
          )

        }, error = function(cond) {
          return(NULL)
        })
        # Error handling to prevent differential gene expression
        # analysis being performed before exploratory data analysis
        if (is.null(gsms) == TRUE) {
          # Update error check
          errorChecks$continueWorkflow <- FALSE
          showNotification(
            "There was an error running differential gene expression
                  analysis. Please ensure you have performed exploratory data
                  analysis first.",
            type = "error"
          )
        } else if (is.null(gsms) == FALSE) {
          # Update error check
          errorChecks$continueWorkflow <- TRUE
        }

        if (errorChecks$continueWorkflow == TRUE) {
          # Error handling to ensure at least one
          # group has two samples and the other group
          # has at least one sample
          if ((lengths(regmatches(gsms, gregexpr("0", gsms))) > 0 &
               lengths(regmatches(gsms, gregexpr("1", gsms))) > 1) |
              (lengths(regmatches(gsms, gregexpr("0", gsms))) > 1 &
               lengths(regmatches(gsms, gregexpr("1", gsms))) > 0)) {
            if (input$dataSource == "GEO") {
              if (input$dataSetType == "Single") {
                results <- tryCatch({
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
              } else if (input$dataSetType == "Combine") {
                results <- tryCatch({
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
            } else if (input$dataSource == "Upload") {
              results <- tryCatch({
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
              if (is.null(results)) {
                if (input$typeOfData == "RNA Sequencing") {
                  # Try again with non-log data
                  knnDataInput <-
                    calculateKnnImpute(all$cpm, input$knnTransformation)
                  results <- tryCatch({
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
            if (is.null(results) == TRUE) {
              # Update Error Check
              errorChecks$continueWorkflow <- FALSE
              # Display notification
              showNotification(
                "There was an error calculating the
                 differential gene expression analysis!",
                type = "error"
              )
            } else if (is.null(results) == FALSE) {
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

      if (errorChecks$continueWorkflow == TRUE) {
        adjustment <- convertAdjustment(input$pValueAdjustment)
        tT <-
          calculateTopDifferentiallyExpressedGenes(results$fit2,
                                                   adjustment)

        dT <-
          calculateDifferentialGeneExpressionSummary(
            results$fit2,
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
            interactiveHistogramPlot(results$fit2, adjustment)
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
              interactiveQQPlot(results$fit2, dT, ct)
            })
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          })


        # Interactive Volcano Plot
        output$iDEVolcano <- tryCatch({
          renderPlotly({
            interactiveVolcanoPlot(results$fit2, dT, ct)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })


        # Interactive Mean Difference Plot
        output$iDEMd <- tryCatch({
          renderPlotly({
            interactiveMeanDifferencePlot(results$fit2, dT, ct)
          })
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        })

        # Update Interactive Heatmap Plot Gene number max
        updateNumericInput(session, "numberOfGenes", max = nrow(tT))

        # Interactive Heatmap Plot
        output$iHeatmap <-
          tryCatch({
            renderPlotly({
              interactiveDGEHeatMapPlot(results$ex,
                                        input$limmaPrecisionWeights,
                                        input$numberOfGenes,
                                        tT)
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

    })

  })
  return(datasetInformationServer)
}
