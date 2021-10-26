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
sourceServer2 <- function(input, output, session) {
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

    # Add Significance Level Cutoff Tool Tips
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

    # Define Variables
    gsetData <- NULL
    geoAccessionCode <- ""
    all <- reactiveValues()
    ct <- 1

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
    output$table <- renderDataTable(expressionData())


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
    output$columnTable <- renderDataTable({
      columnInfo()
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

      # Get a list of all the columns
      columns <- extractSampleNames(expressionData())

      if (input$typeOfData == "RNA Sequencing")
      {
        # Raw counts are converted to counts-per-million (CPM)
        all$cpm <-
          calculateCountsPerMillion(expressionData(), input$cpmTransformation)
      } else if (input$typeOfData == "Microarray")
      {
        all$cpm < -expressionData()
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
      output$logTransformationText <- renderUI({
        helpText(autoLogInformation)
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

      output$knnColumnTable <- renderDataTable(
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

      # Expression dataset table
      output$table <- renderDataTable({
        all$knnDataInput
      })

      # Interactive Box-and-Whisker Plot
      output$interactiveBoxAndWhiskerPlot <- renderPlotly({
        interactiveBoxAndWhiskerPlot(all$knnDataInput,
                                     geoAccessionCode,
                                     input$platform)
      })

      # Interactive Density Plot
      output$interactiveDensityPlot <- renderPlotly({
        interactiveDensityPlot(naOmitInput,
                               geoAccessionCode,
                               input$platform)
      })

      # 3D Interactive Density Plot
      output$interactiveThreeDDensityPlot <- renderPlotly({
        interactiveThreeDDensityPlot(naOmitInput,
                                     geoAccessionCode,
                                     input$platform)
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
        output$interactiveUmapPlot <- renderPlotly({
          interactiveUmapPlot(naOmitInput,
                              input$knn,
                              geoAccessionCode)
        })

        # Heatmap Plot
        output$interactiveHeatMapPlot <- renderPlotly({
          interactiveHeatMapPlot(naOmitInput)
        })

        # Interactive Mean Variance Plot
        output$interactiveMeanVariancePlot <- renderPlotly({
          interactiveMeanVariancePlot(naOmitInput,
                                      geoAccessionCode,
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
                                          geoAccessionCode)
          })

          # Interactive PCA Individual Plot
          output$interactivePcaIndividualsPlot <-
            renderPlotly({
              interactivePrcompPcaIndividualsPlot(pcaPrcompDataInput,
                                                  geoAccessionCode,
                                                  gsetData)
            })

          # Interactive PCA Variables Plot
          output$interactivePcaVariablesPlot <-
            renderPlotly({
              interactivePrcompPcaVariablesPlot(pcaPrcompDataInput,
                                                geoAccessionCode)
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
        calculateEachGroupsSamplesFromDataFrame(
          as.data.frame(
            sapply(
              seq_len(
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

          fit2 <-
            tryCatch({
            calculateDifferentialGeneExpressionRnaSeq(
              all$knnDataInput, gsms,
              input$limmaPrecisionWeights,
              input$forceNormalization)
          }
          , error = function(cond) {
            return(NULL)
          })
          if (is.null(fit2)) {
            # Try again with non-log data
            knnDataInput <- calculateKnnImpute(all$cpm, input$knnTransformation)
            fit22 <- tryCatch({
              calculateDifferentialGeneExpressionRnaSeq(
                knnDataInput, gsms,
                input$limmaPrecisionWeights,
                input$forceNormalization)
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
          if (is.null(fit2) == FALSE) {
            # Convert the UI adjustment into the value needed for the backend
            adjustment <- convertAdjustment(input$pValueAdjustment)

            # Calculate the top differentially expressed genes
            tT <-
              calculateTopDifferentiallyExpressedGenes(fit2,
                                                       adjustment)

            # Calculate genes that are upregulated and downregulated
            dT <-
              calculateDifferentialGeneExpressionSummary(
                fit2,
                adjustment,
                input$significanceLevelCutOff)

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
