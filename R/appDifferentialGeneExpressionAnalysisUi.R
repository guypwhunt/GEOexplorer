differentialGeneExpressionAnalysisUi <- tabPanel("Differential Gene Expression Analysis",
         tabsetPanel(type = "tabs",
                     tabPanel("Set Parameters",
                              mainPanel(
                                DT::dataTableOutput('knnColumnTable'),
                                fluidRow(
                                  column(6,
                                         br(),
                                         uiOutput("dyncolumns"),
                                         selectInput("pValueAdjustment", "Apply adjustment to the P-values:",
                                                     choices = c("Benjamini & Hochberg (False discovery rate)", "Benjamini & Yekutieli", "Bonferroni", "Holm", "None")), # "Hochberg" and "Hommel" were removed
                                         radioButtons("limmaPrecisionWeights",
                                                      label="Apply limma precision weights (vooma):",
                                                      choices=list("Yes","No"),
                                                      selected="No")
                                  ),
                                  br(),
                                  column(6,
                                         radioButtons("forceNormalization",
                                                      label="Force normalization:",
                                                      choices=list("Yes","No"),
                                                      selected="No"),
                                         sliderInput("significanceLevelCutOff", "Significance level cut-off:",
                                                     min = 0, max = 1,
                                                     value = 0.05),
                                         actionButton("differentialExpressionButton", "Analyse")
                                  )
                                )
                              )
                     ),
                     tabPanel("Top Differentially Expressed Genes", br(), span("The table below displays the top differentially expressed genes between the groups selected."), br(), br(), span("adj.P.Val is the P-value after adjustment for multiple testing. This column is generally recommended as the primary statistic by which to interpret results. Genes with the smallest P-values will be the most reliable."), br(), span("P.Value is the Raw P-value"), br(), span("t is the Moderated t-statistic (only available when two groups of Samples are defined)"), br(), span("B is the B-statistic or log-odds that the gene is differentially expressed (only available when two groups of Samples are defined)"), br(), span("logFC is the Log2-fold change between two experimental conditions (only available when two groups of Samples are defined)"), br(), span("F is the moderated F-statistic which combines the t-statistics for all the pair-wise comparisons into an overall test of significance for that gene (only available when more than two groups of Samples are defined)"),  br(), br(), downloadButton("downloadData", "Download"), br(), br(), dataTableOutput('dETable')),
                     tabPanel("Histogram Plot", br(), span("Generated using hist. Use to view the distribution of the P-values in the analysis results. The P-value here is the same as in the Top differentially expressed genes table and computed using all selected contrasts. While the displayed table is limited by size this plot allows you to see the 'big picture' by showing the P-value distribution for all analyzed genes."), br(), br(), plotlyOutput('iDEHistogram')),
                     tabPanel("Venn Diagram Plot", br(), span("Generated using limma (vennDiagram). Use to explore the overlap in significant genes between multiple contrasts.") , plotOutput('dEVennDiagram')),
                     tabPanel("Q-Q Plot", br(), span("Generated using limma (qqt) and R plotly. Plots the quantiles of a data sample against the theoretical quantiles of a Student's t distribution. This plot helps to assess the quality of the limma test results. Ideally the points should lie along a straight line, meaning that the values for moderated t-statistic computed during the test follow their theoretically predicted distribution."), br(), br(), plotlyOutput('iDEQQ')),
                     tabPanel("Volcano Plot", br(), span("Generated using R plotly. A volcano plot displays statistical significance (-log10 P value) versus magnitude of change (log2 fold change) and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff value") , br(), br(), plotlyOutput('iDEVolcano')),
                     tabPanel("Mean Difference Plot", br(), span("Generated using R plotly. A mean difference (MD) plot displays log2 fold change versus average log2 expression values and is useful for visualizing differentially expressed genes. Highlighted genes are significantly differentially expressed at the selected adjusted p-value cutoff."), br(), br(), plotlyOutput('iDEMd'))

         )

)
