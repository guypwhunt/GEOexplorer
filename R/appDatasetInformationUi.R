datasetInformationUi <- tabPanel("Dataset Information",
         tabsetPanel(type = "tabs",
                     tabPanel("Experiment Information", br(), htmlOutput('experimentInfo')),
                     tabPanel("Column Details", dataTableOutput('columnTable')),
                     tabPanel("Dataset", dataTableOutput('table'))
         ))
