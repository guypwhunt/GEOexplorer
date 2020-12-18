# Load Packages
library(shiny)
library(maps)
library(mapproj)

# Load Data
expressiondata <- read.csv(file = 'analysis-output.csv')

print(expressiondata)

# Source Helper Functions
source("helpers.R")


# Define UI ----
ui <- fluidPage(
  titlePanel("GEO2R Visulisation"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Create demographic maps with information from the 2010 US Census"),
      
      selectInput("var", 
                  strong("Choose a variable to display"), 
                  choices = list("Percentage White", "Percentage Black","Percentage Hispanic","Percentage Asian"), 
                  selected = "Percentage White"),
      sliderInput("range", 
                  strong("Range of Interest"),
                  min = 0, max = 100, value = c(0,100))),
    mainPanel(plotOutput("map")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  output$map <- renderPlot({
    data <- switch(input$var, 
                   "Percentage White" = counties$white,
                   "Percentage Black" = counties$black,
                   "Percentage Hispanic" = counties$hispanic,
                   "Percentage Asian" = counties$asian)
    
    color <- switch(input$var, 
                    "Percentage White" = "darkgreen",
                    "Percentage Black" = "blue",
                    "Percentage Hispanic" = "red",
                    "Percentage Asian" = "orange")
    
    legend <- switch(input$var, 
                     "Percentage White" = "% White",
                     "Percentage Black" = "% Black",
                     "Percentage Hispanic" = "% Hispanic",
                     "Percentage Asian" = "% Asian")
    
    percent_map(data, color, legend, input$range[1], input$range[2])
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)