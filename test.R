library(shiny)
library(DT)

ui <- fluidPage(
  title = 'Selectinput column in a table',
  h3("Source:", tags$a("Yihui Xie", href = "https://yihui.shinyapps.io/DT-radio/")),
  DT::dataTableOutput('foo'),
  verbatimTextOutput('sel')
)

server <- function(input, output, session) {
  data <- head(iris, 5)

  for (i in 1:nrow(data)) {
    data$species_selector[i] <- as.character(selectInput(paste0("sel", i), "", choices = unique(iris$Species), width = "100px"))
  }

  output$foo = DT::renderDataTable(
    data, escape = FALSE, selection = 'none', server = FALSE,
    options = list(dom = 't', paging = FALSE, ordering = FALSE),
    callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());")
  )
  output$sel = renderPrint({
    str(sapply(1:nrow(data), function(i) input[[paste0("sel", i)]]))
  })
}

shinyApp(ui, server)
