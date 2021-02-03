library(plotly)

interactiveBoxAndWhiskerPlot <- function(data) {
  data <- na.omit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(data, type = "box", quartilemethod="linear")
  i = 1
  for(col in names(data)) {
    fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
    i <- i+1
  }
  }