library(plotly)
library(ggplot2)

interactiveBoxAndWhiskerPlot <- function(data) {
  data <- na.omit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(data, type = "box", quartilemethod="linear")
  i = 1
  for(col in names(data)) {
    fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
    i <- i+1
  }
  return(fig)
}

interactiveDesnityPlot <- function(data) {
  data <- as.data.frame(data)
  data <- na.omit(data)
  fig <- plot_ly(type = 'scatter', mode = 'lines', name = 'Density Plot')
  i <- 1
  for(col in names(data)) {
    density <- density(data[,i])
    fig <- fig %>% add_trace(x = density$x, y = density$y, name = col)
    i <- i+1
  }
  
  fig <- fig %>% layout(xaxis = list(title = 'Intensity'),
                        yaxis = list(title = 'Density'))
  
  return(fig)
}

interactiveThreeDDesnityPlot <- function(data) {
  data <- as.data.frame(data)
  data <- na.omit(data)
  fig <- plot_ly(type = 'scatter3d', mode = 'lines', name = 'Density Plot')
  i <- 1
  for(col in names(data)) {
    density <- density(data[,i])
    fig <- fig %>% add_trace(x = density$x, y = i, z = density$y, name = col)
    i <- i+1
  }
  
  fig <- fig %>% layout(xaxis = list(title = 'Intensity'),
                        yaxis = list(title = 'Density'))
  
  return(fig)
}