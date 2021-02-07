library(plotly)
library(ggplot2)

interactiveBoxAndWhiskerPlot <- function(data, geoAccessionCode, platform) {
  data <- na.omit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(data, type = "box", quartilemethod="linear")
  i = 1
  for(col in names(data)) {
    fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
    i <- i+1
  }
  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode, "/"), platform)))
  return(fig)
}

interactiveDesnityPlot <- function(data, geoAccessionCode, platform) {
  data <- as.data.frame(data)
  data <- na.omit(data)
  fig <- plot_ly(type = 'scatter', mode = 'lines', name = (paste(paste(geoAccessionCode,platform),'value distribution')))
  i <- 1
  for(col in names(data)) {
    density <- density(data[,i])
    fig <- fig %>% add_trace(x = density$x, y = density$y, name = col)
    i <- i+1
  }
  
  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode,platform),'value distribution')),
                        xaxis = list(title = 'Intensity'),
                        yaxis = list(title = 'Density'))
  return(fig)
}

interactiveThreeDDesnityPlot <- function(data, geoAccessionCode, platform) {
  data <- as.data.frame(data)
  data <- na.omit(data)
  fig <- plot_ly(type = 'scatter3d', mode = 'lines', name = (paste(paste(geoAccessionCode,platform),'value distribution')))
  i <- 1
  for(col in names(data)) {
    density <- density(data[,i])
    fig <- fig %>% add_trace(x = density$x, y = i, z = density$y, name = col)
    i <- i+1
  }
  
  fig <- fig %>% layout(
    title = (paste(paste(geoAccessionCode,platform),'value distribution')),
    scene = list(
      xaxis = list(title = "Intensity"),
      yaxis = list(title = ""),
      zaxis = list(title = "Density")
    ))
  fig
}