library(plotly)
library(ggplot2)
library(limma)

interactiveBoxAndWhiskerPlot <- function(data, geoAccessionCode, platform) {
  data <- as.data.frame(data)
  fig <- plot_ly(type = "box", quartilemethod="linear")
  i = 1
  for(col in names(data)) {
    fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
    i <- i+1
  }
  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode, "/"), platform)))
  fig
}

interactiveDesnityPlot <- function(data, geoAccessionCode, platform) {
  data <- as.data.frame(data)
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
  fig
}

interactiveThreeDDesnityPlot <- function(data, geoAccessionCode, platform) {
  data <- as.data.frame(data)
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

# Look into consolidating data <- na.omit(data) function
interactiveUmapPlot <- function(data) {
  data <- data[!duplicated(data), ]  # remove duplicates
  nNeighbors <- 5
  ump <- umap(t(ex), n_neighbors = nNeighbors, random_state = 123)
  i <- 1
  fig <- plot_ly(type = 'scatter', mode = 'markers')
  for(row in row.names(ump$layout)){
    fig <- fig %>% add_trace(x = ump$layout[i,][1], y = ump$layout[i,][2], name = row)
    i <- i+1
  }
  fig <- fig %>% layout(
    title = (paste('UMAP plot, number of nearest neighbors used =',nNeighbors)))
  fig
}

interactiveMeanVariancePlot <- function(data, geoAccessionCode) {
  data <- lmFit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(data = data, x = ~Amean, y = ~sigma, type = 'scatter', mode = 'markers', marker = list(
    color = 'rgb(17, 157, 255)',
    size = 3,
    line = list(
      color = 'rgb(0, 0, 0)',
      width = 1
    )))
  fig <- fig %>% layout(
    title = (paste('Mean variance trend, ',geoAccessionCode)))
  fig
  }
  