library(plotly)
library(ggplot2)

interactiveBoxAndWhiskerPlot <- function(data, geoAccessionCode, platform) {
  data <- na.omit(data)
  data <- as.data.frame(data)
  fig <- plot_ly(type = "box", quartilemethod="linear")
  i = 1
  for(col in names(data)) {
    fig <- fig %>% add_trace(x = names(data)[i], y = data[,i], quartilemethod="linear", name=names(data)[i])
    i <- i+1
  }
  fig <- fig %>% layout(title = (paste(paste(geoAccessionCode, "/"), platform)))
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
}

# Look into consolidating data <- na.omit(data) function
interactiveUmapPlot <- function(data) {
  data <- na.omit(data)
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
    
}