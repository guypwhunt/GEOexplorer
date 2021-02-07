library(plotly)

# load the diamonds dataset from the ggplot2 package
data(diamonds, package = "ggplot2")
summary(diamonds)
names(diamonds)

new_diamonds <- diamonds[,1]
new_diamonds["depth"] <- diamonds[,5]
new_diamonds["table"] <- diamonds[,6]
new_diamonds["price"] <- diamonds[,7]
new_diamonds["x"] <- diamonds[,8]
new_diamonds["y"] <- diamonds[,9]
new_diamonds["z"] <- diamonds[,10]

new_diamonds <- as.data.frame(new_diamonds)

summary(new_diamonds)


fig <- plot_ly(type = 'scatter', mode = 'lines', name = 'Density Plot')
i = 1
for(col in names(new_diamonds)) {
  density <- density(new_diamonds[,i])
  print(density)
  fig <- fig %>% add_trace(x = ~density$x, y = ~density$y, name = col, fill = 'tozeroy')
  i <- i+1
}

fig <- fig %>% layout(xaxis = list(title = 'Carat'),
                      yaxis = list(title = 'Density'))

fig