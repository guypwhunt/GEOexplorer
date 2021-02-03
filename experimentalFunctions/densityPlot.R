library(plotly)

# load the diamonds dataset from the ggplot2 package
data(diamonds, package = "ggplot2")
summary(diamonds)

df = as.data.frame(diamonds)
plot_ly(df, x = ~depth, y = ~carat, z = ~table, group_by = ~cut, type = "scatter3d", mode = "lines") 

help('plotly_data')