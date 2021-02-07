# 3D Plot Idea

#library(plotly)

# load the diamonds dataset from the ggplot2 package
#data(diamonds, package = "ggplot2")
#summary(diamonds)

#df = as.data.frame(diamonds)
#plot_ly(df, x = ~depth, y = ~carat, z = ~table, group_by = ~cut, type = "scatter3d", mode = "lines") 

library(plotly)

library(ggplot2)
set.seed(1234)

# load the diamonds dataset from the ggplot2 package
data(diamonds, package = "ggplot2")

dfGamma = data.frame(nu75 = rgamma(100, 0.75),
                     nu1 = rgamma(100, 1),
                     nu2 = rgamma(100, 2))

dfGamma = stack(dfGamma)
dfGamma = diamonds
summary(dfGamma)


p <- ggplot(dfGamma, aes(x = depth)) +
  stat_density(aes(group = cut, color = cut),position="identity",geom="line")

fig <- ggplotly(p)

fig