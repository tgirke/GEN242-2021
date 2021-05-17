## pre {

##   max-height: 300px;

##   overflow-y: auto;

## }

## 
## pre[class] {

##   max-height: 300px;

## }


## .scroll-300 {

##   max-height: 300px;

##   overflow-y: auto;

##   background-color: inherit;

## }


## ----ggplot_scatter_plot1, eval=TRUE------------------------------------------
library(ggplot2)
set.seed(1410)  
dsmall <- as.data.frame(diamonds[sample(nrow(diamonds), 1000), ])
p <- ggplot(dsmall, aes(carat, price, color=color)) + 
            geom_point(size=4)
print(p) 


## ----ggplot_scatter_plot_interactive1, eval=TRUE, message=FALSE---------------
library(plotly)
ggplotly(p)

