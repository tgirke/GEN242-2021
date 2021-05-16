---
title: Graphics and Data Visualization in R
author: "First/last name (first.last@ucr.edu)"
date: "Last update: 16 May, 2021" 
output:
  html_document:
    toc: true
    toc_float:
        collapsed: true
        smooth_scroll: true
    toc_depth: 3
    fig_caption: yes
    code_folding: show
    number_sections: true

fontsize: 14pt
bibliography: bibtex.bib
weight: 14
type: docs
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('rgraphics.Rmd', c('html_document'), clean=F); knitr::knit('rgraphics.Rmd', tangle=TRUE)"
-->


<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  document.querySelector("h1").className = "title";
});
</script>
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
  var links = document.links;  
  for (var i = 0, linksLength = links.length; i < linksLength; i++)
    if (links[i].hostname != window.location.hostname)
      links[i].target = '_blank';
});
</script>

<div style="text-align: right"> 
Source code downloads: &nbsp; &nbsp;
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rgraphics/rgraphics.Rmd) ] &nbsp; &nbsp; 
[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rgraphics/rgraphics.R) ]
</div>


## Overview

### Graphics in R

-   Powerful environment for visualizing scientific data
-   Integrated graphics and statistics infrastructure
-   Publication quality graphics
-   Fully programmable
-   Highly reproducible
-   Full [LaTeX](http://www.latex-project.org/), [Sweave](http://www.stat.auckland.ac.nz/~dscott/782/Sweave-manual-20060104.pdf), [knitr](http://yihui.name/knitr/) and [R Markdown](http://rmarkdown.rstudio.com/) support.
-   Vast number of R packages with graphics utilities


### Documentation on Graphics in R

- General 
    - [Graphics Task Page](http://cran.r-project.org/web/views/Graphics.html)
    - [R Graph Gallery](http://www.r-graph-gallery.com/)
    - [R Graphical Manual](https://www.imsbio.co.jp/RGM/R_image_list?page=575&init=true)
    - [Paul Murrell’s book R (Grid) Graphics](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)

- Interactive graphics
    - [rggobi (GGobi)](http://www.ggobi.org/)
    - [iplots](http://www.iplots.org/)
    - [Open GL (`rgl`)](http://rgl.neoscientists.org/gallery.shtml)
    - [plotly](https://plotly.com/r/)


### Graphics Environments

- Viewing and savings graphics in R
    - On-screen graphics
    - postscript, pdf, svg
    - jpeg/png/wmf/tiff/...

- Four major graphics environments
    - Low-level infrastructure
        - R Base Graphics (low- and high-level)
        - `grid`: [Manual](http://www.stat.auckland.ac.nz/~paul/grid/grid.html), [Book](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)
    - High-level infrastructure
        - `lattice`: [Manual](http://lmdvr.r-forge.r-project.org), [Intro](http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm), [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
        - `ggplot2`: [Manual](http://docs.ggplot2.org/current/), [Intro](http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html), [Book](http://had.co.nz/ggplot2/book/)


## Base Graphics

### Overview

- Important high-level plotting functions
    - `plot`: generic x-y plotting
    - `barplot`: bar plots
    - `boxplot`: box-and-whisker plot
    - `hist`: histograms
    - `pie`: pie charts
    - `dotchart`: cleveland dot plots
    - `image, heatmap, contour, persp`: functions to generate image-like plots
    - `qqnorm, qqline, qqplot`: distribution comparison plots
    - `pairs, coplot`: display of multivariant data

- Help on these functions
    - `?myfct`
    - `?plot`
    - `?par`

### Preferred Input Data Objects

- Matrices and data frames
- Vectors
- Named vectors

### Scatter Plots

#### Basic scatter plots

Sample data set for subsequent plots

```r
set.seed(1410)
y <- matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))
y
```

```
##            A          B            C
## a 0.26904539 0.47439030 0.4427788756
## b 0.53178658 0.31128960 0.3233293493
## c 0.93379571 0.04576263 0.0004628517
## d 0.14314802 0.12066723 0.4104402000
## e 0.57627063 0.83251909 0.9884746270
## f 0.49001235 0.38298651 0.8235850153
## g 0.66562596 0.70857731 0.7490944304
## h 0.50089252 0.24772695 0.2117313873
## i 0.57033245 0.06044799 0.8776291364
## j 0.04087422 0.85814118 0.1061618729
```

```r
plot(y[,1], y[,2]) 
```

<img src="rgraphics_files/figure-html/scatter_plot_basic-1.png" width="672" />

#### All pairs


```r
pairs(y) 
```

<img src="rgraphics_files/figure-html/scatter_plot_allpairs-1.png" width="672" />

#### Plot labels


```r
plot(y[,1], y[,2], pch=20, col="red", main="Symbols and Labels")
text(y[,1]+0.03, y[,2], rownames(y))
```

<img src="rgraphics_files/figure-html/scatter_plot_labels-1.png" width="672" />

#### More examples

Print instead of symbols the row names


```r
plot(y[,1], y[,2], type="n", main="Plot of Labels")
text(y[,1], y[,2], rownames(y)) 
```

Usage of important plotting parameters


```r
grid(5, 5, lwd = 2) 
op <- par(mar=c(8,8,8,8), bg="lightblue")
plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2, 
     cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label", 
     ylab="y label", main="My Main", sub="My Sub")
par(op)
```

Important arguments}
    - `mar`: specifies the margin sizes around the plotting area in order: `c(bottom, left, top, right)` 
    - `col`: color of symbols
    - `pch`: type of symbols, samples: `example(points)`
    - `lwd`: size of symbols
    - `cex.*`: control font sizes
    - For details see `?par`


Add a regression line to a plot


```r
plot(y[,1], y[,2])
myline <- lm(y[,2]~y[,1]); abline(myline, lwd=2) 
```

<img src="rgraphics_files/figure-html/scatter_plot_regress-1.png" width="672" />

```r
summary(myline) 
```

```
## 
## Call:
## lm(formula = y[, 2] ~ y[, 1])
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.40357 -0.17912 -0.04299  0.22147  0.46623 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   0.5764     0.2110   2.732   0.0258 *
## y[, 1]       -0.3647     0.3959  -0.921   0.3839  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3095 on 8 degrees of freedom
## Multiple R-squared:  0.09589,	Adjusted R-squared:  -0.01712 
## F-statistic: 0.8485 on 1 and 8 DF,  p-value: 0.3839
```

Same plot as above, but on log scale


```r
plot(y[,1], y[,2], log="xy") 
```

<img src="rgraphics_files/figure-html/scatter_plot_log_scale-1.png" width="672" />

Add a mathematical expression to a plot


```r
plot(y[,1], y[,2]); text(y[1,1], y[1,2], 
     expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) 
```

<img src="rgraphics_files/figure-html/scatter_plot_math-1.png" width="672" />

#### Exercise 1

- __Task 1__: Generate scatter plot for first two columns in `iris` data frame and color dots by its `Species` column.
- __Task 2__: Use the `xlim/ylim` arguments to set limits on the x- and y-axes so that all data points are restricted to the left bottom quadrant of the plot. 

Structure of iris data set:

```r
class(iris)
```

```
## [1] "data.frame"
```

```r
iris[1:4,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
```

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```



### Line Plots

#### Single Data Set


```r
plot(y[,1], type="l", lwd=2, col="blue") 
```

<img src="rgraphics_files/figure-html/line_plot_single-1.png" width="672" />

#### Many Data Sets

Plots line graph for all columns in data frame `y`. The `split.screen` function is used in this example in a for loop to overlay several line graphs in the same plot. 


```r
split.screen(c(1,1)) 
```

```
## [1] 1
```

```r
plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1)
for(i in 2:length(y[1,])) { 
	screen(1, new=FALSE)
	plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", 
             xlab="", main="", bty="n") 
}
```

<img src="rgraphics_files/figure-html/line_plot_many-1.png" width="672" />

```r
close.screen(all=TRUE) 
```

### Bar Plots

#### Basics


```r
barplot(y[1:4,], ylim=c(0, max(y[1:4,])+0.3), beside=TRUE, 
        legend=letters[1:4]) 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1)
     +sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.04) 
```

<img src="rgraphics_files/figure-html/bar_plot_basic-1.png" width="672" />
	
#### Error bars


```r
bar <- barplot(m <- rowMeans(y) * 10, ylim=c(0, 10))
stdev <- sd(t(y))
arrows(bar, m, bar, m + stdev, length=0.15, angle = 90)
```

<img src="rgraphics_files/figure-html/bar_plot_error_bar-1.png" width="672" />

#### Mirrored bar plot


```r
df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))
plot(c(0,12), range(df$y), type = "n")
barplot(height = df$y[df$group == "Above"], add = TRUE, axes = FALSE)
barplot(height = df$y[df$group == "Below"], add = TRUE, axes = FALSE)
```

<img src="rgraphics_files/figure-html/bar_plot_mirrored-1.png" width="672" />

### Histograms


```r
hist(y, freq=TRUE, breaks=10)
```

<img src="rgraphics_files/figure-html/hist_plot-1.png" width="672" />

### Density Plots}


```r
plot(density(y), col="red")
```

<img src="rgraphics_files/figure-html/density_plot-1.png" width="672" />

### Pie Charts


```r
pie(y[,1], col=rainbow(length(y[,1]), start=0.1, end=0.8), clockwise=TRUE)
legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, 
col=rainbow(length(y[,1]), start=0.1, end=0.8), ncol=1) 
```

<img src="rgraphics_files/figure-html/pie_chart-1.png" width="672" />

### Color Selection Utilities

Default color palette and how to change it


```r
palette()
```

```
## [1] "black"   "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "gray62"
```

```r
palette(rainbow(5, start=0.1, end=0.2))
palette()
```

```
## [1] "#FF9900" "#FFBF00" "#FFE600" "#F2FF00" "#CCFF00"
```

```r
palette("default")
```

The `gray` function allows to select any type of gray shades by providing values from 0 to 1


```r
gray(seq(0.1, 1, by= 0.2))
```

```
## [1] "#1A1A1A" "#4D4D4D" "#808080" "#B3B3B3" "#E6E6E6"
```

Color gradients with `colorpanel` function from `gplots` library


```r
library(gplots)
colorpanel(5, "darkblue", "yellow", "white")
```

Much more on colors in R see Earl Glynn's [color chart](http://research.stowers-institute.org/efg/R/Color/Chart/)


### Arranging Several Plots on Single Page

With `par(mfrow=c(nrow, ncol))` one can define how several plots are arranged next to each other.


```r
par(mfrow=c(2,3)) 
for(i in 1:6) plot(1:10) 
```

<img src="rgraphics_files/figure-html/par_mfrow-1.png" width="672" />

#### Arranging Plots with Variable Width

The `layout` function allows to divide the plotting device into variable numbers of rows and columns with the column-widths and the row-heights specified in the respective arguments.


```r
nf <- layout(matrix(c(1,2,3,3), 2, 2, byrow=TRUE), c(3,7), c(5,5), 
             respect=TRUE)
# layout.show(nf)
for(i in 1:3) barplot(1:10) 
```

<img src="rgraphics_files/figure-html/layout_plot-1.png" width="672" />

### Saving Graphics to Files

After the `pdf()` command all graphs are redirected to file `test.pdf`. Works for all common formats similarly: jpeg, png, ps, tiff, ...


```r
pdf("test.pdf"); plot(1:10, 1:10); dev.off() 
```

Generates Scalable Vector Graphics (SVG) files that can be edited in vector graphics programs, such as InkScape.


```r
svg("test.svg"); plot(1:10, 1:10); dev.off() 
```

#### Exercise 2

Bar plots
        
- __Task 1__: Calculate the mean values for the `Species` components of the first four columns in the `iris` data set. Organize the results in a matrix where the row names are the unique values from the `iris Species` column and the column names are the same as in the first four `iris` columns. 
- __Task 2__: Generate two bar plots: one with stacked bars and one with horizontally arranged bars. 

Structure of iris data set:


```r
class(iris)
```

```
## [1] "data.frame"
```

```r
iris[1:4,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
```

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```



## Grid Graphics

- What is `grid`?
    - Low-level graphics system 
    - Highly flexible and controllable system
    - Does not provide high-level functions 
    - Intended as development environment for custom plotting functions 
    - Pre-installed on new R distributions

- Documentation and Help
    - [Manual](http://www.stat.auckland.ac.nz/~paul/grid/grid.html)
    - [Book](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)

## lattice Graphics

- What is `lattice`?
    - High-level graphics system 
    - Developed by Deepayan Sarkar 
    - Implements Trellis graphics system from S-Plus
    - Simplifies high-level plotting tasks: arranging complex graphical features 
    - Syntax similar to R's base graphics

- Documentation and Help
    - [Manual](http://lmdvr.r-forge.r-project.org)
    - [Intro](https://www.statmethods.net/advgraphs/trellis.html)
    - [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
		
Open a list of all functions available in the lattice package


```r
library(lattice) 
library(help=lattice) 
```

Accessing and changing global parameters:


```r
?lattice.options
?trellis.device
```

### Scatter Plot Sample


```r
library(lattice)
p1 <- xyplot(1:8 ~ 1:8 | rep(LETTERS[1:4], each=2), as.table=TRUE) 
plot(p1)
```

<img src="rgraphics_files/figure-html/scatter_plot_lattice-1.png" width="672" />

### Line Plot Sample


```r
library(lattice)
p2 <- parallelplot(~iris[1:4] | Species, iris, horizontal.axis = FALSE, 
              layout = c(1, 3, 1))  
plot(p2)
```

<img src="rgraphics_files/figure-html/line_plot_lattice-1.png" width="672" />

## ggplot2 Graphics

- What is `ggplot2`?
    - High-level graphics system developed by Hadley Wickham
    - Implements grammar of graphics from [Leland Wilkinson](http://www.amazon.com/Grammar-Graphics-Leland-Wilkinson/dp/0387987746) 
    - Streamlines many graphics workflows for complex plots
    - Syntax centered around main `ggplot` function 
    - Simpler `qplot` function provides many shortcuts
        
- Documentation and Help
    - [Manual](https://ggplot2.tidyverse.org/reference/)
    - [Book](https://ggplot2-book.org/)
    - [Cookbook for R](http://www.cookbook-r.com/Graphs/)

### `ggplot2` Usage
	
- `ggplot` function accepts two main arguments
    - Data set to be plotted 
	- Aesthetic mappings provided by `aes` function
- Additional parameters such as geometric objects (_e.g._ points, lines, bars) are passed on by appending them with `+` as separator. 
- List of available `geom_*` functions see [here](https://ggplot2.tidyverse.org/reference/) 
- Settings of plotting theme can be accessed with the command `theme_get()` and its settings can be changed with `theme()`. 
- Preferred input data object 
    - `qplot`: `data.frame` or `tibble` (support for `vector`, `matrix`, `...`)
    - `ggplot`: `data.frame` or `tibble`
- Packages with convenience utilities to create expected inputs
    - `dplyr` (`plyr`)
    - `reshape`

### `qplot` Function

The syntax of `qplot` is similar as R's basic `plot` function

- Arguments
    - `x`: x-coordinates (_e.g._ `col1`)
    - `y`: y-coordinates (_e.g._ `col2`)
	- `data`: `data.frame` or `tibble` with corresponding column names
	- `xlim, ylim`: _e.g._ `xlim=c(0,10)` 
    - `log`: _e.g._ `log="x"` or `log="xy"`
	- `main`: main title; see `?plotmath` for mathematical formula
	- `xlab, ylab`: labels for the x- and y-axes
	- `color`, `shape`, `size`
	- `...`: many arguments accepted by `plot` function

### `qplot`: scatter plot basics

Create sample data, here 3 vectors: `x`, `y` and `cat`

```r
library(ggplot2)
x <- sample(1:10, 10); y <- sample(1:10, 10); cat <- rep(c("A", "B"), 5)
```

Simple scatter plot


```r
qplot(x, y, geom="point")
```

<img src="rgraphics_files/figure-html/qplot_scatter_plot-1.png" width="672" />

Prints dots with different sizes and colors


```r
qplot(x, y, geom="point", size=x, color=cat, 
      main="Dot Size and Color Relative to Some Values")
```

<img src="rgraphics_files/figure-html/qplot_scatter_plot_dot_param-1.png" width="672" />

Drops legend


```r
qplot(x, y, geom="point", size=x, color=cat) + 
      theme(legend.position = "none")
```

<img src="rgraphics_files/figure-html/qplot_scatter_plot_no_legend-1.png" width="672" />

Plot different shapes


```r
qplot(x, y, geom="point", size=5, shape=cat)
```

<img src="rgraphics_files/figure-html/qplot_scatter_plot_shapes-1.png" width="672" />

#### Colored groups


```r
p <- qplot(x, y, geom="point", size=x, color=cat, 
            main="Dot Size and Color Relative to Some Values") + 
     theme(legend.position = "none")
print(p)
```

<img src="rgraphics_files/figure-html/qplot_scatter_plot_colored_groups-1.png" width="672" />

#### Regression line


```r
set.seed(1410)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
p <- qplot(carat, price, data = dsmall) +
           geom_smooth(method="lm")
print(p)
```

<img src="rgraphics_files/figure-html/qplot_scatter_regression_line-1.png" width="672" />

#### Local regression curve (loess)


```r
p <- qplot(carat, price, data=dsmall, geom=c("point", "smooth")) 
print(p) # Setting se=FALSE removes error shade
```

<img src="rgraphics_files/figure-html/qplot_scatter_regression_loess-1.png" width="672" />

### `ggplot` Function

- More important than `qplot` to access full functionality of `ggplot2`
- Main arguments
    - data set, usually a `data.frame` or `tibble`
	- aesthetic mappings provided by `aes` function
- General `ggplot` syntax
    - `ggplot(data, aes(...)) + geom() + ... + stat() + ...`
- Layer specifications
    - `geom(mapping, data, ..., geom, position)`
    - `stat(mapping, data, ..., stat, position)`
- Additional components
    - `scales`
	- `coordinates`
	- `facet`
- `aes()` mappings can be passed on to all components (`ggplot, geom`, etc.). Effects are global when passed on to `ggplot()` and local for other components.
    - `x, y`
	- `color`: grouping vector (factor) 
	- `group`: grouping vector (factor)

#### Changing Plotting Themes in `ggplot`

- Theme settings can be accessed with `theme_get()`
- Their settings can be changed with `theme()`

Example how to change background color to white
		

```r
... + theme(panel.background=element_rect(fill = "white", colour = "black")) 
```

#### Storing `ggplot` Specifications

Plots and layers can be stored in variables


```r
p <- ggplot(dsmall, aes(carat, price)) + geom_point() 
p # or print(p)
```

Returns information about data and aesthetic mappings followed by each layer


```r
summary(p) 
```

Print dots with different sizes and colors


```r
bestfit <- geom_smooth(method = "lm", se = F, color = alpha("steelblue", 0.5), size = 2)
p + bestfit # Plot with custom regression line
```

Syntax to pass on other data sets


```r
p %+% diamonds[sample(nrow(diamonds), 100),] 
```

Saves plot stored in variable `p` to file


```r
ggsave(p, file="myplot.pdf") 
```

Standard R export functons for graphics work as well (see [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rgraphics/rgraphics/#saving-graphics-to-files)). 

### `ggplot`: scatter plots

#### Basic example


```r
set.seed(1410)  
dsmall <- as.data.frame(diamonds[sample(nrow(diamonds), 1000), ])
p <- ggplot(dsmall, aes(carat, price, color=color)) + 
            geom_point(size=4)
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_scatter_plot1-1.png" width="672" />

Interactive version of above plot can be generated with the `ggplotly` function from 
the `plotly` package.


```r
library(plotly)
ggplotly(p)
```

```{=html}
<div id="htmlwidget-7d670bfea5cc1c5d088c" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-7d670bfea5cc1c5d088c">{"x":{"data":[{"x":[0.32,0.34,1.54,1.52,1.51,0.52,0.32,0.52,0.66,1,1,0.5,0.47,1.09,0.35,1.01,0.53,0.42,0.5,2.03,1,0.61,0.9,0.56,0.9,0.28,0.33,1.34,0.52,0.79,0.5,0.7,0.27,0.9,0.31,1.5,0.57,0.57,0.61,0.37,0.9,0.32,0.77,0.3,0.33,0.65,0.57,0.9,1.52,0.42,1.08,0.23,0.32,1.18,0.31,0.48,0.9,0.53,1.01,0.71,1.02,0.9,0.79,0.85,0.72,0.41,0.5,0.93,0.44,0.31,0.43,0.39,0.33,0.24,0.52,0.32,1,0.43,0.7,0.37,0.26,0.7,0.74,1.2,1.03,1.05,1.01,1.04,1.71,0.31,1.13,0.51,0.44,0.38,0.46,0.41,0.41,1.11,0.36,0.31,0.24,0.25,0.32,1.58,1.5,2.03,0.23,0.54,0.58,0.36,0.35,1,0.36,0.41,0.53,0.5,0.79,0.31,1.03,0.52,0.8,1,0.36,0.38,0.6,0.36,1.03,0.7,1.02,0.31,0.74,0.31,0.92,0.34,0.91,0.7,0.23],"y":[589,805,16889,7702,12738,1605,756,1919,3915,3965,16073,1272,1270,5164,829,5543,1956,1621,958,13858,5286,3397,4435,1243,3205,452,492,6303,1774,3328,1589,2355,622,4068,1075,15240,2015,1262,2266,857,3812,756,3428,670,492,4275,2047,6872,9789,1031,5758,530,758,5168,766,2329,3755,2266,5939,2516,4912,8239,2823,2737,2951,920,1354,4140,1235,942,774,1170,492,367,3284,657,5775,919,2596,1124,554,2218,3138,4203,4845,7965,7303,6115,11737,982,4406,1431,1294,871,1262,1295,1018,6241,1064,732,552,548,972,13329,12060,13229,425,1786,1879,689,1026,6762,851,791,1412,1436,3270,734,5804,2217,3312,4786,713,1257,2010,1042,8645,1770,5918,593,2906,571,3916,803,5535,2804,680],"text":["carat: 0.32<br />price:   589<br />color: D","carat: 0.34<br />price:   805<br />color: D","carat: 1.54<br />price: 16889<br />color: D","carat: 1.52<br />price:  7702<br />color: D","carat: 1.51<br />price: 12738<br />color: D","carat: 0.52<br />price:  1605<br />color: D","carat: 0.32<br />price:   756<br />color: D","carat: 0.52<br />price:  1919<br />color: D","carat: 0.66<br />price:  3915<br />color: D","carat: 1.00<br />price:  3965<br />color: D","carat: 1.00<br />price: 16073<br />color: D","carat: 0.50<br />price:  1272<br />color: D","carat: 0.47<br />price:  1270<br />color: D","carat: 1.09<br />price:  5164<br />color: D","carat: 0.35<br />price:   829<br />color: D","carat: 1.01<br />price:  5543<br />color: D","carat: 0.53<br />price:  1956<br />color: D","carat: 0.42<br />price:  1621<br />color: D","carat: 0.50<br />price:   958<br />color: D","carat: 2.03<br />price: 13858<br />color: D","carat: 1.00<br />price:  5286<br />color: D","carat: 0.61<br />price:  3397<br />color: D","carat: 0.90<br />price:  4435<br />color: D","carat: 0.56<br />price:  1243<br />color: D","carat: 0.90<br />price:  3205<br />color: D","carat: 0.28<br />price:   452<br />color: D","carat: 0.33<br />price:   492<br />color: D","carat: 1.34<br />price:  6303<br />color: D","carat: 0.52<br />price:  1774<br />color: D","carat: 0.79<br />price:  3328<br />color: D","carat: 0.50<br />price:  1589<br />color: D","carat: 0.70<br />price:  2355<br />color: D","carat: 0.27<br />price:   622<br />color: D","carat: 0.90<br />price:  4068<br />color: D","carat: 0.31<br />price:  1075<br />color: D","carat: 1.50<br />price: 15240<br />color: D","carat: 0.57<br />price:  2015<br />color: D","carat: 0.57<br />price:  1262<br />color: D","carat: 0.61<br />price:  2266<br />color: D","carat: 0.37<br />price:   857<br />color: D","carat: 0.90<br />price:  3812<br />color: D","carat: 0.32<br />price:   756<br />color: D","carat: 0.77<br />price:  3428<br />color: D","carat: 0.30<br />price:   670<br />color: D","carat: 0.33<br />price:   492<br />color: D","carat: 0.65<br />price:  4275<br />color: D","carat: 0.57<br />price:  2047<br />color: D","carat: 0.90<br />price:  6872<br />color: D","carat: 1.52<br />price:  9789<br />color: D","carat: 0.42<br />price:  1031<br />color: D","carat: 1.08<br />price:  5758<br />color: D","carat: 0.23<br />price:   530<br />color: D","carat: 0.32<br />price:   758<br />color: D","carat: 1.18<br />price:  5168<br />color: D","carat: 0.31<br />price:   766<br />color: D","carat: 0.48<br />price:  2329<br />color: D","carat: 0.90<br />price:  3755<br />color: D","carat: 0.53<br />price:  2266<br />color: D","carat: 1.01<br />price:  5939<br />color: D","carat: 0.71<br />price:  2516<br />color: D","carat: 1.02<br />price:  4912<br />color: D","carat: 0.90<br />price:  8239<br />color: D","carat: 0.79<br />price:  2823<br />color: D","carat: 0.85<br />price:  2737<br />color: D","carat: 0.72<br />price:  2951<br />color: D","carat: 0.41<br />price:   920<br />color: D","carat: 0.50<br />price:  1354<br />color: D","carat: 0.93<br />price:  4140<br />color: D","carat: 0.44<br />price:  1235<br />color: D","carat: 0.31<br />price:   942<br />color: D","carat: 0.43<br />price:   774<br />color: D","carat: 0.39<br />price:  1170<br />color: D","carat: 0.33<br />price:   492<br />color: D","carat: 0.24<br />price:   367<br />color: D","carat: 0.52<br />price:  3284<br />color: D","carat: 0.32<br />price:   657<br />color: D","carat: 1.00<br />price:  5775<br />color: D","carat: 0.43<br />price:   919<br />color: D","carat: 0.70<br />price:  2596<br />color: D","carat: 0.37<br />price:  1124<br />color: D","carat: 0.26<br />price:   554<br />color: D","carat: 0.70<br />price:  2218<br />color: D","carat: 0.74<br />price:  3138<br />color: D","carat: 1.20<br />price:  4203<br />color: D","carat: 1.03<br />price:  4845<br />color: D","carat: 1.05<br />price:  7965<br />color: D","carat: 1.01<br />price:  7303<br />color: D","carat: 1.04<br />price:  6115<br />color: D","carat: 1.71<br />price: 11737<br />color: D","carat: 0.31<br />price:   982<br />color: D","carat: 1.13<br />price:  4406<br />color: D","carat: 0.51<br />price:  1431<br />color: D","carat: 0.44<br />price:  1294<br />color: D","carat: 0.38<br />price:   871<br />color: D","carat: 0.46<br />price:  1262<br />color: D","carat: 0.41<br />price:  1295<br />color: D","carat: 0.41<br />price:  1018<br />color: D","carat: 1.11<br />price:  6241<br />color: D","carat: 0.36<br />price:  1064<br />color: D","carat: 0.31<br />price:   732<br />color: D","carat: 0.24<br />price:   552<br />color: D","carat: 0.25<br />price:   548<br />color: D","carat: 0.32<br />price:   972<br />color: D","carat: 1.58<br />price: 13329<br />color: D","carat: 1.50<br />price: 12060<br />color: D","carat: 2.03<br />price: 13229<br />color: D","carat: 0.23<br />price:   425<br />color: D","carat: 0.54<br />price:  1786<br />color: D","carat: 0.58<br />price:  1879<br />color: D","carat: 0.36<br />price:   689<br />color: D","carat: 0.35<br />price:  1026<br />color: D","carat: 1.00<br />price:  6762<br />color: D","carat: 0.36<br />price:   851<br />color: D","carat: 0.41<br />price:   791<br />color: D","carat: 0.53<br />price:  1412<br />color: D","carat: 0.50<br />price:  1436<br />color: D","carat: 0.79<br />price:  3270<br />color: D","carat: 0.31<br />price:   734<br />color: D","carat: 1.03<br />price:  5804<br />color: D","carat: 0.52<br />price:  2217<br />color: D","carat: 0.80<br />price:  3312<br />color: D","carat: 1.00<br />price:  4786<br />color: D","carat: 0.36<br />price:   713<br />color: D","carat: 0.38<br />price:  1257<br />color: D","carat: 0.60<br />price:  2010<br />color: D","carat: 0.36<br />price:  1042<br />color: D","carat: 1.03<br />price:  8645<br />color: D","carat: 0.70<br />price:  1770<br />color: D","carat: 1.02<br />price:  5918<br />color: D","carat: 0.31<br />price:   593<br />color: D","carat: 0.74<br />price:  2906<br />color: D","carat: 0.31<br />price:   571<br />color: D","carat: 0.92<br />price:  3916<br />color: D","carat: 0.34<br />price:   803<br />color: D","carat: 0.91<br />price:  5535<br />color: D","carat: 0.70<br />price:  2804<br />color: D","carat: 0.23<br />price:   680<br />color: D"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(68,1,84,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(68,1,84,1)"}},"hoveron":"points","name":"D","legendgroup":"D","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.25,0.9,0.75,1.01,1.01,0.54,0.39,0.7,1.01,0.34,0.51,0.41,1.02,0.38,0.43,0.6,1.71,0.71,0.31,0.33,1.51,0.31,2.01,0.91,1.09,0.39,0.5,1.51,0.46,0.35,0.36,0.56,1.26,1.12,0.56,0.63,0.35,0.33,1.66,0.32,0.41,0.37,1.02,0.26,0.74,0.72,0.36,1.63,0.34,0.41,0.52,0.31,2,0.7,0.43,0.38,1.52,0.95,0.41,0.51,0.55,1.01,0.78,0.31,0.41,0.33,1.51,0.34,0.41,0.63,0.34,1.09,0.53,0.29,0.38,0.33,0.89,0.41,1,0.3,0.75,0.32,1.21,1,1.12,0.24,0.5,0.34,0.73,1.11,0.46,0.58,0.41,0.5,0.5,0.71,0.78,0.32,0.5,0.57,0.31,0.33,0.43,0.9,0.23,0.55,0.32,1,0.51,1,0.28,0.75,0.33,0.5,1.57,0.71,2,2.01,1.05,0.26,0.3,0.9,1,0.74,0.76,0.55,0.3,0.52,0.71,0.77,0.32,0.4,1.7,0.31,0.3,1.42,0.37,1.02,0.72,1.01,0.28,0.59,1.01,0.72,1,1,1.01,1.36,1,0.28,0.48,1.6,0.52,0.31,0.7,0.73,0.51,0.23,0.52,0.3,1.22,0.31,0.33,1.6,0.79,0.35,1.03,0.7,0.26,1.01,0.54,1.02,0.7,0.71,0.31,0.32,1.11,0.52,0.34,0.7,1.01,1.01],"y":[783,3239,3231,4912,8976,1250,601,3504,5279,556,1656,1249,7277,1267,1433,2774,14319,3360,816,804,12530,879,16733,3743,10111,985,1217,9302,1661,767,933,1915,6085,5546,1915,2765,921,957,10553,768,1079,1041,4465,556,2999,2498,738,17825,765,1153,1601,872,15840,3419,716,866,12418,4738,638,1340,1551,6128,2379,872,981,693,9343,612,1431,2565,1033,5127,1438,619,1053,928,3872,1187,4435,895,2940,1095,6708,6315,7589,362,1240,1014,2409,4532,1070,1641,818,1351,1410,2190,3444,1095,1746,1613,872,1002,669,3662,530,1448,1080,5330,1831,3763,642,2674,928,1259,10362,2891,16064,15618,5494,599,526,4092,9139,2464,2580,1675,564,1229,2329,3291,954,1125,9596,936,895,8019,999,11128,3235,4608,470,2177,5751,3550,9224,4916,6268,7677,4805,688,1323,10798,1646,958,2952,2886,1443,530,1793,578,7395,734,743,12467,3539,767,5010,2757,554,6968,1163,4540,2415,2921,1133,528,6800,1689,881,2952,5366,4751],"text":["carat: 0.25<br />price:   783<br />color: E","carat: 0.90<br />price:  3239<br />color: E","carat: 0.75<br />price:  3231<br />color: E","carat: 1.01<br />price:  4912<br />color: E","carat: 1.01<br />price:  8976<br />color: E","carat: 0.54<br />price:  1250<br />color: E","carat: 0.39<br />price:   601<br />color: E","carat: 0.70<br />price:  3504<br />color: E","carat: 1.01<br />price:  5279<br />color: E","carat: 0.34<br />price:   556<br />color: E","carat: 0.51<br />price:  1656<br />color: E","carat: 0.41<br />price:  1249<br />color: E","carat: 1.02<br />price:  7277<br />color: E","carat: 0.38<br />price:  1267<br />color: E","carat: 0.43<br />price:  1433<br />color: E","carat: 0.60<br />price:  2774<br />color: E","carat: 1.71<br />price: 14319<br />color: E","carat: 0.71<br />price:  3360<br />color: E","carat: 0.31<br />price:   816<br />color: E","carat: 0.33<br />price:   804<br />color: E","carat: 1.51<br />price: 12530<br />color: E","carat: 0.31<br />price:   879<br />color: E","carat: 2.01<br />price: 16733<br />color: E","carat: 0.91<br />price:  3743<br />color: E","carat: 1.09<br />price: 10111<br />color: E","carat: 0.39<br />price:   985<br />color: E","carat: 0.50<br />price:  1217<br />color: E","carat: 1.51<br />price:  9302<br />color: E","carat: 0.46<br />price:  1661<br />color: E","carat: 0.35<br />price:   767<br />color: E","carat: 0.36<br />price:   933<br />color: E","carat: 0.56<br />price:  1915<br />color: E","carat: 1.26<br />price:  6085<br />color: E","carat: 1.12<br />price:  5546<br />color: E","carat: 0.56<br />price:  1915<br />color: E","carat: 0.63<br />price:  2765<br />color: E","carat: 0.35<br />price:   921<br />color: E","carat: 0.33<br />price:   957<br />color: E","carat: 1.66<br />price: 10553<br />color: E","carat: 0.32<br />price:   768<br />color: E","carat: 0.41<br />price:  1079<br />color: E","carat: 0.37<br />price:  1041<br />color: E","carat: 1.02<br />price:  4465<br />color: E","carat: 0.26<br />price:   556<br />color: E","carat: 0.74<br />price:  2999<br />color: E","carat: 0.72<br />price:  2498<br />color: E","carat: 0.36<br />price:   738<br />color: E","carat: 1.63<br />price: 17825<br />color: E","carat: 0.34<br />price:   765<br />color: E","carat: 0.41<br />price:  1153<br />color: E","carat: 0.52<br />price:  1601<br />color: E","carat: 0.31<br />price:   872<br />color: E","carat: 2.00<br />price: 15840<br />color: E","carat: 0.70<br />price:  3419<br />color: E","carat: 0.43<br />price:   716<br />color: E","carat: 0.38<br />price:   866<br />color: E","carat: 1.52<br />price: 12418<br />color: E","carat: 0.95<br />price:  4738<br />color: E","carat: 0.41<br />price:   638<br />color: E","carat: 0.51<br />price:  1340<br />color: E","carat: 0.55<br />price:  1551<br />color: E","carat: 1.01<br />price:  6128<br />color: E","carat: 0.78<br />price:  2379<br />color: E","carat: 0.31<br />price:   872<br />color: E","carat: 0.41<br />price:   981<br />color: E","carat: 0.33<br />price:   693<br />color: E","carat: 1.51<br />price:  9343<br />color: E","carat: 0.34<br />price:   612<br />color: E","carat: 0.41<br />price:  1431<br />color: E","carat: 0.63<br />price:  2565<br />color: E","carat: 0.34<br />price:  1033<br />color: E","carat: 1.09<br />price:  5127<br />color: E","carat: 0.53<br />price:  1438<br />color: E","carat: 0.29<br />price:   619<br />color: E","carat: 0.38<br />price:  1053<br />color: E","carat: 0.33<br />price:   928<br />color: E","carat: 0.89<br />price:  3872<br />color: E","carat: 0.41<br />price:  1187<br />color: E","carat: 1.00<br />price:  4435<br />color: E","carat: 0.30<br />price:   895<br />color: E","carat: 0.75<br />price:  2940<br />color: E","carat: 0.32<br />price:  1095<br />color: E","carat: 1.21<br />price:  6708<br />color: E","carat: 1.00<br />price:  6315<br />color: E","carat: 1.12<br />price:  7589<br />color: E","carat: 0.24<br />price:   362<br />color: E","carat: 0.50<br />price:  1240<br />color: E","carat: 0.34<br />price:  1014<br />color: E","carat: 0.73<br />price:  2409<br />color: E","carat: 1.11<br />price:  4532<br />color: E","carat: 0.46<br />price:  1070<br />color: E","carat: 0.58<br />price:  1641<br />color: E","carat: 0.41<br />price:   818<br />color: E","carat: 0.50<br />price:  1351<br />color: E","carat: 0.50<br />price:  1410<br />color: E","carat: 0.71<br />price:  2190<br />color: E","carat: 0.78<br />price:  3444<br />color: E","carat: 0.32<br />price:  1095<br />color: E","carat: 0.50<br />price:  1746<br />color: E","carat: 0.57<br />price:  1613<br />color: E","carat: 0.31<br />price:   872<br />color: E","carat: 0.33<br />price:  1002<br />color: E","carat: 0.43<br />price:   669<br />color: E","carat: 0.90<br />price:  3662<br />color: E","carat: 0.23<br />price:   530<br />color: E","carat: 0.55<br />price:  1448<br />color: E","carat: 0.32<br />price:  1080<br />color: E","carat: 1.00<br />price:  5330<br />color: E","carat: 0.51<br />price:  1831<br />color: E","carat: 1.00<br />price:  3763<br />color: E","carat: 0.28<br />price:   642<br />color: E","carat: 0.75<br />price:  2674<br />color: E","carat: 0.33<br />price:   928<br />color: E","carat: 0.50<br />price:  1259<br />color: E","carat: 1.57<br />price: 10362<br />color: E","carat: 0.71<br />price:  2891<br />color: E","carat: 2.00<br />price: 16064<br />color: E","carat: 2.01<br />price: 15618<br />color: E","carat: 1.05<br />price:  5494<br />color: E","carat: 0.26<br />price:   599<br />color: E","carat: 0.30<br />price:   526<br />color: E","carat: 0.90<br />price:  4092<br />color: E","carat: 1.00<br />price:  9139<br />color: E","carat: 0.74<br />price:  2464<br />color: E","carat: 0.76<br />price:  2580<br />color: E","carat: 0.55<br />price:  1675<br />color: E","carat: 0.30<br />price:   564<br />color: E","carat: 0.52<br />price:  1229<br />color: E","carat: 0.71<br />price:  2329<br />color: E","carat: 0.77<br />price:  3291<br />color: E","carat: 0.32<br />price:   954<br />color: E","carat: 0.40<br />price:  1125<br />color: E","carat: 1.70<br />price:  9596<br />color: E","carat: 0.31<br />price:   936<br />color: E","carat: 0.30<br />price:   895<br />color: E","carat: 1.42<br />price:  8019<br />color: E","carat: 0.37<br />price:   999<br />color: E","carat: 1.02<br />price: 11128<br />color: E","carat: 0.72<br />price:  3235<br />color: E","carat: 1.01<br />price:  4608<br />color: E","carat: 0.28<br />price:   470<br />color: E","carat: 0.59<br />price:  2177<br />color: E","carat: 1.01<br />price:  5751<br />color: E","carat: 0.72<br />price:  3550<br />color: E","carat: 1.00<br />price:  9224<br />color: E","carat: 1.00<br />price:  4916<br />color: E","carat: 1.01<br />price:  6268<br />color: E","carat: 1.36<br />price:  7677<br />color: E","carat: 1.00<br />price:  4805<br />color: E","carat: 0.28<br />price:   688<br />color: E","carat: 0.48<br />price:  1323<br />color: E","carat: 1.60<br />price: 10798<br />color: E","carat: 0.52<br />price:  1646<br />color: E","carat: 0.31<br />price:   958<br />color: E","carat: 0.70<br />price:  2952<br />color: E","carat: 0.73<br />price:  2886<br />color: E","carat: 0.51<br />price:  1443<br />color: E","carat: 0.23<br />price:   530<br />color: E","carat: 0.52<br />price:  1793<br />color: E","carat: 0.30<br />price:   578<br />color: E","carat: 1.22<br />price:  7395<br />color: E","carat: 0.31<br />price:   734<br />color: E","carat: 0.33<br />price:   743<br />color: E","carat: 1.60<br />price: 12467<br />color: E","carat: 0.79<br />price:  3539<br />color: E","carat: 0.35<br />price:   767<br />color: E","carat: 1.03<br />price:  5010<br />color: E","carat: 0.70<br />price:  2757<br />color: E","carat: 0.26<br />price:   554<br />color: E","carat: 1.01<br />price:  6968<br />color: E","carat: 0.54<br />price:  1163<br />color: E","carat: 1.02<br />price:  4540<br />color: E","carat: 0.70<br />price:  2415<br />color: E","carat: 0.71<br />price:  2921<br />color: E","carat: 0.31<br />price:  1133<br />color: E","carat: 0.32<br />price:   528<br />color: E","carat: 1.11<br />price:  6800<br />color: E","carat: 0.52<br />price:  1689<br />color: E","carat: 0.34<br />price:   881<br />color: E","carat: 0.70<br />price:  2952<br />color: E","carat: 1.01<br />price:  5366<br />color: E","carat: 1.01<br />price:  4751<br />color: E"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(68,58,131,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(68,58,131,1)"}},"hoveron":"points","name":"E","legendgroup":"E","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1.25,1.04,0.9,0.69,0.58,1.11,0.73,0.8,1.16,0.71,0.9,1.02,0.91,0.3,0.57,0.91,1.59,2.07,0.32,0.79,0.39,1.04,1.29,0.83,0.51,0.32,0.3,1.02,1.05,0.72,1.52,0.43,1.5,1.21,0.42,0.9,0.59,0.3,0.9,0.35,1.2,0.4,0.4,0.4,0.52,0.51,0.41,1.07,0.7,0.54,0.56,0.7,0.43,0.31,1,0.74,1.23,0.7,0.3,0.74,1.62,0.32,1.01,0.42,0.3,0.41,1.23,0.81,0.5,0.3,0.64,0.51,1.24,0.34,0.52,1.11,0.3,0.72,0.41,0.25,0.41,0.37,1.27,0.7,1.27,1.5,0.31,1.01,1.01,1.33,0.73,0.33,0.57,0.31,1.29,1.14,0.75,0.71,0.3,0.5,0.3,1,0.54,0.54,0.41,0.7,1.01,0.35,0.9,0.29,2.02,1.51,2.05,1.01,0.27,0.28,1.74,2,1.22,1.01,0.95,0.3,1.05,0.3,0.32,0.24,1.22,0.7,0.56,0.38,0.32,0.71,0.26,1.09,0.36,0.42,1.18,0.55,0.71,0.71,0.53,0.31,0.48,0.9,0.27,1.32,0.52,1,2.02,1.21,0.55,0.33,0.31,0.33,0.4,0.71,0.6,1.03,0.24,1.01,0.7,0.42,1.47,0.71,0.31,0.24,0.3,0.36,0.51,1,0.73,0.61,0.73,0.41,0.37,1.01,1,2.02,0.55,0.31,0.55,0.31,0.61,0.33,0.23,0.57],"y":[8764,5650,3470,1951,1401,10141,2355,3293,10809,2999,4441,4162,4423,945,1760,3489,9972,18611,828,3806,1183,4321,7126,1968,1574,828,776,4043,6933,2607,11405,1169,12389,8879,910,3619,1195,658,2964,872,9586,982,1404,945,1227,2041,1356,6828,2662,1385,2554,3154,1422,625,5221,2319,11141,3300,886,3454,11733,828,6702,926,886,755,9660,2761,1849,730,2074,1334,10033,784,1250,7823,795,2829,1007,633,870,957,5308,2373,13823,10962,917,4791,6083,8163,2527,1040,1397,625,4977,7293,2648,2809,873,1333,694,6927,1937,1836,1323,1134,5662,906,4232,555,13453,11302,15017,4416,820,825,6821,16580,6518,5345,2867,776,7276,886,977,540,10100,1806,1950,870,645,2409,597,8422,781,773,12175,1692,2416,2590,1610,921,1279,3898,560,12736,1600,4037,18274,8691,1698,957,394,841,1323,2130,1367,5356,370,8552,3270,1031,8115,2103,802,678,601,803,1100,4155,2951,2634,2801,1107,957,7122,6424,15261,1580,901,1692,802,1547,801,402,1746],"text":["carat: 1.25<br />price:  8764<br />color: F","carat: 1.04<br />price:  5650<br />color: F","carat: 0.90<br />price:  3470<br />color: F","carat: 0.69<br />price:  1951<br />color: F","carat: 0.58<br />price:  1401<br />color: F","carat: 1.11<br />price: 10141<br />color: F","carat: 0.73<br />price:  2355<br />color: F","carat: 0.80<br />price:  3293<br />color: F","carat: 1.16<br />price: 10809<br />color: F","carat: 0.71<br />price:  2999<br />color: F","carat: 0.90<br />price:  4441<br />color: F","carat: 1.02<br />price:  4162<br />color: F","carat: 0.91<br />price:  4423<br />color: F","carat: 0.30<br />price:   945<br />color: F","carat: 0.57<br />price:  1760<br />color: F","carat: 0.91<br />price:  3489<br />color: F","carat: 1.59<br />price:  9972<br />color: F","carat: 2.07<br />price: 18611<br />color: F","carat: 0.32<br />price:   828<br />color: F","carat: 0.79<br />price:  3806<br />color: F","carat: 0.39<br />price:  1183<br />color: F","carat: 1.04<br />price:  4321<br />color: F","carat: 1.29<br />price:  7126<br />color: F","carat: 0.83<br />price:  1968<br />color: F","carat: 0.51<br />price:  1574<br />color: F","carat: 0.32<br />price:   828<br />color: F","carat: 0.30<br />price:   776<br />color: F","carat: 1.02<br />price:  4043<br />color: F","carat: 1.05<br />price:  6933<br />color: F","carat: 0.72<br />price:  2607<br />color: F","carat: 1.52<br />price: 11405<br />color: F","carat: 0.43<br />price:  1169<br />color: F","carat: 1.50<br />price: 12389<br />color: F","carat: 1.21<br />price:  8879<br />color: F","carat: 0.42<br />price:   910<br />color: F","carat: 0.90<br />price:  3619<br />color: F","carat: 0.59<br />price:  1195<br />color: F","carat: 0.30<br />price:   658<br />color: F","carat: 0.90<br />price:  2964<br />color: F","carat: 0.35<br />price:   872<br />color: F","carat: 1.20<br />price:  9586<br />color: F","carat: 0.40<br />price:   982<br />color: F","carat: 0.40<br />price:  1404<br />color: F","carat: 0.40<br />price:   945<br />color: F","carat: 0.52<br />price:  1227<br />color: F","carat: 0.51<br />price:  2041<br />color: F","carat: 0.41<br />price:  1356<br />color: F","carat: 1.07<br />price:  6828<br />color: F","carat: 0.70<br />price:  2662<br />color: F","carat: 0.54<br />price:  1385<br />color: F","carat: 0.56<br />price:  2554<br />color: F","carat: 0.70<br />price:  3154<br />color: F","carat: 0.43<br />price:  1422<br />color: F","carat: 0.31<br />price:   625<br />color: F","carat: 1.00<br />price:  5221<br />color: F","carat: 0.74<br />price:  2319<br />color: F","carat: 1.23<br />price: 11141<br />color: F","carat: 0.70<br />price:  3300<br />color: F","carat: 0.30<br />price:   886<br />color: F","carat: 0.74<br />price:  3454<br />color: F","carat: 1.62<br />price: 11733<br />color: F","carat: 0.32<br />price:   828<br />color: F","carat: 1.01<br />price:  6702<br />color: F","carat: 0.42<br />price:   926<br />color: F","carat: 0.30<br />price:   886<br />color: F","carat: 0.41<br />price:   755<br />color: F","carat: 1.23<br />price:  9660<br />color: F","carat: 0.81<br />price:  2761<br />color: F","carat: 0.50<br />price:  1849<br />color: F","carat: 0.30<br />price:   730<br />color: F","carat: 0.64<br />price:  2074<br />color: F","carat: 0.51<br />price:  1334<br />color: F","carat: 1.24<br />price: 10033<br />color: F","carat: 0.34<br />price:   784<br />color: F","carat: 0.52<br />price:  1250<br />color: F","carat: 1.11<br />price:  7823<br />color: F","carat: 0.30<br />price:   795<br />color: F","carat: 0.72<br />price:  2829<br />color: F","carat: 0.41<br />price:  1007<br />color: F","carat: 0.25<br />price:   633<br />color: F","carat: 0.41<br />price:   870<br />color: F","carat: 0.37<br />price:   957<br />color: F","carat: 1.27<br />price:  5308<br />color: F","carat: 0.70<br />price:  2373<br />color: F","carat: 1.27<br />price: 13823<br />color: F","carat: 1.50<br />price: 10962<br />color: F","carat: 0.31<br />price:   917<br />color: F","carat: 1.01<br />price:  4791<br />color: F","carat: 1.01<br />price:  6083<br />color: F","carat: 1.33<br />price:  8163<br />color: F","carat: 0.73<br />price:  2527<br />color: F","carat: 0.33<br />price:  1040<br />color: F","carat: 0.57<br />price:  1397<br />color: F","carat: 0.31<br />price:   625<br />color: F","carat: 1.29<br />price:  4977<br />color: F","carat: 1.14<br />price:  7293<br />color: F","carat: 0.75<br />price:  2648<br />color: F","carat: 0.71<br />price:  2809<br />color: F","carat: 0.30<br />price:   873<br />color: F","carat: 0.50<br />price:  1333<br />color: F","carat: 0.30<br />price:   694<br />color: F","carat: 1.00<br />price:  6927<br />color: F","carat: 0.54<br />price:  1937<br />color: F","carat: 0.54<br />price:  1836<br />color: F","carat: 0.41<br />price:  1323<br />color: F","carat: 0.70<br />price:  1134<br />color: F","carat: 1.01<br />price:  5662<br />color: F","carat: 0.35<br />price:   906<br />color: F","carat: 0.90<br />price:  4232<br />color: F","carat: 0.29<br />price:   555<br />color: F","carat: 2.02<br />price: 13453<br />color: F","carat: 1.51<br />price: 11302<br />color: F","carat: 2.05<br />price: 15017<br />color: F","carat: 1.01<br />price:  4416<br />color: F","carat: 0.27<br />price:   820<br />color: F","carat: 0.28<br />price:   825<br />color: F","carat: 1.74<br />price:  6821<br />color: F","carat: 2.00<br />price: 16580<br />color: F","carat: 1.22<br />price:  6518<br />color: F","carat: 1.01<br />price:  5345<br />color: F","carat: 0.95<br />price:  2867<br />color: F","carat: 0.30<br />price:   776<br />color: F","carat: 1.05<br />price:  7276<br />color: F","carat: 0.30<br />price:   886<br />color: F","carat: 0.32<br />price:   977<br />color: F","carat: 0.24<br />price:   540<br />color: F","carat: 1.22<br />price: 10100<br />color: F","carat: 0.70<br />price:  1806<br />color: F","carat: 0.56<br />price:  1950<br />color: F","carat: 0.38<br />price:   870<br />color: F","carat: 0.32<br />price:   645<br />color: F","carat: 0.71<br />price:  2409<br />color: F","carat: 0.26<br />price:   597<br />color: F","carat: 1.09<br />price:  8422<br />color: F","carat: 0.36<br />price:   781<br />color: F","carat: 0.42<br />price:   773<br />color: F","carat: 1.18<br />price: 12175<br />color: F","carat: 0.55<br />price:  1692<br />color: F","carat: 0.71<br />price:  2416<br />color: F","carat: 0.71<br />price:  2590<br />color: F","carat: 0.53<br />price:  1610<br />color: F","carat: 0.31<br />price:   921<br />color: F","carat: 0.48<br />price:  1279<br />color: F","carat: 0.90<br />price:  3898<br />color: F","carat: 0.27<br />price:   560<br />color: F","carat: 1.32<br />price: 12736<br />color: F","carat: 0.52<br />price:  1600<br />color: F","carat: 1.00<br />price:  4037<br />color: F","carat: 2.02<br />price: 18274<br />color: F","carat: 1.21<br />price:  8691<br />color: F","carat: 0.55<br />price:  1698<br />color: F","carat: 0.33<br />price:   957<br />color: F","carat: 0.31<br />price:   394<br />color: F","carat: 0.33<br />price:   841<br />color: F","carat: 0.40<br />price:  1323<br />color: F","carat: 0.71<br />price:  2130<br />color: F","carat: 0.60<br />price:  1367<br />color: F","carat: 1.03<br />price:  5356<br />color: F","carat: 0.24<br />price:   370<br />color: F","carat: 1.01<br />price:  8552<br />color: F","carat: 0.70<br />price:  3270<br />color: F","carat: 0.42<br />price:  1031<br />color: F","carat: 1.47<br />price:  8115<br />color: F","carat: 0.71<br />price:  2103<br />color: F","carat: 0.31<br />price:   802<br />color: F","carat: 0.24<br />price:   678<br />color: F","carat: 0.30<br />price:   601<br />color: F","carat: 0.36<br />price:   803<br />color: F","carat: 0.51<br />price:  1100<br />color: F","carat: 1.00<br />price:  4155<br />color: F","carat: 0.73<br />price:  2951<br />color: F","carat: 0.61<br />price:  2634<br />color: F","carat: 0.73<br />price:  2801<br />color: F","carat: 0.41<br />price:  1107<br />color: F","carat: 0.37<br />price:   957<br />color: F","carat: 1.01<br />price:  7122<br />color: F","carat: 1.00<br />price:  6424<br />color: F","carat: 2.02<br />price: 15261<br />color: F","carat: 0.55<br />price:  1580<br />color: F","carat: 0.31<br />price:   901<br />color: F","carat: 0.55<br />price:  1692<br />color: F","carat: 0.31<br />price:   802<br />color: F","carat: 0.61<br />price:  1547<br />color: F","carat: 0.33<br />price:   801<br />color: F","carat: 0.23<br />price:   402<br />color: F","carat: 0.57<br />price:  1746<br />color: F"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(49,104,142,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(49,104,142,1)"}},"hoveron":"points","name":"F","legendgroup":"F","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.3,0.9,0.57,0.3,1.5,0.4,0.3,0.41,0.32,1.16,1.62,0.7,1.01,0.71,2.02,1.24,0.7,0.31,1.12,1,1.06,0.3,0.7,1.5,1.2,0.9,1.05,0.32,1.01,0.31,0.53,0.41,1.53,1,1.07,0.5,1.74,0.42,1.05,0.35,0.31,0.3,0.7,1.01,0.26,1.02,1.55,0.3,1.28,1.2,0.31,2,1,1.03,1.25,0.3,2.03,0.9,0.5,1.03,1.54,0.3,1.15,2.18,1.18,0.59,0.73,0.59,0.31,0.33,1.01,0.54,0.38,0.33,0.7,0.33,0.55,1.25,0.4,1.02,0.45,1.51,0.5,1.01,2.14,0.92,0.96,1.04,0.72,0.5,0.73,0.71,0.49,1.51,1.09,1.26,0.51,0.42,0.53,1.03,1.32,1,0.29,1.21,1.3,0.56,1.04,0.3,0.9,0.62,0.32,0.41,0.71,0.54,2.04,0.73,0.5,0.31,1.05,0.42,1,0.43,0.71,0.42,1.06,0.9,0.7,0.32,0.37,0.46,0.5,0.32,1.06,0.39,0.59,0.9,0.41,0.71,0.43,0.82,0.31,1.07,1.01,1.2,1.52,0.24,0.35,0.39,0.54,1.01,0.4,0.43,0.32,0.52,0.25,0.55,1.71,0.33,0.35,0.53,0.41,1.01,0.31,0.31,1.13,0.5,0.33,0.5,0.32,0.41,0.71,0.51,1,0.31,1.53,0.5,1.04,0.74,0.51,0.31,0.31,0.25,1,0.71,0.9,1.7,0.76,0.4,0.32,0.52,2.32,0.44,0.53,1.17,0.31,1.01,0.55,0.32,0.71,1.07,0.72],"y":[878,4314,1617,684,14067,822,673,801,936,9901,15318,2010,4588,2386,17099,6895,2282,544,8072,4816,6325,863,1420,13998,7380,3976,9013,607,6734,605,1273,899,9572,5500,6779,1935,16316,882,7028,659,891,655,2478,8029,394,8489,14106,956,9505,6494,391,12654,6699,5518,11089,650,18680,4133,1787,5814,14438,638,7175,18003,9366,1782,3183,1355,500,1052,6295,2274,983,594,2056,702,2028,9173,1035,5398,1297,14674,1449,5049,12792,3377,3858,6831,2720,1659,2825,3299,1865,11153,3601,8113,1825,1048,1598,6479,9715,3965,712,8774,9540,1819,4200,776,3820,933,842,671,2334,1090,14900,2919,1449,629,3599,1087,7623,929,2443,737,5435,5978,2562,720,746,1560,1020,828,4238,855,1437,5855,839,2492,729,3652,542,7858,4672,5996,12161,559,756,1095,1985,6239,765,1264,698,2012,460,1667,11227,1114,669,1655,1068,4004,583,802,4707,1153,868,1354,645,986,2649,1438,5483,698,13994,1966,9465,3402,1438,816,816,454,4095,2515,3303,18718,2972,870,720,1244,18532,990,1692,9289,789,5741,2016,918,2613,4419,3179],"text":["carat: 0.30<br />price:   878<br />color: G","carat: 0.90<br />price:  4314<br />color: G","carat: 0.57<br />price:  1617<br />color: G","carat: 0.30<br />price:   684<br />color: G","carat: 1.50<br />price: 14067<br />color: G","carat: 0.40<br />price:   822<br />color: G","carat: 0.30<br />price:   673<br />color: G","carat: 0.41<br />price:   801<br />color: G","carat: 0.32<br />price:   936<br />color: G","carat: 1.16<br />price:  9901<br />color: G","carat: 1.62<br />price: 15318<br />color: G","carat: 0.70<br />price:  2010<br />color: G","carat: 1.01<br />price:  4588<br />color: G","carat: 0.71<br />price:  2386<br />color: G","carat: 2.02<br />price: 17099<br />color: G","carat: 1.24<br />price:  6895<br />color: G","carat: 0.70<br />price:  2282<br />color: G","carat: 0.31<br />price:   544<br />color: G","carat: 1.12<br />price:  8072<br />color: G","carat: 1.00<br />price:  4816<br />color: G","carat: 1.06<br />price:  6325<br />color: G","carat: 0.30<br />price:   863<br />color: G","carat: 0.70<br />price:  1420<br />color: G","carat: 1.50<br />price: 13998<br />color: G","carat: 1.20<br />price:  7380<br />color: G","carat: 0.90<br />price:  3976<br />color: G","carat: 1.05<br />price:  9013<br />color: G","carat: 0.32<br />price:   607<br />color: G","carat: 1.01<br />price:  6734<br />color: G","carat: 0.31<br />price:   605<br />color: G","carat: 0.53<br />price:  1273<br />color: G","carat: 0.41<br />price:   899<br />color: G","carat: 1.53<br />price:  9572<br />color: G","carat: 1.00<br />price:  5500<br />color: G","carat: 1.07<br />price:  6779<br />color: G","carat: 0.50<br />price:  1935<br />color: G","carat: 1.74<br />price: 16316<br />color: G","carat: 0.42<br />price:   882<br />color: G","carat: 1.05<br />price:  7028<br />color: G","carat: 0.35<br />price:   659<br />color: G","carat: 0.31<br />price:   891<br />color: G","carat: 0.30<br />price:   655<br />color: G","carat: 0.70<br />price:  2478<br />color: G","carat: 1.01<br />price:  8029<br />color: G","carat: 0.26<br />price:   394<br />color: G","carat: 1.02<br />price:  8489<br />color: G","carat: 1.55<br />price: 14106<br />color: G","carat: 0.30<br />price:   956<br />color: G","carat: 1.28<br />price:  9505<br />color: G","carat: 1.20<br />price:  6494<br />color: G","carat: 0.31<br />price:   391<br />color: G","carat: 2.00<br />price: 12654<br />color: G","carat: 1.00<br />price:  6699<br />color: G","carat: 1.03<br />price:  5518<br />color: G","carat: 1.25<br />price: 11089<br />color: G","carat: 0.30<br />price:   650<br />color: G","carat: 2.03<br />price: 18680<br />color: G","carat: 0.90<br />price:  4133<br />color: G","carat: 0.50<br />price:  1787<br />color: G","carat: 1.03<br />price:  5814<br />color: G","carat: 1.54<br />price: 14438<br />color: G","carat: 0.30<br />price:   638<br />color: G","carat: 1.15<br />price:  7175<br />color: G","carat: 2.18<br />price: 18003<br />color: G","carat: 1.18<br />price:  9366<br />color: G","carat: 0.59<br />price:  1782<br />color: G","carat: 0.73<br />price:  3183<br />color: G","carat: 0.59<br />price:  1355<br />color: G","carat: 0.31<br />price:   500<br />color: G","carat: 0.33<br />price:  1052<br />color: G","carat: 1.01<br />price:  6295<br />color: G","carat: 0.54<br />price:  2274<br />color: G","carat: 0.38<br />price:   983<br />color: G","carat: 0.33<br />price:   594<br />color: G","carat: 0.70<br />price:  2056<br />color: G","carat: 0.33<br />price:   702<br />color: G","carat: 0.55<br />price:  2028<br />color: G","carat: 1.25<br />price:  9173<br />color: G","carat: 0.40<br />price:  1035<br />color: G","carat: 1.02<br />price:  5398<br />color: G","carat: 0.45<br />price:  1297<br />color: G","carat: 1.51<br />price: 14674<br />color: G","carat: 0.50<br />price:  1449<br />color: G","carat: 1.01<br />price:  5049<br />color: G","carat: 2.14<br />price: 12792<br />color: G","carat: 0.92<br />price:  3377<br />color: G","carat: 0.96<br />price:  3858<br />color: G","carat: 1.04<br />price:  6831<br />color: G","carat: 0.72<br />price:  2720<br />color: G","carat: 0.50<br />price:  1659<br />color: G","carat: 0.73<br />price:  2825<br />color: G","carat: 0.71<br />price:  3299<br />color: G","carat: 0.49<br />price:  1865<br />color: G","carat: 1.51<br />price: 11153<br />color: G","carat: 1.09<br />price:  3601<br />color: G","carat: 1.26<br />price:  8113<br />color: G","carat: 0.51<br />price:  1825<br />color: G","carat: 0.42<br />price:  1048<br />color: G","carat: 0.53<br />price:  1598<br />color: G","carat: 1.03<br />price:  6479<br />color: G","carat: 1.32<br />price:  9715<br />color: G","carat: 1.00<br />price:  3965<br />color: G","carat: 0.29<br />price:   712<br />color: G","carat: 1.21<br />price:  8774<br />color: G","carat: 1.30<br />price:  9540<br />color: G","carat: 0.56<br />price:  1819<br />color: G","carat: 1.04<br />price:  4200<br />color: G","carat: 0.30<br />price:   776<br />color: G","carat: 0.90<br />price:  3820<br />color: G","carat: 0.62<br />price:   933<br />color: G","carat: 0.32<br />price:   842<br />color: G","carat: 0.41<br />price:   671<br />color: G","carat: 0.71<br />price:  2334<br />color: G","carat: 0.54<br />price:  1090<br />color: G","carat: 2.04<br />price: 14900<br />color: G","carat: 0.73<br />price:  2919<br />color: G","carat: 0.50<br />price:  1449<br />color: G","carat: 0.31<br />price:   629<br />color: G","carat: 1.05<br />price:  3599<br />color: G","carat: 0.42<br />price:  1087<br />color: G","carat: 1.00<br />price:  7623<br />color: G","carat: 0.43<br />price:   929<br />color: G","carat: 0.71<br />price:  2443<br />color: G","carat: 0.42<br />price:   737<br />color: G","carat: 1.06<br />price:  5435<br />color: G","carat: 0.90<br />price:  5978<br />color: G","carat: 0.70<br />price:  2562<br />color: G","carat: 0.32<br />price:   720<br />color: G","carat: 0.37<br />price:   746<br />color: G","carat: 0.46<br />price:  1560<br />color: G","carat: 0.50<br />price:  1020<br />color: G","carat: 0.32<br />price:   828<br />color: G","carat: 1.06<br />price:  4238<br />color: G","carat: 0.39<br />price:   855<br />color: G","carat: 0.59<br />price:  1437<br />color: G","carat: 0.90<br />price:  5855<br />color: G","carat: 0.41<br />price:   839<br />color: G","carat: 0.71<br />price:  2492<br />color: G","carat: 0.43<br />price:   729<br />color: G","carat: 0.82<br />price:  3652<br />color: G","carat: 0.31<br />price:   542<br />color: G","carat: 1.07<br />price:  7858<br />color: G","carat: 1.01<br />price:  4672<br />color: G","carat: 1.20<br />price:  5996<br />color: G","carat: 1.52<br />price: 12161<br />color: G","carat: 0.24<br />price:   559<br />color: G","carat: 0.35<br />price:   756<br />color: G","carat: 0.39<br />price:  1095<br />color: G","carat: 0.54<br />price:  1985<br />color: G","carat: 1.01<br />price:  6239<br />color: G","carat: 0.40<br />price:   765<br />color: G","carat: 0.43<br />price:  1264<br />color: G","carat: 0.32<br />price:   698<br />color: G","carat: 0.52<br />price:  2012<br />color: G","carat: 0.25<br />price:   460<br />color: G","carat: 0.55<br />price:  1667<br />color: G","carat: 1.71<br />price: 11227<br />color: G","carat: 0.33<br />price:  1114<br />color: G","carat: 0.35<br />price:   669<br />color: G","carat: 0.53<br />price:  1655<br />color: G","carat: 0.41<br />price:  1068<br />color: G","carat: 1.01<br />price:  4004<br />color: G","carat: 0.31<br />price:   583<br />color: G","carat: 0.31<br />price:   802<br />color: G","carat: 1.13<br />price:  4707<br />color: G","carat: 0.50<br />price:  1153<br />color: G","carat: 0.33<br />price:   868<br />color: G","carat: 0.50<br />price:  1354<br />color: G","carat: 0.32<br />price:   645<br />color: G","carat: 0.41<br />price:   986<br />color: G","carat: 0.71<br />price:  2649<br />color: G","carat: 0.51<br />price:  1438<br />color: G","carat: 1.00<br />price:  5483<br />color: G","carat: 0.31<br />price:   698<br />color: G","carat: 1.53<br />price: 13994<br />color: G","carat: 0.50<br />price:  1966<br />color: G","carat: 1.04<br />price:  9465<br />color: G","carat: 0.74<br />price:  3402<br />color: G","carat: 0.51<br />price:  1438<br />color: G","carat: 0.31<br />price:   816<br />color: G","carat: 0.31<br />price:   816<br />color: G","carat: 0.25<br />price:   454<br />color: G","carat: 1.00<br />price:  4095<br />color: G","carat: 0.71<br />price:  2515<br />color: G","carat: 0.90<br />price:  3303<br />color: G","carat: 1.70<br />price: 18718<br />color: G","carat: 0.76<br />price:  2972<br />color: G","carat: 0.40<br />price:   870<br />color: G","carat: 0.32<br />price:   720<br />color: G","carat: 0.52<br />price:  1244<br />color: G","carat: 2.32<br />price: 18532<br />color: G","carat: 0.44<br />price:   990<br />color: G","carat: 0.53<br />price:  1692<br />color: G","carat: 1.17<br />price:  9289<br />color: G","carat: 0.31<br />price:   789<br />color: G","carat: 1.01<br />price:  5741<br />color: G","carat: 0.55<br />price:  2016<br />color: G","carat: 0.32<br />price:   918<br />color: G","carat: 0.71<br />price:  2613<br />color: G","carat: 1.07<br />price:  4419<br />color: G","carat: 0.72<br />price:  3179<br />color: G"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(33,144,140,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(33,144,140,1)"}},"hoveron":"points","name":"G","legendgroup":"G","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1.76,0.9,0.54,1.19,1.02,0.75,2.08,0.57,1.11,1.5,0.9,0.41,1.1,1.56,1.2,0.42,0.72,1.21,0.86,0.43,2,0.35,0.83,0.53,1.3,1.19,0.31,0.53,0.32,1.41,1.5,0.57,1,0.41,1.52,1.1,0.3,0.36,0.43,0.32,0.41,0.7,0.35,0.7,0.33,1.73,1,0.44,0.53,1.41,1.29,1,0.7,2.77,0.8,1.51,0.35,1.45,0.9,0.32,1.28,0.3,0.4,0.36,1.51,1.21,0.62,0.9,0.32,1.37,0.37,0.42,0.71,0.3,0.42,0.51,1.28,0.34,0.63,0.3,0.35,0.32,1.23,2.01,1.21,0.41,0.71,0.9,2.02,1.2,0.5,1.16,1.02,0.4,1,0.7,1.2,0.77,0.7,0.32,1.01,1.18,0.32,1.56,1.09,0.3,0.31,2.01,1.71,0.37,1,0.9,1.71,1.63,0.5,1.17,1.09,0.34,1.06,1.01,0.3,0.44,0.77,1.5,1.5,0.54,0.33,0.41,0.31,0.51,0.32,0.96,0.85,1.7,0.52,0.26,0.36,0.34,0.71,0.9,1.01,0.3,1.21,2.19,1.71,2,0.92,0.97,1],"y":[13626,4194,1439,7147,5335,3011,15065,1077,4456,10584,4657,1115,4276,12098,5016,1087,2056,5936,3115,716,12554,647,2201,1020,5824,7152,435,2031,936,9752,9987,1558,5000,1043,9527,5729,378,794,919,561,1011,2508,788,2137,854,6007,3248,1124,2154,9254,6588,4702,2010,10424,2867,10763,798,8550,4234,580,3355,450,1017,475,10411,7426,1616,4082,708,6337,666,1108,2918,540,898,1443,6093,852,1196,655,859,505,5717,17458,6974,876,2246,3989,15996,6129,1844,6954,4703,800,5069,2066,5082,3996,2815,461,3999,5674,461,11922,4700,675,660,6592,9193,749,4732,4309,13376,12620,1068,4717,4639,526,4401,4965,878,1196,3208,10206,9576,1560,787,827,558,1438,936,4456,1274,14275,1639,482,980,477,2096,2934,4242,776,4249,17016,14893,12765,3854,4943,4042],"text":["carat: 1.76<br />price: 13626<br />color: H","carat: 0.90<br />price:  4194<br />color: H","carat: 0.54<br />price:  1439<br />color: H","carat: 1.19<br />price:  7147<br />color: H","carat: 1.02<br />price:  5335<br />color: H","carat: 0.75<br />price:  3011<br />color: H","carat: 2.08<br />price: 15065<br />color: H","carat: 0.57<br />price:  1077<br />color: H","carat: 1.11<br />price:  4456<br />color: H","carat: 1.50<br />price: 10584<br />color: H","carat: 0.90<br />price:  4657<br />color: H","carat: 0.41<br />price:  1115<br />color: H","carat: 1.10<br />price:  4276<br />color: H","carat: 1.56<br />price: 12098<br />color: H","carat: 1.20<br />price:  5016<br />color: H","carat: 0.42<br />price:  1087<br />color: H","carat: 0.72<br />price:  2056<br />color: H","carat: 1.21<br />price:  5936<br />color: H","carat: 0.86<br />price:  3115<br />color: H","carat: 0.43<br />price:   716<br />color: H","carat: 2.00<br />price: 12554<br />color: H","carat: 0.35<br />price:   647<br />color: H","carat: 0.83<br />price:  2201<br />color: H","carat: 0.53<br />price:  1020<br />color: H","carat: 1.30<br />price:  5824<br />color: H","carat: 1.19<br />price:  7152<br />color: H","carat: 0.31<br />price:   435<br />color: H","carat: 0.53<br />price:  2031<br />color: H","carat: 0.32<br />price:   936<br />color: H","carat: 1.41<br />price:  9752<br />color: H","carat: 1.50<br />price:  9987<br />color: H","carat: 0.57<br />price:  1558<br />color: H","carat: 1.00<br />price:  5000<br />color: H","carat: 0.41<br />price:  1043<br />color: H","carat: 1.52<br />price:  9527<br />color: H","carat: 1.10<br />price:  5729<br />color: H","carat: 0.30<br />price:   378<br />color: H","carat: 0.36<br />price:   794<br />color: H","carat: 0.43<br />price:   919<br />color: H","carat: 0.32<br />price:   561<br />color: H","carat: 0.41<br />price:  1011<br />color: H","carat: 0.70<br />price:  2508<br />color: H","carat: 0.35<br />price:   788<br />color: H","carat: 0.70<br />price:  2137<br />color: H","carat: 0.33<br />price:   854<br />color: H","carat: 1.73<br />price:  6007<br />color: H","carat: 1.00<br />price:  3248<br />color: H","carat: 0.44<br />price:  1124<br />color: H","carat: 0.53<br />price:  2154<br />color: H","carat: 1.41<br />price:  9254<br />color: H","carat: 1.29<br />price:  6588<br />color: H","carat: 1.00<br />price:  4702<br />color: H","carat: 0.70<br />price:  2010<br />color: H","carat: 2.77<br />price: 10424<br />color: H","carat: 0.80<br />price:  2867<br />color: H","carat: 1.51<br />price: 10763<br />color: H","carat: 0.35<br />price:   798<br />color: H","carat: 1.45<br />price:  8550<br />color: H","carat: 0.90<br />price:  4234<br />color: H","carat: 0.32<br />price:   580<br />color: H","carat: 1.28<br />price:  3355<br />color: H","carat: 0.30<br />price:   450<br />color: H","carat: 0.40<br />price:  1017<br />color: H","carat: 0.36<br />price:   475<br />color: H","carat: 1.51<br />price: 10411<br />color: H","carat: 1.21<br />price:  7426<br />color: H","carat: 0.62<br />price:  1616<br />color: H","carat: 0.90<br />price:  4082<br />color: H","carat: 0.32<br />price:   708<br />color: H","carat: 1.37<br />price:  6337<br />color: H","carat: 0.37<br />price:   666<br />color: H","carat: 0.42<br />price:  1108<br />color: H","carat: 0.71<br />price:  2918<br />color: H","carat: 0.30<br />price:   540<br />color: H","carat: 0.42<br />price:   898<br />color: H","carat: 0.51<br />price:  1443<br />color: H","carat: 1.28<br />price:  6093<br />color: H","carat: 0.34<br />price:   852<br />color: H","carat: 0.63<br />price:  1196<br />color: H","carat: 0.30<br />price:   655<br />color: H","carat: 0.35<br />price:   859<br />color: H","carat: 0.32<br />price:   505<br />color: H","carat: 1.23<br />price:  5717<br />color: H","carat: 2.01<br />price: 17458<br />color: H","carat: 1.21<br />price:  6974<br />color: H","carat: 0.41<br />price:   876<br />color: H","carat: 0.71<br />price:  2246<br />color: H","carat: 0.90<br />price:  3989<br />color: H","carat: 2.02<br />price: 15996<br />color: H","carat: 1.20<br />price:  6129<br />color: H","carat: 0.50<br />price:  1844<br />color: H","carat: 1.16<br />price:  6954<br />color: H","carat: 1.02<br />price:  4703<br />color: H","carat: 0.40<br />price:   800<br />color: H","carat: 1.00<br />price:  5069<br />color: H","carat: 0.70<br />price:  2066<br />color: H","carat: 1.20<br />price:  5082<br />color: H","carat: 0.77<br />price:  3996<br />color: H","carat: 0.70<br />price:  2815<br />color: H","carat: 0.32<br />price:   461<br />color: H","carat: 1.01<br />price:  3999<br />color: H","carat: 1.18<br />price:  5674<br />color: H","carat: 0.32<br />price:   461<br />color: H","carat: 1.56<br />price: 11922<br />color: H","carat: 1.09<br />price:  4700<br />color: H","carat: 0.30<br />price:   675<br />color: H","carat: 0.31<br />price:   660<br />color: H","carat: 2.01<br />price:  6592<br />color: H","carat: 1.71<br />price:  9193<br />color: H","carat: 0.37<br />price:   749<br />color: H","carat: 1.00<br />price:  4732<br />color: H","carat: 0.90<br />price:  4309<br />color: H","carat: 1.71<br />price: 13376<br />color: H","carat: 1.63<br />price: 12620<br />color: H","carat: 0.50<br />price:  1068<br />color: H","carat: 1.17<br />price:  4717<br />color: H","carat: 1.09<br />price:  4639<br />color: H","carat: 0.34<br />price:   526<br />color: H","carat: 1.06<br />price:  4401<br />color: H","carat: 1.01<br />price:  4965<br />color: H","carat: 0.30<br />price:   878<br />color: H","carat: 0.44<br />price:  1196<br />color: H","carat: 0.77<br />price:  3208<br />color: H","carat: 1.50<br />price: 10206<br />color: H","carat: 1.50<br />price:  9576<br />color: H","carat: 0.54<br />price:  1560<br />color: H","carat: 0.33<br />price:   787<br />color: H","carat: 0.41<br />price:   827<br />color: H","carat: 0.31<br />price:   558<br />color: H","carat: 0.51<br />price:  1438<br />color: H","carat: 0.32<br />price:   936<br />color: H","carat: 0.96<br />price:  4456<br />color: H","carat: 0.85<br />price:  1274<br />color: H","carat: 1.70<br />price: 14275<br />color: H","carat: 0.52<br />price:  1639<br />color: H","carat: 0.26<br />price:   482<br />color: H","carat: 0.36<br />price:   980<br />color: H","carat: 0.34<br />price:   477<br />color: H","carat: 0.71<br />price:  2096<br />color: H","carat: 0.90<br />price:  2934<br />color: H","carat: 1.01<br />price:  4242<br />color: H","carat: 0.30<br />price:   776<br />color: H","carat: 1.21<br />price:  4249<br />color: H","carat: 2.19<br />price: 17016<br />color: H","carat: 1.71<br />price: 14893<br />color: H","carat: 2.00<br />price: 12765<br />color: H","carat: 0.92<br />price:  3854<br />color: H","carat: 0.97<br />price:  4943<br />color: H","carat: 1.00<br />price:  4042<br />color: H"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(53,183,121,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(53,183,121,1)"}},"hoveron":"points","name":"H","legendgroup":"H","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1.61,1.21,1,1.5,1.5,0.57,0.76,1.04,1.7,1.27,1.61,1.11,2.02,0.4,0.7,1.06,1.22,1.01,0.34,2.11,0.51,0.71,2,0.73,1.01,1.5,2.01,1.27,0.88,2.01,0.26,1.28,0.72,0.7,0.71,0.71,1.83,1.33,0.27,0.72,1.18,0.72,1.21,1.02,0.5,0.7,0.27,1.21,0.47,0.32,1.56,1.01,1.21,1.52,1.01,1.21,0.71,1.61,1.33,0.37,0.33,0.3,1.22,1,1.58,2.04,0.55,0.79,2.01,1.26,1.55,2.04,0.37,0.4,0.4,2.48,0.38,0.7,1.02,1.02,1.04,1.5,1.42,0.33,1.5,0.9,0.9,0.57,1,1.07,1.26,0.77,0.3,0.34],"y":[9341,5124,3276,9702,6115,1180,2175,4013,11888,5249,10098,3183,16128,765,1724,4745,5368,4759,574,15282,1140,2218,16407,2351,4118,7912,15475,7312,3328,10711,445,6762,2148,2083,2400,2161,13711,5094,534,2650,5551,2441,5421,3656,1054,2332,432,7086,884,648,7326,4872,5866,7506,3299,5962,2376,8978,7325,754,608,405,4099,4538,9457,15543,1281,2539,15952,5207,8056,16874,666,849,667,16820,746,2335,3856,4478,3774,8920,7157,891,12148,3274,2948,1212,3360,4451,6566,2375,422,765],"text":["carat: 1.61<br />price:  9341<br />color: I","carat: 1.21<br />price:  5124<br />color: I","carat: 1.00<br />price:  3276<br />color: I","carat: 1.50<br />price:  9702<br />color: I","carat: 1.50<br />price:  6115<br />color: I","carat: 0.57<br />price:  1180<br />color: I","carat: 0.76<br />price:  2175<br />color: I","carat: 1.04<br />price:  4013<br />color: I","carat: 1.70<br />price: 11888<br />color: I","carat: 1.27<br />price:  5249<br />color: I","carat: 1.61<br />price: 10098<br />color: I","carat: 1.11<br />price:  3183<br />color: I","carat: 2.02<br />price: 16128<br />color: I","carat: 0.40<br />price:   765<br />color: I","carat: 0.70<br />price:  1724<br />color: I","carat: 1.06<br />price:  4745<br />color: I","carat: 1.22<br />price:  5368<br />color: I","carat: 1.01<br />price:  4759<br />color: I","carat: 0.34<br />price:   574<br />color: I","carat: 2.11<br />price: 15282<br />color: I","carat: 0.51<br />price:  1140<br />color: I","carat: 0.71<br />price:  2218<br />color: I","carat: 2.00<br />price: 16407<br />color: I","carat: 0.73<br />price:  2351<br />color: I","carat: 1.01<br />price:  4118<br />color: I","carat: 1.50<br />price:  7912<br />color: I","carat: 2.01<br />price: 15475<br />color: I","carat: 1.27<br />price:  7312<br />color: I","carat: 0.88<br />price:  3328<br />color: I","carat: 2.01<br />price: 10711<br />color: I","carat: 0.26<br />price:   445<br />color: I","carat: 1.28<br />price:  6762<br />color: I","carat: 0.72<br />price:  2148<br />color: I","carat: 0.70<br />price:  2083<br />color: I","carat: 0.71<br />price:  2400<br />color: I","carat: 0.71<br />price:  2161<br />color: I","carat: 1.83<br />price: 13711<br />color: I","carat: 1.33<br />price:  5094<br />color: I","carat: 0.27<br />price:   534<br />color: I","carat: 0.72<br />price:  2650<br />color: I","carat: 1.18<br />price:  5551<br />color: I","carat: 0.72<br />price:  2441<br />color: I","carat: 1.21<br />price:  5421<br />color: I","carat: 1.02<br />price:  3656<br />color: I","carat: 0.50<br />price:  1054<br />color: I","carat: 0.70<br />price:  2332<br />color: I","carat: 0.27<br />price:   432<br />color: I","carat: 1.21<br />price:  7086<br />color: I","carat: 0.47<br />price:   884<br />color: I","carat: 0.32<br />price:   648<br />color: I","carat: 1.56<br />price:  7326<br />color: I","carat: 1.01<br />price:  4872<br />color: I","carat: 1.21<br />price:  5866<br />color: I","carat: 1.52<br />price:  7506<br />color: I","carat: 1.01<br />price:  3299<br />color: I","carat: 1.21<br />price:  5962<br />color: I","carat: 0.71<br />price:  2376<br />color: I","carat: 1.61<br />price:  8978<br />color: I","carat: 1.33<br />price:  7325<br />color: I","carat: 0.37<br />price:   754<br />color: I","carat: 0.33<br />price:   608<br />color: I","carat: 0.30<br />price:   405<br />color: I","carat: 1.22<br />price:  4099<br />color: I","carat: 1.00<br />price:  4538<br />color: I","carat: 1.58<br />price:  9457<br />color: I","carat: 2.04<br />price: 15543<br />color: I","carat: 0.55<br />price:  1281<br />color: I","carat: 0.79<br />price:  2539<br />color: I","carat: 2.01<br />price: 15952<br />color: I","carat: 1.26<br />price:  5207<br />color: I","carat: 1.55<br />price:  8056<br />color: I","carat: 2.04<br />price: 16874<br />color: I","carat: 0.37<br />price:   666<br />color: I","carat: 0.40<br />price:   849<br />color: I","carat: 0.40<br />price:   667<br />color: I","carat: 2.48<br />price: 16820<br />color: I","carat: 0.38<br />price:   746<br />color: I","carat: 0.70<br />price:  2335<br />color: I","carat: 1.02<br />price:  3856<br />color: I","carat: 1.02<br />price:  4478<br />color: I","carat: 1.04<br />price:  3774<br />color: I","carat: 1.50<br />price:  8920<br />color: I","carat: 1.42<br />price:  7157<br />color: I","carat: 0.33<br />price:   891<br />color: I","carat: 1.50<br />price: 12148<br />color: I","carat: 0.90<br />price:  3274<br />color: I","carat: 0.90<br />price:  2948<br />color: I","carat: 0.57<br />price:  1212<br />color: I","carat: 1.00<br />price:  3360<br />color: I","carat: 1.07<br />price:  4451<br />color: I","carat: 1.26<br />price:  6566<br />color: I","carat: 0.77<br />price:  2375<br />color: I","carat: 0.30<br />price:   422<br />color: I","carat: 0.34<br />price:   765<br />color: I"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(143,215,68,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(143,215,68,1)"}},"hoveron":"points","name":"I","legendgroup":"I","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[1.38,1.05,0.51,0.85,1.02,0.61,4.5,0.74,2.03,1.6,0.3,2.01,1.11,1.12,1.65,1.5,1.07,0.4,0.76,0.77,1.35,1.22,1.74,0.34,2,1.6,1.25,1.07,0.4,2.01,0.31,0.72,1.08,1.25,1.55,0.35,1.21,1.37,1.17,1.25,0.56,2.31,1.01,0.53,0.41,2,1.5,1.27,0.9,2.01,1.76],"y":[6494,4145,1071,2295,3879,1233,18531,2154,12268,8904,394,12607,4497,4066,5914,8760,4686,641,1809,2086,4974,4676,8690,469,6521,7550,4690,4422,631,12209,453,2136,4633,5110,8749,591,6299,5063,2825,4331,1119,12369,3906,1122,637,11189,6750,5690,3024,9725,9053],"text":["carat: 1.38<br />price:  6494<br />color: J","carat: 1.05<br />price:  4145<br />color: J","carat: 0.51<br />price:  1071<br />color: J","carat: 0.85<br />price:  2295<br />color: J","carat: 1.02<br />price:  3879<br />color: J","carat: 0.61<br />price:  1233<br />color: J","carat: 4.50<br />price: 18531<br />color: J","carat: 0.74<br />price:  2154<br />color: J","carat: 2.03<br />price: 12268<br />color: J","carat: 1.60<br />price:  8904<br />color: J","carat: 0.30<br />price:   394<br />color: J","carat: 2.01<br />price: 12607<br />color: J","carat: 1.11<br />price:  4497<br />color: J","carat: 1.12<br />price:  4066<br />color: J","carat: 1.65<br />price:  5914<br />color: J","carat: 1.50<br />price:  8760<br />color: J","carat: 1.07<br />price:  4686<br />color: J","carat: 0.40<br />price:   641<br />color: J","carat: 0.76<br />price:  1809<br />color: J","carat: 0.77<br />price:  2086<br />color: J","carat: 1.35<br />price:  4974<br />color: J","carat: 1.22<br />price:  4676<br />color: J","carat: 1.74<br />price:  8690<br />color: J","carat: 0.34<br />price:   469<br />color: J","carat: 2.00<br />price:  6521<br />color: J","carat: 1.60<br />price:  7550<br />color: J","carat: 1.25<br />price:  4690<br />color: J","carat: 1.07<br />price:  4422<br />color: J","carat: 0.40<br />price:   631<br />color: J","carat: 2.01<br />price: 12209<br />color: J","carat: 0.31<br />price:   453<br />color: J","carat: 0.72<br />price:  2136<br />color: J","carat: 1.08<br />price:  4633<br />color: J","carat: 1.25<br />price:  5110<br />color: J","carat: 1.55<br />price:  8749<br />color: J","carat: 0.35<br />price:   591<br />color: J","carat: 1.21<br />price:  6299<br />color: J","carat: 1.37<br />price:  5063<br />color: J","carat: 1.17<br />price:  2825<br />color: J","carat: 1.25<br />price:  4331<br />color: J","carat: 0.56<br />price:  1119<br />color: J","carat: 2.31<br />price: 12369<br />color: J","carat: 1.01<br />price:  3906<br />color: J","carat: 0.53<br />price:  1122<br />color: J","carat: 0.41<br />price:   637<br />color: J","carat: 2.00<br />price: 11189<br />color: J","carat: 1.50<br />price:  6750<br />color: J","carat: 1.27<br />price:  5690<br />color: J","carat: 0.90<br />price:  3024<br />color: J","carat: 2.01<br />price:  9725<br />color: J","carat: 1.76<br />price:  9053<br />color: J"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(253,231,37,1)","opacity":1,"size":15.1181102362205,"symbol":"circle","line":{"width":1.88976377952756,"color":"rgba(253,231,37,1)"}},"hoveron":"points","name":"J","legendgroup":"J","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":40.1826484018265,"l":54.7945205479452},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[0.0165,4.7135],"tickmode":"array","ticktext":["1","2","3","4"],"tickvals":[1,2,3,4],"categoryorder":"array","categoryarray":["1","2","3","4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":{"text":"carat","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-555.8,19635.8],"tickmode":"array","ticktext":["0","5000","10000","15000"],"tickvals":[0,5000,10000,15000],"categoryorder":"array","categoryarray":["0","5000","10000","15000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":{"text":"price","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"y":0.96751968503937},"annotations":[{"text":"color","x":1.02,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"left","yanchor":"bottom","legendTitle":true}],"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","showSendToCloud":false},"source":"A","attrs":{"249d494e9865":{"x":{},"y":{},"colour":{},"type":"scatter"}},"cur_data":"249d494e9865","visdat":{"249d494e9865":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

#### Regression line


```r
p <- ggplot(dsmall, aes(carat, price)) + geom_point() + 
            geom_smooth(method="lm", se=FALSE) +
    	    theme(panel.background=element_rect(fill = "white", colour = "black"))
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_regression_line-1.png" width="672" />

#### Several regression lines


```r
p <- ggplot(dsmall, aes(carat, price, group=color)) + 
            geom_point(aes(color=color), size=2) + 
            geom_smooth(aes(color=color), method = "lm", se=FALSE) 
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_many_regression_lines-1.png" width="672" />

#### Local regression curve (loess)


```r
p <- ggplot(dsmall, aes(carat, price)) + geom_point() + geom_smooth() 
print(p) # Setting se=FALSE removes error shade
```

<img src="rgraphics_files/figure-html/ggplot_loess_regression-1.png" width="672" />

### `ggplot`: line plot


```r
p <- ggplot(iris, aes(Petal.Length, Petal.Width, group=Species, 
            color=Species)) + geom_line() 
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_line_plot-1.png" width="672" />

### Faceting


```r
p <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) + 
    	    geom_line(aes(color=Species), size=1) + 
            facet_wrap(~Species, ncol=1)
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_line_plot_faceting-1.png" width="672" />

#### Exercise 3

Scatter plots with `ggplot2`

- __Task 1__: Generate scatter plot for first two columns in `iris` data frame and color dots by its `Species` column.
- __Task 2__: Use the `xlim` and `ylim` arguments to set limits on the x- and y-axes so that all data points are restricted to the left bottom quadrant of the plot. 
- __Task 3__: Generate corresponding line plot with faceting presenting the individual data sets in saparate plots. 

Structure of `iris` data set


```r
class(iris)
```

```
## [1] "data.frame"
```

```r
iris[1:4,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
```

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```



### Bar Plots

Sample Set: the following transforms the `iris` data set into a ggplot2-friendly format.

Calculate mean values for aggregates given by `Species` column in `iris` data set


```r
iris_mean <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=mean) 
```

Calculate standard deviations for aggregates given by `Species` column in `iris` data set


```r
iris_sd <- aggregate(iris[,1:4], by=list(Species=iris$Species), FUN=sd) 
```

Reformat `iris_mean` with `melt`


```r
library(reshape2) # Defines melt function
df_mean <- melt(iris_mean, id.vars=c("Species"), variable.name = "Samples", value.name="Values")
```

Reformat `iris_sd` with `melt`


```r
df_sd <- melt(iris_sd, id.vars=c("Species"), variable.name = "Samples", value.name="Values")
```

Define standard deviation limits


```r
limits <- aes(ymax = df_mean[,"Values"] + df_sd[,"Values"], ymin=df_mean[,"Values"] - df_sd[,"Values"])
```

#### Verical orientation


```r
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
            geom_bar(position="dodge", stat="identity")
print(p) 
```

<img src="rgraphics_files/figure-html/iris_mean_bar_plot-1.png" width="672" />

To enforce that the bars are plotted in the order specified in the input data, one can instruct `ggplot` 
to do so by turning the corresponding column (here `Species`) into an ordered factor as follows. 



```r
df_mean$Species <- factor(df_mean$Species, levels=unique(df_mean$Species), ordered=TRUE)
```

In the above example this is not necessary since `ggplot` uses this order already. 

#### Horizontal orientation


```r
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
            geom_bar(position="dodge", stat="identity") + coord_flip() + 
            theme(axis.text.y=element_text(angle=0, hjust=1))
print(p) 
```

<img src="rgraphics_files/figure-html/iris_mean_bar_plot_sideways-1.png" width="672" />

#### Faceting


```r
p <- ggplot(df_mean, aes(Samples, Values)) + geom_bar(aes(fill = Species), stat="identity") +
            facet_wrap(~Species, ncol=1)
print(p) 
```

#### Error bars


```r
p <- ggplot(df_mean, aes(Samples, Values, fill = Species)) + 
            geom_bar(position="dodge", stat="identity") + 
            geom_errorbar(limits, position="dodge") 
print(p) 
```

<img src="rgraphics_files/figure-html/iris_mean_bar_plot_error-1.png" width="672" />

#### Mirrored 


```r
df <- data.frame(group = rep(c("Above", "Below"), each=10), x = rep(1:10, 2), y = c(runif(10, 0, 1), runif(10, -1, 0)))
p <- ggplot(df, aes(x=x, y=y, fill=group)) + 
	        geom_bar(stat="identity", position="identity")
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot2_mirrored_barplot-1.png" width="672" />

### Changing Color Settings


```r
library(RColorBrewer)
# display.brewer.all() 
p <- ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) +
            geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge") + 
            scale_fill_brewer(palette="Blues") + scale_color_brewer(palette = "Greys") 
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot2_color1-1.png" width="672" />

#### Using standard R color theme


```r
p <- ggplot(df_mean, aes(Samples, Values, fill=Species, color=Species)) + 
            geom_bar(position="dodge", stat="identity") + geom_errorbar(limits, position="dodge") + 
            scale_fill_manual(values=c("red", "green3", "blue")) + 
            scale_color_manual(values=c("red", "green3", "blue")) 
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot2_color2-1.png" width="672" />

#### Exercise 4

Bar plots

- __Task 1__: Calculate the mean values for the `Species` components of the first four columns in the `iris` data set. Use the `melt` function from the `reshape2` package to bring the data into the expected format for `ggplot`.
- __Task 2__: Generate two bar plots: one with stacked bars and one with horizontally arranged bars. 

Structure of iris data set


```r
class(iris)
```

```
## [1] "data.frame"
```

```r
iris[1:4,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
## 3          4.7         3.2          1.3         0.2  setosa
## 4          4.6         3.1          1.5         0.2  setosa
```

```r
table(iris$Species)
```

```
## 
##     setosa versicolor  virginica 
##         50         50         50
```



### Data reformatting example 

Here for line plot


```r
y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("Sample", 1:5, sep="")))
y <- data.frame(Position=1:length(y[,1]), y)
y[1:4, ] # First rows of input format expected by melt()
```

```
##    Position    Sample1    Sample2    Sample3    Sample4     Sample5
## g1        1  1.5336975 -1.0365027 -2.0276195 -0.4580396 -0.06460952
## g2        2 -2.0960304  2.1878704  0.7260334  0.8274617  0.24192162
## g3        3 -0.8233125  0.4250477  0.6526331 -0.4509962 -1.06778801
## g4        4  1.0961555  0.8101104 -0.3403762 -0.7222191 -0.72737741
```

```r
df <- melt(y, id.vars=c("Position"), variable.name = "Samples", value.name="Values")
p <- ggplot(df, aes(Position, Values)) + geom_line(aes(color=Samples)) + facet_wrap(~Samples, ncol=1)
print(p)
```

<img src="rgraphics_files/figure-html/ggplot_melt_data-1.png" width="672" />

Same data can be represented in box plot as follows


```r
ggplot(df, aes(Samples, Values, fill=Samples)) + geom_boxplot() + geom_jitter(color="darkgrey")
```

### Jitter Plots


```r
p <- ggplot(dsmall, aes(color, price/carat)) + 
            geom_jitter(alpha = I(1 / 2), aes(color=color))
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_jitter_plot-1.png" width="672" />

### Box plots


```r
p <- ggplot(dsmall, aes(color, price/carat, fill=color)) + geom_boxplot()
print(p) 
```

<img src="rgraphics_files/figure-html/ggplot_box_plot2-1.png" width="672" />

### Violin plots                                                                                                                                                               
                                                                                                                                                                              

```r
p <- ggplot(dsmall, aes(color, price/carat, fill=color)) + geom_violin()                                                                                                      
print(p)                                                                                                                                                                      
```

<img src="rgraphics_files/figure-html/ggplot_box_violin-1.png" width="672" />

Same violin plot as interactive plot generated with `ggplotly`, where the actual data points 
are shown as well by including `geom_jitter()`. 















































