---
title: Graphics and Data Visualization in R
author: "First/last name (first.last@ucr.edu)"
date: "Last update: 17 July, 2021" 
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
weight: 15
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

#### Bar plot of loan payments and amortization tables

The following imports a `mortgage` payment function (from [here](https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rgraphics/scripts/mortgage.R)) that calculates
monthly and annual mortgage/loan payments, generates amortization tables and
plots the results in form of a bar plot. A Shiny App using this function has been created
by Antoine Soetewey [here](https://antoinesoetewey.shinyapps.io/mortgage-calculator/). 


```r
source("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rgraphics/scripts/mortgage.R")
```

```
## The monthly mortgage payments and amortization rates can be calculted with the mortgage() function like this: 
##  
##     m <- mortgage(P=500000, I=6, L=30, plotData=TRUE)
##         P = principal (loan amount)
##         I = annual interest rate
##         L = length of the loan in years
```

```r
m <- mortgage(P=250000, I=6, L=15, plotData=TRUE)
```

```
## 
## The payments for this loan are:
##  
##             Monthly payment: $2109.642 (stored in m$monthPay)
## 
##             Total cost: $379735.6
## 
## The amortization data for each of the 180 months are stored in "m$aDFmonth".
## 
## The amortization data for each of the 15 years are stored in "m$aDFyear".
```

<img src="rgraphics_files/figure-html/bar_plot_mortgage-1.png" width="672" />

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

Much more on colors in R see Earl Glynn's [color chart](http://earlglynn.github.io/color/ColorInR)


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

### Design Concept of `ggplot2`

Plotting formalized and implemented by the grammar of graphics by Leland Wilkinson and Hadley Wickham [@Wickham2010-nu; @Wilkinson2012-wv; @Wickham2009-aq]. The plotting process
in `ggplot2` is devided into layers including:

1. Data: the actual data to be plotted
2. Aesthetics: the scales onto which the data will be mapped
3. Geometries: shapes used to represent data (_e.g._ bar or scatter plot)
4. Facets: row and column layout of sub-plots
5. Statistics: data models and summaries
6. Coordinates: the plotting space
7. Theme: styles to be used, such as fonts, backgrounds, etc.

<center><img title="sphm" src="../images/grammar_ggplot2.png" width="500" /></center>


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
    - `tidyr` and `reshape2`

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











































































































