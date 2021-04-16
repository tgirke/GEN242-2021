---
title: "Introduction to R" 
author: "Author: Thomas Girke"
date: "Last update: 16 April, 2021" 
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
weight: 4
type: docs
---

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('Rbasics.Rmd', c('html_document'), clean=FALSE); knitr::knit('Rbasics.Rmd', tangle=TRUE)"; Rscript ../md2jekyll.R Rbasics.knit.md 8; Rscript -e "rmarkdown::render('Rbasics.Rmd', c('pdf_document'))"
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

Source code downloads:    
\[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rbasics/Rbasics.Rmd) \]    
\[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rbasics/Rbasics.R) \]

</div>

## Overview

## What is R?

[R](http://cran.at.r-project.org) is a powerful statistical environment and
programming language for the analysis and visualization of data. The
associated [Bioconductor](http://bioconductor.org/) and CRAN package
repositories provide many additional R packages for statistical data analysis
for a wide array of research areas. The R software is free and runs on all
common operating systems.

## Why Using R?

-   Complete statistical environment and programming language
-   Efficient functions and data structures for data analysis
-   Powerful graphics
-   Access to fast growing number of analysis packages
-   Most widely used language in bioinformatics
-   Is standard for data mining and biostatistical analysis
-   Technical advantages: free, open-source, available for all OSs

## Books and Documentation

-   simpleR - Using R for Introductory Statistics (John Verzani, 2004) - [URL](http://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf)
-   Bioinformatics and Computational Biology Solutions Using R and Bioconductor (Gentleman et al., 2005) - [URL](http://www.bioconductor.org/help/publications/books/bioinformatics-and-computational-biology-solutions/)
-   More on this see “Finding Help” section in UCR Manual - [URL](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Finding-Help)

## R Working Environments

<center>
<img title="R_Interfaces" src="../images/rinterface.png"/>
</center>
<center>
R Projects and Interfaces
</center>

Some R working environments with support for syntax highlighting and utilities to send code
to the R console:

-   [RStudio](https://www.rstudio.com/products/rstudio/features): excellent choice for beginners ([Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf))
-   Basic R code editors provided by Rguis
-   [gedit](https://wiki.gnome.org/Apps/Gedit), [Rgedit](http://rgedit.sourceforge.net/), [RKWard](https://rkward.kde.org/), [Eclipse](http://www.walware.de/goto/statet), [Tinn-R](http://www.sciviews.org/Tinn-R/), [Notepad++](https://notepad-plus-plus.org/), [NppToR](http://sourceforge.net/projects/npptor/)
-   [Vim-R-Tmux](http://manuals.bioinformatics.ucr.edu/home/programming-in-r/vim-r): R working environment based on vim and tmux
-   [Emacs](http://www.xemacs.org/Download/index.html) ([ESS add-on package](http://ess.r-project.org/))

### Example: RStudio

New integrated development environment (IDE) for [R](http://www.rstudio.com/ide/download/). Highly functional for both beginners and
advanced.

<center>
<img title="RStudio" src="../images/rstudio.png"/>
</center>
<center>
RStudio IDE
</center>

Some userful shortcuts: `Ctrl+Enter` (send code), `Ctrl+Shift+C` (comment/uncomment), `Ctrl+1/2` (switch window focus)

### Example: Nvim-R-Tmux

Terminal-based Working Environment for R: [Nvim-R-Tmux](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/#reading-and-writing-external-data).

<center>
<img title="Nvim-R-Tmux" src="../images/Nvim-R.gif" >
</center>
<center>
Nvim-R-Tmux IDE for R
</center>

## R Package Repositories

-   CRAN (&gt;14,000 packages) general data analysis - [URL](http://cran.at.r-project.org/)
-   Bioconductor (&gt;2,000 packages) bioscience data analysis - [URL](http://www.bioconductor.org/)
-   Omegahat (&gt;90 packages) programming interfaces - [URL](https://github.com/omegahat?tab=repositories)
-   RStudio packages - [URL](https://www.rstudio.com/products/rpackages/)

## Installation of R, RStudio and R Packages

1.  Install R for your operating system from [CRAN](http://cran.at.r-project.org/).

2.  Install RStudio from [RStudio](http://www.rstudio.com/ide/download).

3.  Install CRAN Packages from R console like this:

    ``` r
    install.packages(c("pkg1", "pkg2")) 
    install.packages("pkg.zip", repos=NULL)
    ```

4.  Install Bioconductor packages as follows:

    ``` r
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager") # Installs BiocManager if not available yet
    BiocManager::version() # Reports Bioconductor version
    BiocManager::install(c("pkg1", "pkg2")) # Installs packages specified under "pkg1" 
    ```

5.  For more details consult the [Bioc Install page](http://www.bioconductor.org/install/)
    and [BiocInstaller](http://www.bioconductor.org/packages/release/bioc/html/BiocInstaller.html) package.

## Getting Around

### Startup and Closing Behavior

-   **Starting R**:
    The R GUI versions, including RStudio, under Windows and Mac OS X can be
    opened by double-clicking their icons. Alternatively, one can start it by
    typing `R` in a terminal (default under Linux).

-   **Startup/Closing Behavior**:
    The R environment is controlled by hidden files in the startup directory:
    `.RData`, `.Rhistory` and `.Rprofile` (optional).

-   **Closing R**:

``` r
q()  
```

    Save workspace image? [y/n/c]:

-   **Note**:
    When responding with `y`, then the entire R workspace will be written to
    the `.RData` file which can become very large. Often it is better to select `n` here,
    because a much better working pratice is to save an analysis protocol to an `R` or `Rmd` source file.
    This way one can quickly regenerate all data sets and objects needed in a future session.

## Navigating directories

List objects in current R session

``` r
ls()
```

Return content of current working directory

``` r
dir()
```

Return path of current working directory

``` r
getwd()
```

Change current working directory

``` r
setwd("/home/user")
```

## Basic Syntax

Create an object with the assignment operator `<-` or `=`

``` r
object <- ...
```

General R command syntax

``` r
object <- function_name(arguments) 
object <- object[arguments] 
```

Instead of the assignment operator one can use the `assign` function

``` r
assign("x", function(arguments))
```

Finding help

``` r
?function_name
```

Load a library/package

``` r
library("my_library") 
```

List functions defined by a library

``` r
library(help="my_library")
```

Load library manual (PDF or HTML file)

``` r
vignette("my_library") 
```

Execute an R script from within R

``` r
source("my_script.R")
```

Execute an R script from command-line (the first of the three options is preferred)

``` sh
$ Rscript my_script.R
$ R CMD BATCH my_script.R 
$ R --slave < my_script.R 
```

## Data Types

### Numeric data

Example: `1, 2, 3, ...`

``` r
x <- c(1, 2, 3)
x
```

    ## [1] 1 2 3

``` r
is.numeric(x)
```

    ## [1] TRUE

``` r
as.character(x)
```

    ## [1] "1" "2" "3"

### Character data

Example: `"a", "b", "c", ...`

``` r
x <- c("1", "2", "3")
x
```

    ## [1] "1" "2" "3"

``` r
is.character(x)
```

    ## [1] TRUE

``` r
as.numeric(x)
```

    ## [1] 1 2 3

### Complex data

Example: mix of both

``` r
c(1, "b", 3)
```

    ## [1] "1" "b" "3"

### Logical data

Example: `TRUE` of `FALSE`

``` r
x <- 1:10 < 5
x  
```

    ##  [1]  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

``` r
!x
```

    ##  [1] FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

``` r
which(x) # Returns index for the 'TRUE' values in logical vector
```

    ## [1] 1 2 3 4

## Data Objects

### Object types

#### Vectors (1D)

Definition: `numeric` or `character`

``` r
myVec <- 1:10; names(myVec) <- letters[1:10]  
myVec <- setNames(1:10, letters[1:10]) # Same as above in single step
myVec[1:5]
```

    ## a b c d e 
    ## 1 2 3 4 5

``` r
myVec[c(2,4,6,8)]
```

    ## b d f h 
    ## 2 4 6 8

``` r
myVec[c("b", "d", "f")]
```

    ## b d f 
    ## 2 4 6

#### Factors (1D)

Definition: vectors with grouping information

``` r
factor(c("dog", "cat", "mouse", "dog", "dog", "cat"))
```

    ## [1] dog   cat   mouse dog   dog   cat  
    ## Levels: cat dog mouse

#### Matrices (2D)

Definition: two dimensional structures with data of same type

``` r
myMA <- matrix(1:30, 3, 10, byrow = TRUE) 
class(myMA)
```

    ## [1] "matrix" "array"

``` r
myMA[1:2,]
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## [1,]    1    2    3    4    5    6    7    8    9    10
    ## [2,]   11   12   13   14   15   16   17   18   19    20

``` r
myMA[1, , drop=FALSE]
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ## [1,]    1    2    3    4    5    6    7    8    9    10

#### Data Frames (2D)

Definition: `data.frames` are two dimensional objects with data of variable types

``` r
myDF <- data.frame(Col1=1:10, Col2=10:1) 
myDF[1:2, ]
```

    ##   Col1 Col2
    ## 1    1   10
    ## 2    2    9

#### Tibbles

`Tibbles` are a more modern version of `data.frames`. Among many other advantages, one can see
here that `tibbles` have a nicer printing bahavior. Much more detailed information on this object
class is provided in the [`dplyr/tidyverse`](https://girke.bioinformatics.ucr.edu/GEN242/manuals/rbasics/rbasics/#dplyr-environment) manual section.

``` r
library(tidyverse)
as_tibble(iris)
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <fct>  
    ##  1          5.1         3.5          1.4         0.2 setosa 
    ##  2          4.9         3            1.4         0.2 setosa 
    ##  3          4.7         3.2          1.3         0.2 setosa 
    ##  4          4.6         3.1          1.5         0.2 setosa 
    ##  5          5           3.6          1.4         0.2 setosa 
    ##  6          5.4         3.9          1.7         0.4 setosa 
    ##  7          4.6         3.4          1.4         0.3 setosa 
    ##  8          5           3.4          1.5         0.2 setosa 
    ##  9          4.4         2.9          1.4         0.2 setosa 
    ## 10          4.9         3.1          1.5         0.1 setosa 
    ## # … with 140 more rows

#### Arrays

Definition: data structure with one, two or more dimensions

#### Lists

Definition: containers for any object type

``` r
myL <- list(name="Fred", wife="Mary", no.children=3, child.ages=c(4,7,9)) 
myL
```

    ## $name
    ## [1] "Fred"
    ## 
    ## $wife
    ## [1] "Mary"
    ## 
    ## $no.children
    ## [1] 3
    ## 
    ## $child.ages
    ## [1] 4 7 9

``` r
myL[[4]][1:2] 
```

    ## [1] 4 7

## Functions

Definition: piece of code

``` r
myfct <- function(arg1, arg2, ...) { 
    function_body 
}
```

## Subsetting of data objects

**(1.) Subsetting by positive or negative index/position numbers**

``` r
myVec <- 1:26; names(myVec) <- LETTERS 
myVec[1:4]
```

    ## A B C D 
    ## 1 2 3 4

**(2.) Subsetting by same length logical vectors**

``` r
myLog <- myVec > 10
myVec[myLog] 
```

    ##  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z 
    ## 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26

**(3.) Subsetting by field names**

``` r
myVec[c("B", "K", "M")]
```

    ##  B  K  M 
    ##  2 11 13

**(4.) Subset with `$` sign**: references a single column or list component by its name

``` r
iris$Species[1:8]
```

    ## [1] setosa setosa setosa setosa setosa setosa setosa setosa
    ## Levels: setosa versicolor virginica

## Important Utilities

### Combining Objects

The `c` function combines vectors and lists

``` r
c(1, 2, 3)
```

    ## [1] 1 2 3

``` r
x <- 1:3; y <- 101:103
c(x, y)
```

    ## [1]   1   2   3 101 102 103

The `cbind` and `rbind` functions can be used to append columns and rows, respecively.

``` r
ma <- cbind(x, y)
ma
```

    ##      x   y
    ## [1,] 1 101
    ## [2,] 2 102
    ## [3,] 3 103

``` r
rbind(ma, ma)
```

    ##      x   y
    ## [1,] 1 101
    ## [2,] 2 102
    ## [3,] 3 103
    ## [4,] 1 101
    ## [5,] 2 102
    ## [6,] 3 103

### Accessing Dimensions of Objects

Length and dimension information of objects

``` r
length(iris$Species)
```

    ## [1] 150

``` r
dim(iris)
```

    ## [1] 150   5

### Accessing Name Slots of Objects

Accessing row and column names of 2D objects

``` r
rownames(iris)[1:8]
```

    ## [1] "1" "2" "3" "4" "5" "6" "7" "8"

``` r
colnames(iris)
```

    ## [1] "Sepal.Length" "Sepal.Width"  "Petal.Length" "Petal.Width"  "Species"

Return name field of vectors and lists

``` r
names(myVec)
```

    ##  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X"
    ## [25] "Y" "Z"

``` r
names(myL)
```

    ## [1] "name"        "wife"        "no.children" "child.ages"

### Sorting Objects

The function `sort` returns a vector in ascending or descending order

``` r
sort(10:1)
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

The function `order` returns a sorting index for sorting an object

``` r
sortindex <- order(iris[,1], decreasing = FALSE)
sortindex[1:12]
```

    ##  [1] 14  9 39 43 42  4  7 23 48  3 30 12

``` r
iris[sortindex,][1:2,]
```

    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ## 14          4.3         3.0          1.1         0.1  setosa
    ## 9           4.4         2.9          1.4         0.2  setosa

``` r
sortindex <- order(-iris[,1]) # Same as decreasing=TRUE
```

Sorting multiple columns

``` r
iris[order(iris$Sepal.Length, iris$Sepal.Width),][1:2,]
```

    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ## 14          4.3         3.0          1.1         0.1  setosa
    ## 9           4.4         2.9          1.4         0.2  setosa

## Operators and Calculations

### Comparison Operators

Comparison operators: `==`, `!=`, `<`, `>`, `<=`, `>=`

``` r
1==1
```

    ## [1] TRUE

Logical operators: AND: `&`, OR: `|`, NOT: `!`

``` r
x <- 1:10; y <- 10:1
x > y & x > 5
```

    ##  [1] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE

### Basic Calculations

To look up math functions, see Function Index [here](http://cran.at.r-project.org/doc/manuals/R-intro.html#Function-and-variable-index)

``` r
x + y
```

    ##  [1] 11 11 11 11 11 11 11 11 11 11

``` r
sum(x)
```

    ## [1] 55

``` r
mean(x)
```

    ## [1] 5.5

``` r
apply(iris[1:6,1:3], 1, mean) 
```

    ##        1        2        3        4        5        6 
    ## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667

## Reading and Writing External Data

### Import of tabular data

Import of a tab-delimited tabular file

``` r
myDF <- read.delim("myData.xls", sep="\t")
```

Import of Google Sheets. The following example imports a sample Google Sheet from [here](https://docs.google.com/spreadsheets/d/1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM/edit#gid=472150521).
Detailed instructions for interacting from R with Google Sheets with the required `googlesheets4` package are [here](https://googlesheets4.tidyverse.org/).

``` r
library(googlesheets4)
mysheet <- read_sheet("1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM", skip=4)
myDF <- as.data.frame(mysheet)
myDF
```

Import from Excel sheets works well with `readxl`. For details see the `readxl` package manual [here](https://readxl.tidyverse.org/). Note: working with tab- or comma-delimited files is more flexible and highly preferred for automated analysis workflows.

``` r
library("readxl")
mysheet <- read_excel(targets_path, sheet="Sheet1")
```

Additional import functions are described in the `readr` package section [here](https://girke.bioinformatics.ucr.edu/GEN242/manuals/rbasics/rbasics/#reading-and-writing-tabular-files).

### Export of tabular data

``` r
write.table(myDF, file="myfile.xls", sep="\t", quote=FALSE, col.names=NA)
```

### Line-wise import

``` r
myDF <- readLines("myData.txt")
```

### Line-wise export

``` r
writeLines(month.name, "myData.txt")
```

### Export R object

``` r
mylist <- list(C1=iris[,1], C2=iris[,2]) # Example to export
saveRDS(mylist, "mylist.rds")
```

### Import R object

``` r
mylist <- readRDS("mylist.rds")
```

### Copy and paste into R

On Windows/Linux systems

``` r
read.delim("clipboard") 
```

On Mac OS X systems

``` r
read.delim(pipe("pbpaste")) 
```

### Copy and paste from R

On Windows/Linux systems

``` r
write.table(iris, "clipboard", sep="\t", col.names=NA, quote=FALSE) 
```

On Mac OS X systems

``` r
zz <- pipe('pbcopy', 'w')
write.table(iris, zz, sep="\t", col.names=NA, quote=FALSE)
close(zz) 
```

### Homework 3A

Homework 3A: [Object Subsetting Routines and Import/Export](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw03/hw03/)

## Useful R Functions

### Unique entries

Make vector entries unique with `unique`

``` r
length(iris$Sepal.Length)
```

    ## [1] 150

``` r
length(unique(iris$Sepal.Length))
```

    ## [1] 35

### Count occurrences

Count occurrences of entries with `table`

``` r
table(iris$Species)
```

    ## 
    ##     setosa versicolor  virginica 
    ##         50         50         50

### Aggregate data

Compute aggregate statistics with `aggregate`

``` r
aggregate(iris[,1:4], by=list(iris$Species), FUN=mean, na.rm=TRUE)
```

    ##      Group.1 Sepal.Length Sepal.Width Petal.Length Petal.Width
    ## 1     setosa        5.006       3.428        1.462       0.246
    ## 2 versicolor        5.936       2.770        4.260       1.326
    ## 3  virginica        6.588       2.974        5.552       2.026

### Intersect data

Compute intersect between two vectors with `%in%`

``` r
month.name %in% c("May", "July")
```

    ##  [1] FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE

### Merge data frames

Join two data frames by common field entries with `merge` (here row names `by.x=0`). To obtain only the common rows, change `all=TRUE` to `all=FALSE`. To merge on specific columns, refer to them by their position numbers or their column names.

``` r
frame1 <- iris[sample(1:length(iris[,1]), 30), ]
frame1[1:2,]
```

    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
    ## 55          6.5         2.8          4.6         1.5 versicolor
    ## 49          5.3         3.7          1.5         0.2     setosa

``` r
dim(frame1)
```

    ## [1] 30  5

``` r
my_result <- merge(frame1, iris, by.x = 0, by.y = 0, all = TRUE)
dim(my_result)
```

    ## [1] 150  11

## Graphics in R

### Advantages

-   Powerful environment for visualizing scientific data
-   Integrated graphics and statistics infrastructure
-   Publication quality graphics
-   Fully programmable
-   Highly reproducible
-   Full [LaTeX](http://www.latex-project.org/) and Markdown support via `knitr` and `R markdown`
-   Vast number of R packages with graphics utilities

### Documentation for R Graphics

**General**

-   Graphics Task Page - [URL](http://cran.r-project.org/web/views/Graphics.html)
-   R Graph Gallery - [URL](http://addictedtor.free.fr/graphiques/allgraph.php)
-   R Graphical Manual - [URL](http://cged.genes.nig.ac.jp/RGM2/index.php)
-   Paul Murrell’s book R (Grid) Graphics - [URL](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)

**Interactive graphics**

-   rggobi\` (GGobi) - [URL](http://www.ggobi.org/)
-   `iplots` - [URL](http://www.rosuda.org/iplots/)
-   Open GL (`rgl`) - [URL](http://rgl.neoscientists.org/gallery.shtml)

### Graphics Environments

**Viewing and saving graphics in R**

-   On-screen graphics
-   postscript, pdf, svg
-   jpeg, png, wmf, tiff, …

**Four major graphic environments**

1.  Low-level infrastructure

-   R Base Graphics (low- and high-level)
-   `grid`: [Manual](http://www.stat.auckland.ac.nz/~paul/grid/grid.html)

2.  High-level infrastructure
    \\begin{itemize}

-   `lattice`: [Manual](http://lmdvr.r-forge.r-project.org), [Intro](http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm), [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
-   `ggplot2`: [Manual](http://had.co.nz/ggplot2/), [Intro](http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html), [Book](http://had.co.nz/ggplot2/book/)

### Base Graphics: Overview

**Important high-level plotting functions**

-   `plot`: generic x-y plotting
-   `barplot`: bar plots
-   `boxplot`: box-and-whisker plot
-   `hist`: histograms
-   `pie`: pie charts
-   `dotchart`: cleveland dot plots
-   `image, heatmap, contour, persp`: functions to generate image-like plots
-   `qqnorm, qqline, qqplot`: distribution comparison plots
-   `pairs, coplot`: display of multivariant data

**Help on graphics functions**

-   `?myfct`
-   `?plot`
-   `?par`

#### Preferred Object Types

-   Matrices and data frames
-   Vectors
-   Named vectors

### Scatter Plots

### Basic Scatter Plot

Sample data set for subsequent plots

``` r
set.seed(1410)
y <- matrix(runif(30), ncol=3, dimnames=list(letters[1:10], LETTERS[1:3]))
```

Plot data

``` r
plot(y[,1], y[,2]) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/basic_scatter_plot-1.png" width="672" />

#### All pairs

``` r
pairs(y) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/pairs_scatter_plot-1.png" width="672" />

#### With labels

``` r
plot(y[,1], y[,2], pch=20, col="red", main="Symbols and Labels")
text(y[,1]+0.03, y[,2], rownames(y))
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/labels_scatter_plot-1.png" width="672" />

### More examples

**Print instead of symbols the row names**

``` r
plot(y[,1], y[,2], type="n", main="Plot of Labels")
text(y[,1], y[,2], rownames(y)) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/row_scatter_plot-1.png" width="672" />

**Usage of important plotting parameters**

``` r
grid(5, 5, lwd = 2) 
op <- par(mar=c(8,8,8,8), bg="lightblue")
plot(y[,1], y[,2], type="p", col="red", cex.lab=1.2, cex.axis=1.2, 
     cex.main=1.2, cex.sub=1, lwd=4, pch=20, xlab="x label", 
     ylab="y label", main="My Main", sub="My Sub")
par(op)
```

**Important arguments**

-   `mar`: specifies the margin sizes around the plotting area in order: `c(bottom, left, top, right)`
-   `col`: color of symbols
-   `pch`: type of symbols, samples: `example(points)`
-   `lwd`: size of symbols
-   `cex.*`: control font sizes
-   For details see `?par`

#### Add regression line

``` r
plot(y[,1], y[,2])
myline <- lm(y[,2]~y[,1]); abline(myline, lwd=2) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_regression-1.png" width="672" />

``` r
summary(myline) 
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
    ## Multiple R-squared:  0.09589,    Adjusted R-squared:  -0.01712 
    ## F-statistic: 0.8485 on 1 and 8 DF,  p-value: 0.3839

#### Log scale

Same plot as above, but on log scale

``` r
plot(y[,1], y[,2], log="xy") 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_regression_log-1.png" width="672" />

#### Add a mathematical expression

``` r
plot(y[,1], y[,2]); text(y[1,1], y[1,2], expression(sum(frac(1,sqrt(x^2*pi)))), cex=1.3) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_regression_math-1.png" width="672" />

### Homework 3B

Homework 3B: [Scatter Plots](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw03/hw03/)

### Line Plots

#### Single data set

``` r
plot(y[,1], type="l", lwd=2, col="blue") 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_line_single-1.png" width="672" />

#### Many Data Sets

Plots line graph for all columns in data frame `y`. The `split.screen` function is used in this example in a for loop to overlay several line graphs in the same plot.

``` r
split.screen(c(1,1)) 
```

    ## [1] 1

``` r
plot(y[,1], ylim=c(0,1), xlab="Measurement", ylab="Intensity", type="l", lwd=2, col=1)
for(i in 2:length(y[1,])) { 
    screen(1, new=FALSE)
    plot(y[,i], ylim=c(0,1), type="l", lwd=2, col=i, xaxt="n", yaxt="n", ylab="", xlab="", main="", bty="n") 
}
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_line_many-1.png" width="672" />

``` r
close.screen(all=TRUE) 
```

### Bar Plots

#### Basics

``` r
barplot(y[1:4,], ylim=c(0, max(y[1:4,])+0.3), beside=TRUE, legend=letters[1:4]) 
text(labels=round(as.vector(as.matrix(y[1:4,])),2), x=seq(1.5, 13, by=1) + sort(rep(c(0,1,2), 4)), y=as.vector(as.matrix(y[1:4,]))+0.04) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_bar_simple-1.png" width="672" />

#### Error Bars

``` r
bar <- barplot(m <- rowMeans(y) * 10, ylim=c(0, 10))
stdev <- sd(t(y))
arrows(bar, m, bar, m + stdev, length=0.15, angle = 90)
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_bar_error-1.png" width="672" />

### Histograms

``` r
hist(y, freq=TRUE, breaks=10)
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_hist-1.png" width="672" />

### Density Plots

``` r
plot(density(y), col="red")
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_dens-1.png" width="672" />

### Pie Charts

``` r
pie(y[,1], col=rainbow(length(y[,1]), start=0.1, end=0.8), clockwise=TRUE)
legend("topright", legend=row.names(y), cex=1.3, bty="n", pch=15, pt.cex=1.8, 
col=rainbow(length(y[,1]), start=0.1, end=0.8), ncol=1) 
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_pie-1.png" width="672" />

### Color Selection Utilities

Default color palette and how to change it

``` r
palette()
```

    ## [1] "black"   "#DF536B" "#61D04F" "#2297E6" "#28E2E5" "#CD0BBC" "#F5C710" "gray62"

``` r
palette(rainbow(5, start=0.1, end=0.2))
palette()
```

    ## [1] "#FF9900" "#FFBF00" "#FFE600" "#F2FF00" "#CCFF00"

``` r
palette("default")
```

The `gray` function allows to select any type of gray shades by providing values from 0 to 1

``` r
gray(seq(0.1, 1, by= 0.2))
```

    ## [1] "#1A1A1A" "#4D4D4D" "#808080" "#B3B3B3" "#E6E6E6"

Color gradients with `colorpanel` function from `gplots` library\`

``` r
library(gplots)
colorpanel(5, "darkblue", "yellow", "white")
```

    ## [1] "#00008B" "#808046" "#FFFF00" "#FFFF80" "#FFFFFF"

Much more on colors in R see Earl Glynn’s color chart [here](http://research.stowers-institute.org/efg/R/Color/Chart/)

### Saving Graphics to File

After the `pdf()` command all graphs are redirected to file `test.pdf`. Works for all common formats similarly: jpeg, png, ps, tiff, …

``` r
pdf("test.pdf")
plot(1:10, 1:10)
dev.off() 
```

Generates Scalable Vector Graphics (SVG) files that can be edited in vector graphics programs, such as InkScape.

``` r
library("RSvgDevice")
devSVG("test.svg")
plot(1:10, 1:10)
dev.off() 
```

### Homework 3C

Homework 3C: [Bar Plots](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw03/hw03/)

## Analysis Routine

### Overview

The following exercise introduces a variety of useful data analysis utilities in R.

### Analysis Routine: Data Import

-   **Step 1**: To get started with this exercise, direct your R session to a dedicated workshop directory and download into this directory the following sample tables. Then import the files into Excel and save them as tab delimited text files.

    -   [MolecularWeight\_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls)
    -   [TargetP\_analysis\_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls)

**Import the tables into R**

Import molecular weight table

``` r
my_mw <- read.delim(file="MolecularWeight_tair7.xls", header=TRUE, sep="\t") 
my_mw[1:2,]
```

Import subcelluar targeting table

``` r
my_target <- read.delim(file="TargetP_analysis_tair7.xls", header=TRUE, sep="\t") 
my_target[1:2,]
```

Online import of molecular weight table

``` r
my_mw <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls", header=TRUE, sep="\t") 
my_mw[1:2,]
```

    ##   Sequence.id Molecular.Weight.Da. Residues
    ## 1 AT1G08520.1                83285      760
    ## 2 AT1G08530.1                27015      257

Online import of subcelluar targeting table

``` r
my_target <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls", header=TRUE, sep="\t") 
my_target[1:2,]
```

    ##      GeneName Loc   cTP   mTP    SP other
    ## 1 AT1G08520.1   C 0.822 0.137 0.029 0.039
    ## 2 AT1G08530.1   C 0.817 0.058 0.010 0.100

### Merging Data Frames

-   **Step 2**: Assign uniform gene ID column titles

``` r
colnames(my_target)[1] <- "ID"
colnames(my_mw)[1] <- "ID" 
```

-   **Step 3**: Merge the two tables based on common ID field

``` r
my_mw_target <- merge(my_mw, my_target, by.x="ID", by.y="ID", all.x=T)
```

-   **Step 4**: Shorten one table before the merge and then remove the non-matching rows (NAs) in the merged file

``` r
my_mw_target2a <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all.x=T)  # To remove non-matching rows, use the argument setting 'all=F'.
my_mw_target2 <- na.omit(my_mw_target2a) # Removes rows containing "NAs" (non-matching rows).
```

-   **Homework 3D**: How can the merge function in the previous step be executed so that only the common rows among the two data frames are returned? Prove that both methods - the two step version with `na.omit` and your method - return identical results.
-   **Homework 3E**: Replace all `NAs` in the data frame `my_mw_target2a` with zeros.

### Filtering Data

-   **Step 5**: Retrieve all records with a value of greater than 100,000 in ‘MW’ column and ‘C’ value in ‘Loc’ column (targeted to chloroplast).

``` r
query <- my_mw_target[my_mw_target[, 2] > 100000 & my_mw_target[, 4] == "C", ] 
query[1:4, ]
```

    ##              ID Molecular.Weight.Da. Residues Loc   cTP   mTP    SP other
    ## 219 AT1G02730.1               132588     1181   C 0.972 0.038 0.008 0.045
    ## 243 AT1G02890.1               136825     1252   C 0.748 0.529 0.011 0.013
    ## 281 AT1G03160.1               100732      912   C 0.871 0.235 0.011 0.007
    ## 547 AT1G05380.1               126360     1138   C 0.740 0.099 0.016 0.358

``` r
dim(query)
```

    ## [1] 170   8

-   **Homework 3F**: How many protein entries in the `my`\_mw`_target` data frame have a MW of greater then 4,000 and less then 5,000. Subset the data frame accordingly and sort it by MW to check that your result is correct.

### String Substitutions

-   **Step 6**: Use a regular expression in a substitute function to generate a separate ID column that lacks the gene model extensions.

``` r
my_mw_target3 <- data.frame(loci=gsub("\\..*", "", as.character(my_mw_target[,1]), perl = TRUE), my_mw_target)
my_mw_target3[1:3,1:8]
```

    ##        loci          ID Molecular.Weight.Da. Residues Loc  cTP   mTP    SP
    ## 1 AT1G01010 AT1G01010.1                49426      429   _ 0.10 0.090 0.075
    ## 2 AT1G01020 AT1G01020.1                28092      245   * 0.01 0.636 0.158
    ## 3 AT1G01020 AT1G01020.2                21711      191   * 0.01 0.636 0.158

-   **Homework 3G**: Retrieve those rows in `my_mw_target3` where the second column contains the following identifiers: `c("AT5G52930.1", "AT4G18950.1", "AT1G15385.1", "AT4G36500.1", "AT1G67530.1")`. Use the `%in%` function for this query. As an alternative approach, assign the second column to the row index of the data frame and then perform the same query again using the row index. Explain the difference of the two methods.

### Calculations on Data Frames

-   **Step 7**: Count the number of duplicates in the loci column with the `table` function and append the result to the data frame with the `cbind` function.

``` r
mycounts <- table(my_mw_target3[,1])[my_mw_target3[,1]]
my_mw_target4 <- cbind(my_mw_target3, Freq=mycounts[as.character(my_mw_target3[,1])]) 
```

-   **Step 8**: Perform a vectorized devision of columns 3 and 4 (average AA weight per protein)

``` r
data.frame(my_mw_target4, avg_AA_WT=(my_mw_target4[,3] / my_mw_target4[,4]))[1:2,] 
```

    ##        loci          ID Molecular.Weight.Da. Residues Loc  cTP   mTP    SP other Freq.Var1
    ## 1 AT1G01010 AT1G01010.1                49426      429   _ 0.10 0.090 0.075 0.925 AT1G01010
    ## 2 AT1G01020 AT1G01020.1                28092      245   * 0.01 0.636 0.158 0.448 AT1G01020
    ##   Freq.Freq avg_AA_WT
    ## 1         1  115.2121
    ## 2         2  114.6612

-   **Step 9**: Calculate for each row the mean and standard deviation across several columns

``` r
mymean <- apply(my_mw_target4[,6:9], 1, mean)
mystdev <- apply(my_mw_target4[,6:9], 1, sd, na.rm=TRUE)
data.frame(my_mw_target4, mean=mymean, stdev=mystdev)[1:2,5:12] 
```

    ##   Loc  cTP   mTP    SP other Freq.Var1 Freq.Freq   mean
    ## 1   _ 0.10 0.090 0.075 0.925 AT1G01010         1 0.2975
    ## 2   * 0.01 0.636 0.158 0.448 AT1G01020         2 0.3130

### Plotting Example

-   **Step 10**: Generate scatter plot columns: ‘MW’ and ‘Residues’

``` r
plot(my_mw_target4[1:500,3:4], col="red")
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_example-1.png" width="672" />

### Export Results and Run Entire Exercise as Script

-   **Step 11**: Write the data frame `my_mw_target4` into a tab-delimited text file and inspect it in Excel.

``` r
write.table(my_mw_target4, file="my_file.xls", quote=FALSE, sep="\t", col.names = NA) 
```

-   **Homework 3H**: Write all commands from this exercise into an R script named `exerciseRbasics.R`, or download it from [here](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/exerciseRbasics.R). Then execute the script with the `source` function like this: `source("exerciseRbasics.R")`. This will run all commands of this exercise and generate the corresponding output files in the current working directory.

``` r
source("exerciseRbasics.R")
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.0.3 (2020-10-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Debian GNU/Linux 10 (buster)
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
    ##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
    ## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.30        gplots_3.1.1      RSQLite_2.2.1     data.table_1.13.2 forcats_0.5.0    
    ##  [6] stringr_1.4.0     dplyr_1.0.2       purrr_0.3.4       readr_1.4.0       tidyr_1.1.2      
    ## [11] tibble_3.0.4      tidyverse_1.3.0   ggplot2_3.3.2     limma_3.46.0      BiocStyle_2.18.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.2          viridisLite_0.3.0   bit64_4.0.5         jsonlite_1.7.1     
    ##  [5] modelr_0.1.8        gtools_3.8.2        assertthat_0.2.1    BiocManager_1.30.10
    ##  [9] highr_0.8           blob_1.2.1          cellranger_1.1.0    yaml_2.2.1         
    ## [13] pillar_1.4.7        backports_1.2.0     glue_1.4.2          digest_0.6.27      
    ## [17] rvest_0.3.6         colorspace_2.0-0    htmltools_0.5.1.1   plyr_1.8.6         
    ## [21] pkgconfig_2.0.3     broom_0.7.2         haven_2.3.1         bookdown_0.21      
    ## [25] scales_1.1.1        generics_0.1.0      farver_2.0.3        ellipsis_0.3.1     
    ## [29] withr_2.3.0         cli_2.2.0           magrittr_2.0.1      crayon_1.3.4       
    ## [33] readxl_1.3.1        memoise_1.1.0       evaluate_0.14       ps_1.4.0           
    ## [37] fs_1.5.0            fansi_0.4.1         xml2_1.3.2          blogdown_1.1.7     
    ## [41] tools_4.0.3         hms_0.5.3           lifecycle_0.2.0     munsell_0.5.0      
    ## [45] reprex_0.3.0        compiler_4.0.3      caTools_1.18.1      rlang_0.4.8        
    ## [49] grid_4.0.3          rstudioapi_0.13     bitops_1.0-6        labeling_0.4.2     
    ## [53] rmarkdown_2.5       gtable_0.3.0        codetools_0.2-16    DBI_1.1.0          
    ## [57] reshape2_1.4.4      R6_2.5.0            lubridate_1.7.9.2   bit_4.0.4          
    ## [61] utf8_1.1.4          KernSmooth_2.23-17  stringi_1.5.3       Rcpp_1.0.5         
    ## [65] vctrs_0.3.5         dbplyr_2.0.0        tidyselect_1.1.0    xfun_0.20

## References
