---
title: "Introduction to R" 
author: "Author: Thomas Girke"
date: "Last update: 07 April, 2021" 
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

<script src="/rmarkdown-libs/htmlwidgets/htmlwidgets.js"></script>

<script src="/rmarkdown-libs/jquery/jquery.min.js"></script>

<link href="/rmarkdown-libs/datatables-css/datatables-crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/datatables-binding/datatables.js"></script>

<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.min.css" rel="stylesheet" />
<link href="/rmarkdown-libs/dt-core/css/jquery.dataTables.extra.css" rel="stylesheet" />

<script src="/rmarkdown-libs/dt-core/js/jquery.dataTables.min.js"></script>

<link href="/rmarkdown-libs/nouislider/jquery.nouislider.min.css" rel="stylesheet" />

<script src="/rmarkdown-libs/nouislider/jquery.nouislider.min.js"></script>

<link href="/rmarkdown-libs/selectize/selectize.bootstrap3.css" rel="stylesheet" />

<script src="/rmarkdown-libs/selectize/selectize.min.js"></script>

<link href="/rmarkdown-libs/crosstalk/css/crosstalk.css" rel="stylesheet" />

<script src="/rmarkdown-libs/crosstalk/js/crosstalk.min.js"></script>

<!---
- Compile from command-line
Rscript -e "rmarkdown::render('Rbasics.Rmd', c('html_document'), clean=F); knitr::knit('Rbasics.Rmd', tangle=TRUE)"; Rscript ../md2jekyll.R Rbasics.knit.md 8; Rscript -e "rmarkdown::render('Rbasics.Rmd', c('pdf_document'))"
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
\[ [.Rmd](https://girke.bioinformatics.ucr.edu/GEN242/custom/custom/source_code/Rbasics.Rmd) \],
\[ [.R](https://girke.bioinformatics.ucr.edu/GEN242/custom/custom/source_code/Rbasics.Rmd) \]

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

  - Complete statistical environment and programming language
  - Efficient functions and data structures for data analysis
  - Powerful graphics
  - Access to fast growing number of analysis packages
  - Most widely used language in bioinformatics
  - Is standard for data mining and biostatistical analysis
  - Technical advantages: free, open-source, available for all OSs

## Books and Documentation

  - simpleR - Using R for Introductory Statistics (John Verzani, 2004) - [URL](http://cran.r-project.org/doc/contrib/Verzani-SimpleR.pdf)
  - Bioinformatics and Computational Biology Solutions Using R and Bioconductor (Gentleman et al., 2005) - [URL](http://www.bioconductor.org/help/publications/books/bioinformatics-and-computational-biology-solutions/)
  - More on this see “Finding Help” section in UCR Manual - [URL](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Finding-Help)

## R Working Environments

<center>

<img title="R_Interfaces" src="../images/rinterface.png"/>

</center>

<center>

R Projects and Interfaces

</center>

Some R working environments with support for syntax highlighting and utilities to send code
to the R console:

  - [RStudio](https://www.rstudio.com/products/rstudio/features): excellent choice for beginners ([Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2016/01/rstudio-IDE-cheatsheet.pdf))
  - Basic R code editors provided by Rguis
  - [gedit](https://wiki.gnome.org/Apps/Gedit), [Rgedit](http://rgedit.sourceforge.net/), [RKWard](https://rkward.kde.org/), [Eclipse](http://www.walware.de/goto/statet), [Tinn-R](http://www.sciviews.org/Tinn-R/), [Notepad++](https://notepad-plus-plus.org/), [NppToR](http://sourceforge.net/projects/npptor/)
  - [Vim-R-Tmux](http://manuals.bioinformatics.ucr.edu/home/programming-in-r/vim-r): R working environment based on vim and tmux
  - [Emacs](http://www.xemacs.org/Download/index.html) ([ESS add-on package](http://ess.r-project.org/))

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

  - CRAN (\>14,000 packages) general data analysis - [URL](http://cran.at.r-project.org/)
  - Bioconductor (\>2,000 packages) bioscience data analysis - [URL](http://www.bioconductor.org/)
  - Omegahat (\>90 packages) programming interfaces - [URL](https://github.com/omegahat?tab=repositories)
  - RStudio packages - [URL](https://www.rstudio.com/products/rpackages/)

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

  - **Starting R**:
    The R GUI versions, including RStudio, under Windows and Mac OS X can be
    opened by double-clicking their icons. Alternatively, one can start it by
    typing `R` in a terminal (default under Linux).

  - **Startup/Closing Behavior**:
    The R environment is controlled by hidden files in the startup directory:
    `.RData`, `.Rhistory` and `.Rprofile` (optional).

  - **Closing R**:

<!-- end list -->

``` r
q()  
```

    Save workspace image? [y/n/c]:

  - **Note**:
    When responding with `y`, then the entire R workspace will be written to
    the `.RData` file which can become very large. Often it is better to select `n` here,
    because a much better working pratice is to save an analysis protocol to an `R` or `Rmd` source file.
    This way one can quickly regenerate all data sets and objects needed in a future session.

## Navigating directories

Create an object with the assignment operator `<-` or `=`

``` r
object <- ...
```

Instead of the assignment operator one can use the `assign` function

``` r
assign("x", function(arguments))
```

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

General R command syntax

``` r
object <- function_name(arguments) 
object <- object[arguments] 
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
class is provided in the [`dplyr`](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/#dplyr-environment) section below.

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

``` r
iris$Species[1:8]
```

    ## [1] setosa setosa setosa setosa setosa setosa setosa setosa
    ## Levels: setosa versicolor virginica

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

Import of Excel file. Note: working with tab- or comma-delimited files is more flexible and preferred.

``` r
library(gdata)
myDF <- read.xls"myData.xls")
```

Import of Google Sheets. The following example imports a sample Google Sheet from [here](https://docs.google.com/spreadsheets/d/1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM/edit#gid=472150521).
Detailed instructions for interacting from R with Google Sheets with the required `googlesheets4` package are [here](https://googlesheets4.tidyverse.org/).

``` r
library(googlesheets4)
mysheet <- read_sheet("1U-32UcwZP1k3saKeaH1mbvEAOfZRdNHNkWK2GI1rpPM", skip=4)
myDF <- as.data.frame(mysheet)
myDF
```

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
write.table(iris, "clipboard", sep="\t", col.names=NA, quote=F) 
```

On Mac OS X systems

``` r
zz <- pipe('pbcopy', 'w')
write.table(iris, zz, sep="\t", col.names=NA, quote=F)
close(zz) 
```

### Homework 3A

Homework 3A: [Object Subsetting Routines and Import/Export](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw3/hw03/)

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

## dplyr Environment

Modern object classes and methods for handling `data.frame` like structures
are provided by the `dplyr` and `data.table` packages. The following gives a
short introduction to the usage and functionalities of the `dplyr` package.
More detailed tutorials on this topic can be found here:

  - [dplyr: A Grammar of Data Manipulation](https://rdrr.io/cran/dplyr/)
  - [Introduction to `dplyr`](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)
  - [Tutorial on `dplyr`](http://genomicsclass.github.io/book/pages/dplyr_tutorial.html)
  - [Cheatsheet for Joins from Jenny Bryan](http://stat545.com/bit001_dplyr-cheatsheet.html)
  - [Tibbles](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html)
  - [Intro to `data.table` package](https://www.r-bloggers.com/intro-to-the-data-table-package/)
  - [Big data with `dplyr` and `data.table`](https://www.r-bloggers.com/working-with-large-datasets-with-dplyr-and-data-table/)
  - [Fast lookups with `dplyr` and `data.table`](https://www.r-bloggers.com/fast-data-lookups-in-r-dplyr-vs-data-table/)

### Installation

The `dplyr` environment has evolved into an ecosystem of packages. To simplify
package management, one can install and load the entire collection via the
`tidyverse` package. For more details on `tidyverse` see
[here](http://tidyverse.org/).

``` r
install.packages("tidyverse")
```

### Construct a `tibble` (`tibble`)

``` r
library(tidyverse)
as_tibble(iris) # coerce data.frame to tibble tbl
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

Alternative functions producing the same result include `as_data_frame` and `tbl_df`:

``` r
as_data_frame(iris) 
tbl_df(iris) 
```

### Reading and writing tabular files

While the base R read/write utilities can be used for `data.frames`, best time
performance with the least amount of typing is achieved with the export/import
functions from the `readr` package. For very large files the `fread` function from
the `data.table` package achieves the best time performance.

#### Import with `readr`

Import functions provided by `readr` include:

  - `read_csv()`: comma separated (CSV) files
  - `read_tsv()`: tab separated files
  - `read_delim()`: general delimited files
  - `read_fwf()`: fixed width files
  - `read_table()`: tabular files where colums are separated by white-space.
  - `read_log()`: web log files

Create a sample tab delimited file for import

``` r
write_tsv(iris, "iris.txt") # Creates sample file
```

Import with `read_tsv`

``` r
iris_df <- read_tsv("iris.txt") # Import with read_tbv from readr package
iris_df
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>  
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

To import Google Sheets directly into R, see [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/#reading-and-writing-external-data).

#### Fast table import with `fread`

The `fread` function from the `data.table` package provides the best time performance for reading large
tabular files into R.

``` r
library(data.table)
iris_df <- as_data_frame(fread("iris.txt")) # Import with fread and conversion to tibble
iris_df
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>  
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

Note: to ignore lines starting with comment signs, one can pass on to `fread` a shell
command for preprocessing the file. The following example illustrates this option.

``` r
fread("grep -v '^#' iris.txt") 
```

#### Export with `readr`

Export function provided by `readr` inlcude

  - `write_delim()`: general delimited files
  - `write_csv()`: comma separated (CSV) files
  - `write_excel_csv()`: excel style CSV files
  - `write_tsv()`: tab separated files

For instance, the `write_tsv` function writes a `data.frame` or `tibble` to a tab delimited file with much nicer
default settings than the base R `write.table` function.

``` r
write_tsv(iris_df, "iris.txt")
```

### Column and row binds

The equivalents to base R’s `rbind` and `cbind` are `bind_rows` and `bind_cols`, respectively.

``` r
bind_cols(iris_df, iris_df)
```

    ## # A tibble: 150 x 10
    ##    Sepal.Length...1 Sepal.Width...2 Petal.Length...3 Petal.Width...4 Species...5 Sepal.Length...6
    ##               <dbl>           <dbl>            <dbl>           <dbl> <chr>                  <dbl>
    ##  1              5.1             3.5              1.4             0.2 setosa                   5.1
    ##  2              4.9             3                1.4             0.2 setosa                   4.9
    ##  3              4.7             3.2              1.3             0.2 setosa                   4.7
    ##  4              4.6             3.1              1.5             0.2 setosa                   4.6
    ##  5              5               3.6              1.4             0.2 setosa                   5  
    ##  6              5.4             3.9              1.7             0.4 setosa                   5.4
    ##  7              4.6             3.4              1.4             0.3 setosa                   4.6
    ##  8              5               3.4              1.5             0.2 setosa                   5  
    ##  9              4.4             2.9              1.4             0.2 setosa                   4.4
    ## 10              4.9             3.1              1.5             0.1 setosa                   4.9
    ## # … with 140 more rows, and 4 more variables: Sepal.Width...7 <dbl>, Petal.Length...8 <dbl>,
    ## #   Petal.Width...9 <dbl>, Species...10 <chr>

``` r
bind_rows(iris_df, iris_df)
```

    ## # A tibble: 300 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>  
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
    ## # … with 290 more rows

### Extract column as vector

The subsetting operators `[[` and `$`can be used to extract from a `tibble` single columns as vector.

``` r
iris_df[[5]][1:12]
```

    ##  [1] "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa"
    ## [11] "setosa" "setosa"

``` r
iris_df$Species[1:12]
```

    ##  [1] "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa"
    ## [11] "setosa" "setosa"

### Important `dplyr` functions

1.  `filter()` and `slice()`
2.  `arrange()`
3.  `select()` and `rename()`
4.  `distinct()`
5.  `mutate()` and `transmute()`
6.  `summarise()`
7.  `sample_n()` and `sample_frac()`

### Slice and filter functions

#### Filter function

``` r
filter(iris_df, Sepal.Length > 7.5, Species=="virginica")
```

    ## # A tibble: 6 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species  
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>    
    ## 1          7.6         3            6.6         2.1 virginica
    ## 2          7.7         3.8          6.7         2.2 virginica
    ## 3          7.7         2.6          6.9         2.3 virginica
    ## 4          7.7         2.8          6.7         2   virginica
    ## 5          7.9         3.8          6.4         2   virginica
    ## 6          7.7         3            6.1         2.3 virginica

#### Base R code equivalent

``` r
iris_df[iris_df[, "Sepal.Length"] > 7.5 & iris_df[, "Species"]=="virginica", ]
```

    ## # A tibble: 6 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species  
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>    
    ## 1          7.6         3            6.6         2.1 virginica
    ## 2          7.7         3.8          6.7         2.2 virginica
    ## 3          7.7         2.6          6.9         2.3 virginica
    ## 4          7.7         2.8          6.7         2   virginica
    ## 5          7.9         3.8          6.4         2   virginica
    ## 6          7.7         3            6.1         2.3 virginica

#### Including boolean operators

``` r
filter(iris_df, Sepal.Length > 7.5 | Sepal.Length < 5.5, Species=="virginica")
```

    ## # A tibble: 7 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species  
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>    
    ## 1          7.6         3            6.6         2.1 virginica
    ## 2          4.9         2.5          4.5         1.7 virginica
    ## 3          7.7         3.8          6.7         2.2 virginica
    ## 4          7.7         2.6          6.9         2.3 virginica
    ## 5          7.7         2.8          6.7         2   virginica
    ## 6          7.9         3.8          6.4         2   virginica
    ## 7          7.7         3            6.1         2.3 virginica

#### Subset rows by position

`dplyr` approach

``` r
slice(iris_df, 1:2)
```

    ## # A tibble: 2 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>  
    ## 1          5.1         3.5          1.4         0.2 setosa 
    ## 2          4.9         3            1.4         0.2 setosa

Base R code equivalent

``` r
iris_df[1:2,]
```

    ## # A tibble: 2 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>  
    ## 1          5.1         3.5          1.4         0.2 setosa 
    ## 2          4.9         3            1.4         0.2 setosa

#### Subset rows by names

Since `tibbles` do not contain row names, row wise subsetting via the `[,]` operator cannot be used.
However, the corresponding behavior can be achieved by passing to `select` a row position index
obtained by basic R intersect utilities such as `match`.

Create a suitable test `tibble`

``` r
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1
```

    ## # A tibble: 10 x 5
    ##    ids1    CA1   CA2   CA3   CA4
    ##    <chr> <int> <int> <int> <int>
    ##  1 g1        1    11    21    31
    ##  2 g2        2    12    22    32
    ##  3 g3        3    13    23    33
    ##  4 g4        4    14    24    34
    ##  5 g5        5    15    25    35
    ##  6 g6        6    16    26    36
    ##  7 g7        7    17    27    37
    ##  8 g8        8    18    28    38
    ##  9 g9        9    19    29    39
    ## 10 g10      10    20    30    40

`dplyr` approach

``` r
slice(df1, match(c("g10", "g4", "g4"), df1$ids1))
```

    ## # A tibble: 3 x 5
    ##   ids1    CA1   CA2   CA3   CA4
    ##   <chr> <int> <int> <int> <int>
    ## 1 g10      10    20    30    40
    ## 2 g4        4    14    24    34
    ## 3 g4        4    14    24    34

Base R equivalent

``` r
df1_old <- as.data.frame(df1)
rownames(df1_old) <- df1_old[,1]
df1_old[c("g10", "g4", "g4"),]
```

    ##      ids1 CA1 CA2 CA3 CA4
    ## g10   g10  10  20  30  40
    ## g4     g4   4  14  24  34
    ## g4.1   g4   4  14  24  34

### Sorting with `arrange`

Row-wise ordering based on specific columns

`dplyr` approach

``` r
arrange(iris_df, Species, Sepal.Length, Sepal.Width)
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>  
    ##  1          4.3         3            1.1         0.1 setosa 
    ##  2          4.4         2.9          1.4         0.2 setosa 
    ##  3          4.4         3            1.3         0.2 setosa 
    ##  4          4.4         3.2          1.3         0.2 setosa 
    ##  5          4.5         2.3          1.3         0.3 setosa 
    ##  6          4.6         3.1          1.5         0.2 setosa 
    ##  7          4.6         3.2          1.4         0.2 setosa 
    ##  8          4.6         3.4          1.4         0.3 setosa 
    ##  9          4.6         3.6          1           0.2 setosa 
    ## 10          4.7         3.2          1.3         0.2 setosa 
    ## # … with 140 more rows

For ordering descendingly use `desc()` function

``` r
arrange(iris_df, desc(Species), Sepal.Length, Sepal.Width)
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species  
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>    
    ##  1          4.9         2.5          4.5         1.7 virginica
    ##  2          5.6         2.8          4.9         2   virginica
    ##  3          5.7         2.5          5           2   virginica
    ##  4          5.8         2.7          5.1         1.9 virginica
    ##  5          5.8         2.7          5.1         1.9 virginica
    ##  6          5.8         2.8          5.1         2.4 virginica
    ##  7          5.9         3            5.1         1.8 virginica
    ##  8          6           2.2          5           1.5 virginica
    ##  9          6           3            4.8         1.8 virginica
    ## 10          6.1         2.6          5.6         1.4 virginica
    ## # … with 140 more rows

Base R code equivalent

``` r
iris_df[order(iris_df$Species, iris_df$Sepal.Length, iris_df$Sepal.Width), ]
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>  
    ##  1          4.3         3            1.1         0.1 setosa 
    ##  2          4.4         2.9          1.4         0.2 setosa 
    ##  3          4.4         3            1.3         0.2 setosa 
    ##  4          4.4         3.2          1.3         0.2 setosa 
    ##  5          4.5         2.3          1.3         0.3 setosa 
    ##  6          4.6         3.1          1.5         0.2 setosa 
    ##  7          4.6         3.2          1.4         0.2 setosa 
    ##  8          4.6         3.4          1.4         0.3 setosa 
    ##  9          4.6         3.6          1           0.2 setosa 
    ## 10          4.7         3.2          1.3         0.2 setosa 
    ## # … with 140 more rows

``` r
iris_df[order(iris_df$Species, decreasing=TRUE), ] 
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species  
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>    
    ##  1          6.3         3.3          6           2.5 virginica
    ##  2          5.8         2.7          5.1         1.9 virginica
    ##  3          7.1         3            5.9         2.1 virginica
    ##  4          6.3         2.9          5.6         1.8 virginica
    ##  5          6.5         3            5.8         2.2 virginica
    ##  6          7.6         3            6.6         2.1 virginica
    ##  7          4.9         2.5          4.5         1.7 virginica
    ##  8          7.3         2.9          6.3         1.8 virginica
    ##  9          6.7         2.5          5.8         1.8 virginica
    ## 10          7.2         3.6          6.1         2.5 virginica
    ## # … with 140 more rows

### Select columns with `select`

Select specific columns

``` r
select(iris_df, Species, Petal.Length, Sepal.Length)
```

    ## # A tibble: 150 x 3
    ##    Species Petal.Length Sepal.Length
    ##    <chr>          <dbl>        <dbl>
    ##  1 setosa           1.4          5.1
    ##  2 setosa           1.4          4.9
    ##  3 setosa           1.3          4.7
    ##  4 setosa           1.5          4.6
    ##  5 setosa           1.4          5  
    ##  6 setosa           1.7          5.4
    ##  7 setosa           1.4          4.6
    ##  8 setosa           1.5          5  
    ##  9 setosa           1.4          4.4
    ## 10 setosa           1.5          4.9
    ## # … with 140 more rows

Select range of columns by name

``` r
select(iris_df, Sepal.Length : Petal.Width)
```

    ## # A tibble: 150 x 4
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width
    ##           <dbl>       <dbl>        <dbl>       <dbl>
    ##  1          5.1         3.5          1.4         0.2
    ##  2          4.9         3            1.4         0.2
    ##  3          4.7         3.2          1.3         0.2
    ##  4          4.6         3.1          1.5         0.2
    ##  5          5           3.6          1.4         0.2
    ##  6          5.4         3.9          1.7         0.4
    ##  7          4.6         3.4          1.4         0.3
    ##  8          5           3.4          1.5         0.2
    ##  9          4.4         2.9          1.4         0.2
    ## 10          4.9         3.1          1.5         0.1
    ## # … with 140 more rows

Drop specific columns (here range)

``` r
select(iris_df, -(Sepal.Length : Petal.Width))
```

    ## # A tibble: 150 x 1
    ##    Species
    ##    <chr>  
    ##  1 setosa 
    ##  2 setosa 
    ##  3 setosa 
    ##  4 setosa 
    ##  5 setosa 
    ##  6 setosa 
    ##  7 setosa 
    ##  8 setosa 
    ##  9 setosa 
    ## 10 setosa 
    ## # … with 140 more rows

### Renaming columns with `rename`

`dplyr` approach

``` r
rename(iris_df, new_col_name = Species)
```

    ## # A tibble: 150 x 5
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width new_col_name
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>       
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

Base R code approach

``` r
colnames(iris_df)[colnames(iris_df)=="Species"] <- "new_col_names"
```

### Obtain unique rows with `distinct`

`dplyr` approach

``` r
distinct(iris_df, Species, .keep_all=TRUE)
```

    ## # A tibble: 3 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>     
    ## 1          5.1         3.5          1.4         0.2 setosa    
    ## 2          7           3.2          4.7         1.4 versicolor
    ## 3          6.3         3.3          6           2.5 virginica

Base R code approach

``` r
iris_df[!duplicated(iris_df$Species),]
```

    ## # A tibble: 3 x 5
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species   
    ##          <dbl>       <dbl>        <dbl>       <dbl> <chr>     
    ## 1          5.1         3.5          1.4         0.2 setosa    
    ## 2          7           3.2          4.7         1.4 versicolor
    ## 3          6.3         3.3          6           2.5 virginica

### Add columns

#### `mutate`

The `mutate` function allows to append columns to existing ones.

``` r
mutate(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)
```

    ## # A tibble: 150 x 7
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species Ratio   Sum
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>   <dbl> <dbl>
    ##  1          5.1         3.5          1.4         0.2 setosa   1.46   8.6
    ##  2          4.9         3            1.4         0.2 setosa   1.63   7.9
    ##  3          4.7         3.2          1.3         0.2 setosa   1.47   7.9
    ##  4          4.6         3.1          1.5         0.2 setosa   1.48   7.7
    ##  5          5           3.6          1.4         0.2 setosa   1.39   8.6
    ##  6          5.4         3.9          1.7         0.4 setosa   1.38   9.3
    ##  7          4.6         3.4          1.4         0.3 setosa   1.35   8  
    ##  8          5           3.4          1.5         0.2 setosa   1.47   8.4
    ##  9          4.4         2.9          1.4         0.2 setosa   1.52   7.3
    ## 10          4.9         3.1          1.5         0.1 setosa   1.58   8  
    ## # … with 140 more rows

#### `transmute`

The `transmute` function does the same as `mutate` but drops existing columns

``` r
transmute(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)
```

    ## # A tibble: 150 x 2
    ##    Ratio   Sum
    ##    <dbl> <dbl>
    ##  1  1.46   8.6
    ##  2  1.63   7.9
    ##  3  1.47   7.9
    ##  4  1.48   7.7
    ##  5  1.39   8.6
    ##  6  1.38   9.3
    ##  7  1.35   8  
    ##  8  1.47   8.4
    ##  9  1.52   7.3
    ## 10  1.58   8  
    ## # … with 140 more rows

#### `bind_cols`

The `bind_cols` function is the equivalent of `cbind` in base R. To add rows, use the corresponding
`bind_rows` function.

``` r
bind_cols(iris_df, iris_df)
```

    ## # A tibble: 150 x 10
    ##    Sepal.Length...1 Sepal.Width...2 Petal.Length...3 Petal.Width...4 Species...5 Sepal.Length...6
    ##               <dbl>           <dbl>            <dbl>           <dbl> <chr>                  <dbl>
    ##  1              5.1             3.5              1.4             0.2 setosa                   5.1
    ##  2              4.9             3                1.4             0.2 setosa                   4.9
    ##  3              4.7             3.2              1.3             0.2 setosa                   4.7
    ##  4              4.6             3.1              1.5             0.2 setosa                   4.6
    ##  5              5               3.6              1.4             0.2 setosa                   5  
    ##  6              5.4             3.9              1.7             0.4 setosa                   5.4
    ##  7              4.6             3.4              1.4             0.3 setosa                   4.6
    ##  8              5               3.4              1.5             0.2 setosa                   5  
    ##  9              4.4             2.9              1.4             0.2 setosa                   4.4
    ## 10              4.9             3.1              1.5             0.1 setosa                   4.9
    ## # … with 140 more rows, and 4 more variables: Sepal.Width...7 <dbl>, Petal.Length...8 <dbl>,
    ## #   Petal.Width...9 <dbl>, Species...10 <chr>

### Summarize data

Summary calculation on single column

``` r
summarize(iris_df, mean(Petal.Length))
```

    ## # A tibble: 1 x 1
    ##   `mean(Petal.Length)`
    ##                  <dbl>
    ## 1                 3.76

Summary calculation on many columns

``` r
summarize_all(iris_df[,1:4], mean)
```

    ## # A tibble: 1 x 4
    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width
    ##          <dbl>       <dbl>        <dbl>       <dbl>
    ## 1         5.84        3.06         3.76        1.20

Summarize by grouping column

``` r
summarize(group_by(iris_df, Species), mean(Petal.Length))
```

    ## # A tibble: 3 x 2
    ##   Species    `mean(Petal.Length)`
    ##   <chr>                     <dbl>
    ## 1 setosa                     1.46
    ## 2 versicolor                 4.26
    ## 3 virginica                  5.55

Aggregate summaries

``` r
summarize_all(group_by(iris_df, Species), mean) 
```

    ## # A tibble: 3 x 5
    ##   Species    Sepal.Length Sepal.Width Petal.Length Petal.Width
    ##   <chr>             <dbl>       <dbl>        <dbl>       <dbl>
    ## 1 setosa             5.01        3.43         1.46       0.246
    ## 2 versicolor         5.94        2.77         4.26       1.33 
    ## 3 virginica          6.59        2.97         5.55       2.03

Note: `group_by` does the looping for the user similar to `aggregate` or `tapply`.

### Merging tibbles

The `dplyr` package provides several join functions for merging `tibbles` by a common key column
similar to the `merge` function in base R. These `*_join` functions include:

  - `inner_join()`: returns join only for rows matching among both `tibbles`
  - `full_join()`: returns join for all (matching and non-matching) rows of two `tibbles`
  - `left_join()`: returns join for all rows in first `tibble`
  - `right_join()`: returns join for all rows in second `tibble`
  - `anti_join()`: returns for first `tibble` only those rows that have no match in the second one

Sample `tibbles` to illustrate `*.join` functions.

``` r
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1
```

    ## # A tibble: 10 x 5
    ##    ids1    CA1   CA2   CA3   CA4
    ##    <chr> <int> <int> <int> <int>
    ##  1 g1        1    11    21    31
    ##  2 g2        2    12    22    32
    ##  3 g3        3    13    23    33
    ##  4 g4        4    14    24    34
    ##  5 g5        5    15    25    35
    ##  6 g6        6    16    26    36
    ##  7 g7        7    17    27    37
    ##  8 g8        8    18    28    38
    ##  9 g9        9    19    29    39
    ## 10 g10      10    20    30    40

``` r
df2 <- bind_cols(data_frame(ids2=paste0("g", c(2,5,11,12))), as_data_frame(matrix(1:16, 4, 4, dimnames=list(1:4, paste0("CB", 1:4)))))
df2
```

    ## # A tibble: 4 x 5
    ##   ids2    CB1   CB2   CB3   CB4
    ##   <chr> <int> <int> <int> <int>
    ## 1 g2        1     5     9    13
    ## 2 g5        2     6    10    14
    ## 3 g11       3     7    11    15
    ## 4 g12       4     8    12    16

#### Inner join

``` r
inner_join(df1, df2, by=c("ids1"="ids2"))
```

    ## # A tibble: 2 x 9
    ##   ids1    CA1   CA2   CA3   CA4   CB1   CB2   CB3   CB4
    ##   <chr> <int> <int> <int> <int> <int> <int> <int> <int>
    ## 1 g2        2    12    22    32     1     5     9    13
    ## 2 g5        5    15    25    35     2     6    10    14

#### Left join

``` r
left_join(df1, df2, by=c("ids1"="ids2"))
```

    ## # A tibble: 10 x 9
    ##    ids1    CA1   CA2   CA3   CA4   CB1   CB2   CB3   CB4
    ##    <chr> <int> <int> <int> <int> <int> <int> <int> <int>
    ##  1 g1        1    11    21    31    NA    NA    NA    NA
    ##  2 g2        2    12    22    32     1     5     9    13
    ##  3 g3        3    13    23    33    NA    NA    NA    NA
    ##  4 g4        4    14    24    34    NA    NA    NA    NA
    ##  5 g5        5    15    25    35     2     6    10    14
    ##  6 g6        6    16    26    36    NA    NA    NA    NA
    ##  7 g7        7    17    27    37    NA    NA    NA    NA
    ##  8 g8        8    18    28    38    NA    NA    NA    NA
    ##  9 g9        9    19    29    39    NA    NA    NA    NA
    ## 10 g10      10    20    30    40    NA    NA    NA    NA

#### Right join

``` r
right_join(df1, df2, by=c("ids1"="ids2"))
```

    ## # A tibble: 4 x 9
    ##   ids1    CA1   CA2   CA3   CA4   CB1   CB2   CB3   CB4
    ##   <chr> <int> <int> <int> <int> <int> <int> <int> <int>
    ## 1 g2        2    12    22    32     1     5     9    13
    ## 2 g5        5    15    25    35     2     6    10    14
    ## 3 g11      NA    NA    NA    NA     3     7    11    15
    ## 4 g12      NA    NA    NA    NA     4     8    12    16

#### Full join

``` r
full_join(df1, df2, by=c("ids1"="ids2"))
```

    ## # A tibble: 12 x 9
    ##    ids1    CA1   CA2   CA3   CA4   CB1   CB2   CB3   CB4
    ##    <chr> <int> <int> <int> <int> <int> <int> <int> <int>
    ##  1 g1        1    11    21    31    NA    NA    NA    NA
    ##  2 g2        2    12    22    32     1     5     9    13
    ##  3 g3        3    13    23    33    NA    NA    NA    NA
    ##  4 g4        4    14    24    34    NA    NA    NA    NA
    ##  5 g5        5    15    25    35     2     6    10    14
    ##  6 g6        6    16    26    36    NA    NA    NA    NA
    ##  7 g7        7    17    27    37    NA    NA    NA    NA
    ##  8 g8        8    18    28    38    NA    NA    NA    NA
    ##  9 g9        9    19    29    39    NA    NA    NA    NA
    ## 10 g10      10    20    30    40    NA    NA    NA    NA
    ## 11 g11      NA    NA    NA    NA     3     7    11    15
    ## 12 g12      NA    NA    NA    NA     4     8    12    16

#### Anti join

``` r
anti_join(df1, df2, by=c("ids1"="ids2"))
```

    ## # A tibble: 8 x 5
    ##   ids1    CA1   CA2   CA3   CA4
    ##   <chr> <int> <int> <int> <int>
    ## 1 g1        1    11    21    31
    ## 2 g3        3    13    23    33
    ## 3 g4        4    14    24    34
    ## 4 g6        6    16    26    36
    ## 5 g7        7    17    27    37
    ## 6 g8        8    18    28    38
    ## 7 g9        9    19    29    39
    ## 8 g10      10    20    30    40

For additional join options users want to cosult the `*_join` help pages.

### Chaining

To simplify chaining of serveral operations, `dplyr` provides the `%>%`
operator, where `x %>% f(y)` turns into `f(x, y)`. This way one can pipe
together multiple operations by writing them from left-to-right or
top-to-bottom. This makes for easy to type and readable code.

#### Example 1

Series of data manipulations and export

``` r
read_tsv("iris.txt") %>% # Import with read_tbv from readr package
    as_tibble() %>% # Declare tibble to use
    select(Sepal.Length:Species) %>% # Select columns
    filter(Species=="setosa") %>% # Filter rows by some value
    arrange(Sepal.Length) %>% # Sort by some column
    mutate(Subtract=Petal.Length - Petal.Width) # Calculate and append
```

    ## # A tibble: 50 x 6
    ##    Sepal.Length Sepal.Width Petal.Length Petal.Width Species Subtract
    ##           <dbl>       <dbl>        <dbl>       <dbl> <chr>      <dbl>
    ##  1          4.3         3            1.1         0.1 setosa      1   
    ##  2          4.4         2.9          1.4         0.2 setosa      1.2 
    ##  3          4.4         3            1.3         0.2 setosa      1.1 
    ##  4          4.4         3.2          1.3         0.2 setosa      1.1 
    ##  5          4.5         2.3          1.3         0.3 setosa      1   
    ##  6          4.6         3.1          1.5         0.2 setosa      1.3 
    ##  7          4.6         3.4          1.4         0.3 setosa      1.10
    ##  8          4.6         3.6          1           0.2 setosa      0.8 
    ##  9          4.6         3.2          1.4         0.2 setosa      1.2 
    ## 10          4.7         3.2          1.3         0.2 setosa      1.1 
    ## # … with 40 more rows

``` r
    # write_tsv("iris.txt") # Export to file, omitted here to show result 
```

#### Example 2

Series of summary calculations for grouped data (`group_by`)

``` r
iris_df %>% # Declare tibble to use 
    group_by(Species) %>% # Group by species
    summarize(Mean_Sepal.Length=mean(Sepal.Length), 
              Max_Sepal.Length=max(Sepal.Length),
              Min_Sepal.Length=min(Sepal.Length),
              SD_Sepal.Length=sd(Sepal.Length),
              Total=n()) 
```

    ## # A tibble: 3 x 6
    ##   Species    Mean_Sepal.Length Max_Sepal.Length Min_Sepal.Length SD_Sepal.Length Total
    ##   <chr>                  <dbl>            <dbl>            <dbl>           <dbl> <int>
    ## 1 setosa                  5.01              5.8              4.3           0.352    50
    ## 2 versicolor              5.94              7                4.9           0.516    50
    ## 3 virginica               6.59              7.9              4.9           0.636    50

#### Example 3

Combining `dplyr` chaining with `ggplot`

``` r
iris_df %>% 
    group_by(Species) %>% 
    summarize_all(mean) %>% 
    reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%
    ggplot(aes(Samples, Values, fill = Species)) + 
        geom_bar(position="dodge", stat="identity")
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plyr_chaining3-1.png" width="672" />

## SQLite Databases

`SQLite` is a lightweight relational database solution. The `RSQLite` package provides an easy to use interface to create, manage and query `SQLite` databases directly from R. Basic instructions
for using `SQLite` from the command-line are available [here](https://www.sqlite.org/cli.html). A short introduction to `RSQLite` is available [here](https://github.com/rstats-db/RSQLite/blob/master/vignettes/RSQLite.Rmd).

## Loading data into SQLite databases

The following loads two `data.frames` derived from the `iris` data set (here `mydf1` and `mydf2`)
into an SQLite database (here `test.db`).

``` r
library(RSQLite)
unlink("test.db") # Delete any existing test.db
mydb <- dbConnect(SQLite(), "test.db") # Creates database file test.db
mydf1 <- data.frame(ids=paste0("id", seq_along(iris[,1])), iris)
mydf2 <- mydf1[sample(seq_along(mydf1[,1]), 10),]
dbWriteTable(mydb, "mydf1", mydf1)
dbWriteTable(mydb, "mydf2", mydf2)
```

### List names of tables in database

``` r
dbListTables(mydb)
```

    ## [1] "mydf1" "mydf2"

### Import table into `data.frame`

``` r
dbGetQuery(mydb, 'SELECT * FROM mydf2')
```

    ##      ids Sepal.Length Sepal.Width Petal.Length Petal.Width    Species
    ## 1  id122          5.6         2.8          4.9         2.0  virginica
    ## 2   id68          5.8         2.7          4.1         1.0 versicolor
    ## 3   id26          5.0         3.0          1.6         0.2     setosa
    ## 4  id150          5.9         3.0          5.1         1.8  virginica
    ## 5   id77          6.8         2.8          4.8         1.4 versicolor
    ## 6   id59          6.6         2.9          4.6         1.3 versicolor
    ## 7   id45          5.1         3.8          1.9         0.4     setosa
    ## 8  id143          5.8         2.7          5.1         1.9  virginica
    ## 9   id20          5.1         3.8          1.5         0.3     setosa
    ## 10  id12          4.8         3.4          1.6         0.2     setosa

### Query database

``` r
dbGetQuery(mydb, 'SELECT * FROM mydf1 WHERE "Sepal.Length" < 4.6')
```

    ##    ids Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ## 1  id9          4.4         2.9          1.4         0.2  setosa
    ## 2 id14          4.3         3.0          1.1         0.1  setosa
    ## 3 id39          4.4         3.0          1.3         0.2  setosa
    ## 4 id42          4.5         2.3          1.3         0.3  setosa
    ## 5 id43          4.4         3.2          1.3         0.2  setosa

### Join tables

The two tables can be joined on the shared `ids` column as follows.

``` r
dbGetQuery(mydb, 'SELECT * FROM mydf1, mydf2 WHERE mydf1.ids = mydf2.ids')
```

    ##      ids Sepal.Length Sepal.Width Petal.Length Petal.Width    Species   ids Sepal.Length
    ## 1   id12          4.8         3.4          1.6         0.2     setosa  id12          4.8
    ## 2   id20          5.1         3.8          1.5         0.3     setosa  id20          5.1
    ## 3   id26          5.0         3.0          1.6         0.2     setosa  id26          5.0
    ## 4   id45          5.1         3.8          1.9         0.4     setosa  id45          5.1
    ## 5   id59          6.6         2.9          4.6         1.3 versicolor  id59          6.6
    ## 6   id68          5.8         2.7          4.1         1.0 versicolor  id68          5.8
    ## 7   id77          6.8         2.8          4.8         1.4 versicolor  id77          6.8
    ## 8  id122          5.6         2.8          4.9         2.0  virginica id122          5.6
    ## 9  id143          5.8         2.7          5.1         1.9  virginica id143          5.8
    ## 10 id150          5.9         3.0          5.1         1.8  virginica id150          5.9
    ##    Sepal.Width Petal.Length Petal.Width    Species
    ## 1          3.4          1.6         0.2     setosa
    ## 2          3.8          1.5         0.3     setosa
    ## 3          3.0          1.6         0.2     setosa
    ## 4          3.8          1.9         0.4     setosa
    ## 5          2.9          4.6         1.3 versicolor
    ## 6          2.7          4.1         1.0 versicolor
    ## 7          2.8          4.8         1.4 versicolor
    ## 8          2.8          4.9         2.0  virginica
    ## 9          2.7          5.1         1.9  virginica
    ## 10         3.0          5.1         1.8  virginica

## Graphics in R

### Advantages

  - Powerful environment for visualizing scientific data
  - Integrated graphics and statistics infrastructure
  - Publication quality graphics
  - Fully programmable
  - Highly reproducible
  - Full [LaTeX](http://www.latex-project.org/) and Markdown support via `knitr` and `R markdown`
  - Vast number of R packages with graphics utilities

### Documentation for R Graphics

**General**

  - Graphics Task Page - [URL](http://cran.r-project.org/web/views/Graphics.html)
  - R Graph Gallery - [URL](http://addictedtor.free.fr/graphiques/allgraph.php)
  - R Graphical Manual - [URL](http://cged.genes.nig.ac.jp/RGM2/index.php)
  - Paul Murrell’s book R (Grid) Graphics - [URL](http://www.stat.auckland.ac.nz/~paul/RGraphics/rgraphics.html)

**Interactive graphics**

  - rggobi\` (GGobi) - [URL](http://www.ggobi.org/)
  - `iplots` - [URL](http://www.rosuda.org/iplots/)
  - Open GL (`rgl`) - [URL](http://rgl.neoscientists.org/gallery.shtml)

### Graphics Environments

**Viewing and saving graphics in R**

  - On-screen graphics
  - postscript, pdf, svg
  - jpeg, png, wmf, tiff, …

**Four major graphic environments**

1)  Low-level infrastructure

<!-- end list -->

  - R Base Graphics (low- and high-level)
  - `grid`: [Manual](http://www.stat.auckland.ac.nz/~paul/grid/grid.html)

<!-- end list -->

2)  High-level infrastructure
    \\begin{itemize}

<!-- end list -->

  - `lattice`: [Manual](http://lmdvr.r-forge.r-project.org), [Intro](http://www.his.sunderland.ac.uk/~cs0her/Statistics/UsingLatticeGraphicsInR.htm), [Book](http://www.amazon.com/Lattice-Multivariate-Data-Visualization-Use/dp/0387759689)
  - `ggplot2`: [Manual](http://had.co.nz/ggplot2/), [Intro](http://www.ling.upenn.edu/~joseff/rstudy/summer2010_ggplot2_intro.html), [Book](http://had.co.nz/ggplot2/book/)

### Base Graphics: Overview

**Important high-level plotting functions**

  - `plot`: generic x-y plotting
  - `barplot`: bar plots
  - `boxplot`: box-and-whisker plot
  - `hist`: histograms
  - `pie`: pie charts
  - `dotchart`: cleveland dot plots
  - `image, heatmap, contour, persp`: functions to generate image-like plots
  - `qqnorm, qqline, qqplot`: distribution comparison plots
  - `pairs, coplot`: display of multivariant data

**Help on graphics functions**

  - `?myfct`
  - `?plot`
  - `?par`

#### Preferred Object Types

  - Matrices and data frames
  - Vectors
  - Named vectors

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

\_\_Important arguments\_

  - `mar`: specifies the margin sizes around the plotting area in order: `c(bottom, left, top, right)`
  - `col`: color of symbols
  - `pch`: type of symbols, samples: `example(points)`
  - `lwd`: size of symbols
  - `cex.*`: control font sizes
  - For details see `?par`

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

Homework 3B: [Scatter Plots](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw3/hw03/)

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

Homework 3C: [Bar Plots](https://girke.bioinformatics.ucr.edu/GEN242/assignments/homework/hw3/hw03/)

## Analysis Routine

### Overview

The following exercise introduces a variety of useful data analysis utilities in R.

### Analysis Routine: Data Import

  - **Step 1**: To get started with this exercise, direct your R session to a dedicated workshop directory and download into this directory the following sample tables. Then import the files into Excel and save them as tab delimited text files.
    
      - [MolecularWeight\_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls)
      - [TargetP\_analysis\_tair7.xls](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls)

**Import the tables into R**

Import molecular weight table

``` r
my_mw <- read.delim(file="MolecularWeight_tair7.xls", header=T, sep="\t") 
my_mw[1:2,]
```

Import subcelluar targeting table

``` r
my_target <- read.delim(file="TargetP_analysis_tair7.xls", header=T, sep="\t") 
my_target[1:2,]
```

Online import of molecular weight table

``` r
my_mw <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls", header=T, sep="\t") 
my_mw[1:2,]
```

    ##   Sequence.id Molecular.Weight.Da. Residues
    ## 1 AT1G08520.1                83285      760
    ## 2 AT1G08530.1                27015      257

Online import of subcelluar targeting table

``` r
my_target <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls", header=T, sep="\t") 
my_target[1:2,]
```

    ##      GeneName Loc   cTP   mTP    SP other
    ## 1 AT1G08520.1   C 0.822 0.137 0.029 0.039
    ## 2 AT1G08530.1   C 0.817 0.058 0.010 0.100

### Merging Data Frames

  - **Step 2**: Assign uniform gene ID column titles

<!-- end list -->

``` r
colnames(my_target)[1] <- "ID"
colnames(my_mw)[1] <- "ID" 
```

  - **Step 3**: Merge the two tables based on common ID field

<!-- end list -->

``` r
my_mw_target <- merge(my_mw, my_target, by.x="ID", by.y="ID", all.x=T)
```

  - **Step 4**: Shorten one table before the merge and then remove the non-matching rows (NAs) in the merged file

<!-- end list -->

``` r
my_mw_target2a <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all.x=T)  # To remove non-matching rows, use the argument setting 'all=F'.
my_mw_target2 <- na.omit(my_mw_target2a) # Removes rows containing "NAs" (non-matching rows).
```

  - **Homework 3D**: How can the merge function in the previous step be executed so that only the common rows among the two data frames are returned? Prove that both methods - the two step version with `na.omit` and your method - return identical results.
  - **Homework 3E**: Replace all `NAs` in the data frame `my_mw_target2a` with zeros.

### Filtering Data

  - **Step 5**: Retrieve all records with a value of greater than 100,000 in ‘MW’ column and ‘C’ value in ‘Loc’ column (targeted to chloroplast).

<!-- end list -->

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

  - **Homework 3F**: How many protein entries in the `my`\_mw`_target` data frame have a MW of greater then 4,000 and less then 5,000. Subset the data frame accordingly and sort it by MW to check that your result is correct.

### String Substitutions

  - **Step 6**: Use a regular expression in a substitute function to generate a separate ID column that lacks the gene model extensions.
    \<\<label=Exercise 4.7, eval=TRUE, echo=TRUE, keep.source=TRUE\>\>=

<!-- end list -->

``` r
my_mw_target3 <- data.frame(loci=gsub("\\..*", "", as.character(my_mw_target[,1]), perl = TRUE), my_mw_target)
my_mw_target3[1:3,1:8]
```

    ##        loci          ID Molecular.Weight.Da. Residues Loc  cTP   mTP    SP
    ## 1 AT1G01010 AT1G01010.1                49426      429   _ 0.10 0.090 0.075
    ## 2 AT1G01020 AT1G01020.1                28092      245   * 0.01 0.636 0.158
    ## 3 AT1G01020 AT1G01020.2                21711      191   * 0.01 0.636 0.158

  - **Homework 3G**: Retrieve those rows in `my_mw_target3` where the second column contains the following identifiers: `c("AT5G52930.1", "AT4G18950.1", "AT1G15385.1", "AT4G36500.1", "AT1G67530.1")`. Use the `%in%` function for this query. As an alternative approach, assign the second column to the row index of the data frame and then perform the same query again using the row index. Explain the difference of the two methods.

### Calculations on Data Frames

  - **Step 7**: Count the number of duplicates in the loci column with the `table` function and append the result to the data frame with the `cbind` function.

<!-- end list -->

``` r
mycounts <- table(my_mw_target3[,1])[my_mw_target3[,1]]
my_mw_target4 <- cbind(my_mw_target3, Freq=mycounts[as.character(my_mw_target3[,1])]) 
```

  - **Step 8**: Perform a vectorized devision of columns 3 and 4 (average AA weight per protein)

<!-- end list -->

``` r
data.frame(my_mw_target4, avg_AA_WT=(my_mw_target4[,3] / my_mw_target4[,4]))[1:2,5:11] 
```

    ##   Loc  cTP   mTP    SP other Freq.Var1 Freq.Freq
    ## 1   _ 0.10 0.090 0.075 0.925 AT1G01010         1
    ## 2   * 0.01 0.636 0.158 0.448 AT1G01020         2

  - **Step 9**: Calculate for each row the mean and standard deviation across several columns

<!-- end list -->

``` r
mymean <- apply(my_mw_target4[,6:9], 1, mean)
mystdev <- apply(my_mw_target4[,6:9], 1, sd, na.rm=TRUE)
data.frame(my_mw_target4, mean=mymean, stdev=mystdev)[1:2,5:12] 
```

    ##   Loc  cTP   mTP    SP other Freq.Var1 Freq.Freq   mean
    ## 1   _ 0.10 0.090 0.075 0.925 AT1G01010         1 0.2975
    ## 2   * 0.01 0.636 0.158 0.448 AT1G01020         2 0.3130

### Plotting Example

  - **Step 10**: Generate scatter plot columns: ‘MW’ and ‘Residues’

<!-- end list -->

``` r
plot(my_mw_target4[1:500,3:4], col="red")
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/plot_example-1.png" width="672" />

### Export Results and Run Entire Exercise as Script

  - **Step 11**: Write the data frame `my_mw_target4` into a tab-delimited text file and inspect it in Excel.

<!-- end list -->

``` r
write.table(my_mw_target4, file="my_file.xls", quote=F, sep="\t", col.names = NA) 
```

  - **Homework 3H**: Write all commands from this exercise into an R script named `exerciseRbasics.R`, or download it from [here](http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/exerciseRbasics.R). Then execute the script with the `source` function like this: `source("exerciseRbasics.R")`. This will run all commands of this exercise and generate the corresponding output files in the current working directory.

<!-- end list -->

``` r
source("exerciseRbasics.R")
```

## R Markdown

### Overview

R Markdown combines markdown (an easy to write plain text format) with embedded
R code chunks. When compiling R Markdown documents, the code components can be
evaluated so that both the code and its output can be included in the final
document. This makes analysis reports highly reproducible by allowing to automatically
regenerate them when the underlying R code or data changes. R Markdown
documents (`.Rmd` files) can be rendered to various formats including HTML and
PDF. The R code in an `.Rmd` document is processed by `knitr`, while the
resulting `.md` file is rendered by `pandoc` to the final output formats
(*e.g.* HTML or PDF). Historically, R Markdown is an extension of the older
`Sweave/Latex` environment. Rendering of mathematical expressions and reference
management is also supported by R Markdown using embedded Latex syntax and
Bibtex, respectively.

### Quick Start

#### Install R Markdown

``` r
install.packages("rmarkdown")
```

#### Initialize a new R Markdown (`Rmd`) script

To minimize typing, it can be helful to start with an R Markdown template and
then modify it as needed. Note the file name of an R Markdown scirpt needs to
have the extension `.Rmd`. Template files for the following examples are available
here:

  - R Markdown sample script: [`sample.Rmd`](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/07_Rbasics/sample.Rmd)
  - Bibtex file for handling citations and reference section: [`bibtex.bib`](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/07_Rbasics/bibtex.bib)

Users want to download these files, open the `sample.Rmd` file with their preferred R IDE
(*e.g.* RStudio, vim or emacs), initilize an R session and then direct their R session to
the location of these two files.

#### Metadata section

The metadata section (YAML header) in an R Markdown script defines how it will be processed and
rendered. The metadata section also includes both title, author, and date information as well as
options for customizing the output format. For instance, PDF and HTML output can be defined
with `pdf_document` and `html_document`, respectively. The `BiocStyle::` prefix will use the
formatting style of the [`BiocStyle`](http://bioconductor.org/packages/release/bioc/html/BiocStyle.html)
package from Bioconductor.

``` 
 ---
title: "My First R Markdown Document"
author: "Author: First Last"
date: "Last update: 07 April, 2021"
output:
  BiocStyle::html_document:
    toc: true
    toc_depth: 3
    fig_caption: yes

fontsize: 14pt
bibliography: bibtex.bib
 ---
```

#### Render `Rmd` script

An R Markdown script can be evaluated and rendered with the following `render` command or by pressing the `knit` button in RStudio.
The `output_format` argument defines the format of the output (*e.g.* `html_document`). The setting `output_format="all"` will generate
all supported output formats. Alternatively, one can specify several output formats in the metadata section as shown in the above example.

``` r
rmarkdown::render("sample.Rmd", clean=TRUE, output_format="html_document")
```

The following shows two options how to run the rendering from the command-line.

``` sh
$ Rscript -e "rmarkdown::render('sample.Rmd', clean=TRUE)"
```

Alternatively, one can use a Makefile to evaluate and render an R Markdown
script. A sample Makefile for rendering the above `sample.Rmd` can be
downloaded [`here`](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/07_Rbasics/Makefile).
To apply it to a custom `Rmd` file, one needs open the Makefile in a text
editor and change the value assigned to `MAIN` (line 13) to the base name of
the corresponding `.Rmd` file (*e.g.* assign `systemPipeRNAseq` if the file
name is `systemPipeRNAseq.Rmd`). To execute the `Makefile`, run the following
command from the command-line.

``` sh
$ make -B
```

#### R code chunks

R Code Chunks can be embedded in an R Markdown script by using three backticks
at the beginning of a new line along with arguments enclosed in curly braces
controlling the behavior of the code. The following lines contain the
plain R code. A code chunk is terminated by a new line starting with three backticks.
The following shows an example of such a code chunk. Note the backslashes are
not part of it. They have been added to print the code chunk syntax in this document.

```` 
    ```\{r code_chunk_name, eval=FALSE\}
    x <- 1:10
    ```
````

The following lists the most important arguments to control the behavior of R code chunks:

  - `r`: specifies language for code chunk, here R
  - `chode_chunk_name`: name of code chunk; this name needs to be unique
  - `eval`: if assigned `TRUE` the code will be evaluated
  - `warning`: if assigned `FALSE` warnings will not be shown
  - `message`: if assigned `FALSE` messages will not be shown
  - `cache`: if assigned `TRUE` results will be cached to reuse in future rendering instances
  - `fig.height`: allows to specify height of figures in inches
  - `fig.width`: allows to specify width of figures in inches

For more details on code chunk options see [here](https://www.rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf).

#### Learning Markdown

The basic syntax of Markdown and derivatives like kramdown is extremely easy to learn. Rather
than providing another introduction on this topic, here are some useful sites for learning Markdown:

  - [Markdown Intro on GitHub](https://guides.github.com/features/mastering-markdown/)
  - [Markdown Cheet Sheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
  - [Markdown Basics from RStudio](http://rmarkdown.rstudio.com/authoring_basics.html)
  - [R Markdown Cheat Sheet](http://www.rstudio.com/wp-content/uploads/2015/02/rmarkdown-cheatsheet.pdf)
  - [kramdown Syntax](http://kramdown.gettalong.org/syntax.html)

#### Tables

There are several ways to render tables. First, they can be printed within the R code chunks. Second,
much nicer formatted tables can be generated with the functions `kable`, `pander` or `xtable`. The following
example uses `kable` from the `knitr` package.

``` r
library(knitr)
kable(iris[1:12,])
```

| Sepal.Length | Sepal.Width | Petal.Length | Petal.Width | Species |
| -----------: | ----------: | -----------: | ----------: | :------ |
|          5.1 |         3.5 |          1.4 |         0.2 | setosa  |
|          4.9 |         3.0 |          1.4 |         0.2 | setosa  |
|          4.7 |         3.2 |          1.3 |         0.2 | setosa  |
|          4.6 |         3.1 |          1.5 |         0.2 | setosa  |
|          5.0 |         3.6 |          1.4 |         0.2 | setosa  |
|          5.4 |         3.9 |          1.7 |         0.4 | setosa  |
|          4.6 |         3.4 |          1.4 |         0.3 | setosa  |
|          5.0 |         3.4 |          1.5 |         0.2 | setosa  |
|          4.4 |         2.9 |          1.4 |         0.2 | setosa  |
|          4.9 |         3.1 |          1.5 |         0.1 | setosa  |
|          5.4 |         3.7 |          1.5 |         0.2 | setosa  |
|          4.8 |         3.4 |          1.6 |         0.2 | setosa  |

A much more elegant and powerful solution is to create fully interactive tables
with the [`DT` package](https://rstudio.github.io/DT/). This JavaScirpt based environment provides a wrapper to the DataTables library
using jQuery. The resulting tables can be sorted, queried and resized by the
user. For an example see below.

``` r
library(DT)
datatable(iris, filter = 'top', options = list(
      pageLength = 100, scrollX = TRUE, scrollY = "600px", autoWidth = TRUE
     ))
```

<div id="htmlwidget-1" class="datatables html-widget" style="width:100%;height:auto;">

</div>

<script type="application/json" data-for="htmlwidget-1">{"x":{"filter":"top","filterHTML":"<tr>\n  <td><\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"4.3\" data-max=\"7.9\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"2\" data-max=\"4.4\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"1\" data-max=\"6.9\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none; position: absolute; width: 200px;\">\n      <div data-min=\"0.1\" data-max=\"2.5\" data-scale=\"1\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"factor\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"width: 100%; display: none;\">\n      <select multiple=\"multiple\" style=\"width: 100%;\" data-options=\"[&quot;setosa&quot;,&quot;versicolor&quot;,&quot;virginica&quot;]\"><\/select>\n    <\/div>\n  <\/td>\n<\/tr>","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150"],[5.1,4.9,4.7,4.6,5,5.4,4.6,5,4.4,4.9,5.4,4.8,4.8,4.3,5.8,5.7,5.4,5.1,5.7,5.1,5.4,5.1,4.6,5.1,4.8,5,5,5.2,5.2,4.7,4.8,5.4,5.2,5.5,4.9,5,5.5,4.9,4.4,5.1,5,4.5,4.4,5,5.1,4.8,5.1,4.6,5.3,5,7,6.4,6.9,5.5,6.5,5.7,6.3,4.9,6.6,5.2,5,5.9,6,6.1,5.6,6.7,5.6,5.8,6.2,5.6,5.9,6.1,6.3,6.1,6.4,6.6,6.8,6.7,6,5.7,5.5,5.5,5.8,6,5.4,6,6.7,6.3,5.6,5.5,5.5,6.1,5.8,5,5.6,5.7,5.7,6.2,5.1,5.7,6.3,5.8,7.1,6.3,6.5,7.6,4.9,7.3,6.7,7.2,6.5,6.4,6.8,5.7,5.8,6.4,6.5,7.7,7.7,6,6.9,5.6,7.7,6.3,6.7,7.2,6.2,6.1,6.4,7.2,7.4,7.9,6.4,6.3,6.1,7.7,6.3,6.4,6,6.9,6.7,6.9,5.8,6.8,6.7,6.7,6.3,6.5,6.2,5.9],[3.5,3,3.2,3.1,3.6,3.9,3.4,3.4,2.9,3.1,3.7,3.4,3,3,4,4.4,3.9,3.5,3.8,3.8,3.4,3.7,3.6,3.3,3.4,3,3.4,3.5,3.4,3.2,3.1,3.4,4.1,4.2,3.1,3.2,3.5,3.6,3,3.4,3.5,2.3,3.2,3.5,3.8,3,3.8,3.2,3.7,3.3,3.2,3.2,3.1,2.3,2.8,2.8,3.3,2.4,2.9,2.7,2,3,2.2,2.9,2.9,3.1,3,2.7,2.2,2.5,3.2,2.8,2.5,2.8,2.9,3,2.8,3,2.9,2.6,2.4,2.4,2.7,2.7,3,3.4,3.1,2.3,3,2.5,2.6,3,2.6,2.3,2.7,3,2.9,2.9,2.5,2.8,3.3,2.7,3,2.9,3,3,2.5,2.9,2.5,3.6,3.2,2.7,3,2.5,2.8,3.2,3,3.8,2.6,2.2,3.2,2.8,2.8,2.7,3.3,3.2,2.8,3,2.8,3,2.8,3.8,2.8,2.8,2.6,3,3.4,3.1,3,3.1,3.1,3.1,2.7,3.2,3.3,3,2.5,3,3.4,3],[1.4,1.4,1.3,1.5,1.4,1.7,1.4,1.5,1.4,1.5,1.5,1.6,1.4,1.1,1.2,1.5,1.3,1.4,1.7,1.5,1.7,1.5,1,1.7,1.9,1.6,1.6,1.5,1.4,1.6,1.6,1.5,1.5,1.4,1.5,1.2,1.3,1.4,1.3,1.5,1.3,1.3,1.3,1.6,1.9,1.4,1.6,1.4,1.5,1.4,4.7,4.5,4.9,4,4.6,4.5,4.7,3.3,4.6,3.9,3.5,4.2,4,4.7,3.6,4.4,4.5,4.1,4.5,3.9,4.8,4,4.9,4.7,4.3,4.4,4.8,5,4.5,3.5,3.8,3.7,3.9,5.1,4.5,4.5,4.7,4.4,4.1,4,4.4,4.6,4,3.3,4.2,4.2,4.2,4.3,3,4.1,6,5.1,5.9,5.6,5.8,6.6,4.5,6.3,5.8,6.1,5.1,5.3,5.5,5,5.1,5.3,5.5,6.7,6.9,5,5.7,4.9,6.7,4.9,5.7,6,4.8,4.9,5.6,5.8,6.1,6.4,5.6,5.1,5.6,6.1,5.6,5.5,4.8,5.4,5.6,5.1,5.1,5.9,5.7,5.2,5,5.2,5.4,5.1],[0.2,0.2,0.2,0.2,0.2,0.4,0.3,0.2,0.2,0.1,0.2,0.2,0.1,0.1,0.2,0.4,0.4,0.3,0.3,0.3,0.2,0.4,0.2,0.5,0.2,0.2,0.4,0.2,0.2,0.2,0.2,0.4,0.1,0.2,0.2,0.2,0.2,0.1,0.2,0.2,0.3,0.3,0.2,0.6,0.4,0.3,0.2,0.2,0.2,0.2,1.4,1.5,1.5,1.3,1.5,1.3,1.6,1,1.3,1.4,1,1.5,1,1.4,1.3,1.4,1.5,1,1.5,1.1,1.8,1.3,1.5,1.2,1.3,1.4,1.4,1.7,1.5,1,1.1,1,1.2,1.6,1.5,1.6,1.5,1.3,1.3,1.3,1.2,1.4,1.2,1,1.3,1.2,1.3,1.3,1.1,1.3,2.5,1.9,2.1,1.8,2.2,2.1,1.7,1.8,1.8,2.5,2,1.9,2.1,2,2.4,2.3,1.8,2.2,2.3,1.5,2.3,2,2,1.8,2.1,1.8,1.8,1.8,2.1,1.6,1.9,2,2.2,1.5,1.4,2.3,2.4,1.8,1.8,2.1,2.4,2.3,1.9,2.3,2.5,2.3,1.9,2,2.3,1.8],["setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","setosa","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","versicolor","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica","virginica"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Sepal.Length<\/th>\n      <th>Sepal.Width<\/th>\n      <th>Petal.Length<\/th>\n      <th>Petal.Width<\/th>\n      <th>Species<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":100,"scrollX":true,"scrollY":"600px","autoWidth":true,"columnDefs":[{"className":"dt-right","targets":[1,2,3,4]},{"orderable":false,"targets":0}],"order":[],"orderClasses":false,"orderCellsTop":true}},"evals":[],"jsHooks":[]}</script>

#### Figures

Plots generated by the R code chunks in an R Markdown document can be automatically
inserted in the output file. The size of the figure can be controlled with the `fig.height`
and `fig.width` arguments.

``` r
library(ggplot2)
dsmall <- diamonds[sample(nrow(diamonds), 1000), ]
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
```

<img src="/en/tutorials/rbasics/Rbasics_files/figure-html/some_jitter_plot-1.png" width="672" />

Sometimes it can be useful to explicitly write an image to a file and then insert that
image into the final document by referencing its file name in the R Markdown source. For
instance, this can be useful for time consuming analyses. The following code will generate a
file named `myplot.png`. To insert the file in the final document, one can use standard
Markdown or HTML syntax, *e.g.*: `<img src="myplot.png"/>`.

``` r
png("images/myplot.png")
ggplot(dsmall, aes(color, price/carat)) + geom_jitter(alpha = I(1 / 2), aes(color=color))
dev.off()
```

    ## png 
    ##   2

<center>

<img title="some_title" src="../images/myplot.png"/>

</center>

#### Inline R code

To evaluate R code inline, one can enclose an R expression with a single back-tick
followed by `r` and then the actual expression. For instance, the back-ticked version
of ‘r 1 + 1’ evaluates to 2 and ‘r pi’ evaluates to 3.1415927.

#### Mathematical equations

To render mathematical equations, one can use standard Latex syntax. When expressions are
enclosed with single `$` signs then they will be shown inline, while
enclosing them with double `$$` signs will show them in display mode. For instance, the following
Latex syntax `d(X,Y) = \sqrt[]{ \sum_{i=1}^{n}{(x_{i}-y_{i})^2} }` renders in display mode as follows:

`$$d(X,Y) = \sqrt[]{ \sum_{i=1}^{n}{(x_{i}-y_{i})^2} }$$`

#### Citations and bibliographies

Citations and bibliographies can be autogenerated in R Markdown in a similar
way as in Latex/Bibtex. Reference collections should be stored in a separate
file in Bibtex or other supported formats. To cite a publication in an R Markdown
script, one uses the syntax `[@<id1>]` where `<id1>` needs to be replaced with a
reference identifier present in the Bibtex database listed in the metadata section
of the R Markdown script (*e.g.* `bibtex.bib`). For instance, to cite Lawrence et al.
(2013), one uses its reference identifier (*e.g.* `Lawrence2013-kt`) as `<id1>` (Lawrence et al. 2013).
This will place the citation inline in the text and add the corresponding
reference to a reference list at the end of the output document. For the latter a
special section called `References` needs to be specified at the end of the R Markdown script.
To fine control the formatting of citations and reference lists, users want to consult this
the corresponding [R Markdown page](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html).
Also, for general reference management and outputting references in Bibtex format [Paperpile](https://paperpile.com/features)
can be very helpful.

#### Viewing R Markdown report on HPCC cluster

R Markdown reports located on UCR’s HPCC Cluster can be viewed locally in a web browser (without moving
the source HTML) by creating a symbolic link from a user’s `.html` directory. This way any updates to
the report can be viewed immediately without creating another copy of the HTML file. For instance, if user
`ttest` has generated an R Markdown report under `~/bigdata/today/rmarkdown/sample.html`, then the
proper symbolic link to this file can be created as follows:

``` r
cd ~/.html
ln -s ~/bigdata/today/rmarkdown/sample.html sample.html
```

After this one can view the report in a web browser using this URL <http://biocluster.ucr.edu/~ttest/rmarkdown/sample.html>.
If necessary access to the URL can be restricted with a password following the instructions [here](http://hpcc.ucr.edu/manuals_linux-cluster_sharing.html#sharing-files-on-the-web).

A sample R Markdown report for an RNA-Seq project is given here:

  - [RNASeq.html](https://girke.bioinformatics.ucr.edu/GEN242/pages/mydoc/demo_files/systemPipeRNAseq.html){:target="\_blank"}
  - [RNASeq.Rmd](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/pages/mydoc/demo_files/systemPipeRNAseq.Rmd){:target="\_blank"}

## Shiny Web Apps

### What is Shiny?

[Shiny](https://shiny.rstudio.com/gallery/) is an R-based environment for building interactive web applications for
data analysis and exploration. Since most JavaScript code is autogenerated by
the environment, basic R knowledge is sufficient for developing Shiny apps.
They can be deployed on local computers or web servers including custom and cloud-based servers (
e.g. AWS, GCP, [shinyapp.io](http://www.shinyapps.io/) service). The basic structure of a Shiny app is an
`app.R` script containing the following components:

1.  User interface
    
    ``` r
    ui <- fluidPage()
    ```

2.  Server function
    
    ``` r
    server <- function(input, output) {}
    ```

3.  Statement to run shiny app
    
    ``` r
    shinyApp(ui = ui, server = server)
    ```

Alternatively, the `ui` and `server` functions can be organized in two script files, a `ui.R` and a `server.R` script, respectively.

### Develop and test Shiny app locally

Open R and set session to parent directory (here `myappdir`) containing shiny script `app.R`, and the
run it with the `runApp()` function. A sample `app.R` script for testing can be downloaded from [here](https://raw.githubusercontent.com/tgirke/GEN242/gh-pages/_vignettes/07_Rbasics/shinyapp/app.R).

``` r
library(shiny)
runApp("myappdir") # To show code in app, add argument: display.mode="showcase" 
```

This will open the app in a web browser.

### Deploy on web server

This can be done on local or cloud systems. An easy solution is to get an account on [shinyapps.io](http://www.shinyapps.io/)
and then deploy Shiny apps there. For details, see [here](https://shiny.rstudio.com/deploy/).

``` r
setwd("myappdir")
library(rsconnect)
deployApp()
```

### Example Shiny app

The following Shiny app is hosted on `shinyapps.io` and embedded into the markdown (or html) source of this page
using the following iframe syntax:

``` r
<iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px"></iframe>
```

<iframe src="https://tgirke.shinyapps.io/diamonds/" style="border: none; width: 880px; height: 900px">

</iframe>

### Learning Shiny

The Shiny section on the Rstudio site contains excellent [tutorials](https://shiny.rstudio.com/tutorial/).
In addition, users may want to explore the example apps included in the `shiny` package. This can be
done by loading the individual examples (see [here](https://shiny.rstudio.com/tutorial/written-tutorial/lesson1/)) or saving
the code to a user writable directory like this:

``` r
mydir <- system.file("examples", package="shiny")
dir.create('my_shiny_test_dir')
file.copy(mydir, "my_shiny_test_dir", recursive=TRUE)
setwd("my_shiny_test_dir/examples")
runApp("01_hello") # Runs first example app in directory 
dir() # Lists available Shiny examples (directories). 
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

<div id="refs" class="references">

<div id="ref-Lawrence2013-kt">

Lawrence, Michael, Wolfgang Huber, Hervé Pagès, Patrick Aboyoun, Marc Carlson, Robert Gentleman, Martin T Morgan, and Vincent J Carey. 2013. “Software for Computing and Annotating Genomic Ranges.” *PLoS Comput. Biol.* 9 (8): e1003118. <https://doi.org/10.1371/journal.pcbi.1003118>.

</div>

</div>
