---
title: "Environments dplyr, tidyr and some SQLite"
author: Thomas Girke
date: May 25, 2021
output: 
  ioslides_presentation:
    keep_md: yes
    logo: ./images/ucr_logo.png
    widescreen: yes
    df_print: paged
    smaller: true
subtitle: "Introduction to dplyr, tidyr, readr and some SQLite" 
bibliography: bibtex.bib
---

<!---
- ioslides manual: 
   https://bookdown.org/yihui/rmarkdown/ioslides-presentation.html

- Compile from command-line
Rscript -e "rmarkdown::render('dplyr_slides.Rmd'); knitr::knit('dplyr_slides.Rmd', tangle=TRUE)"
-->

<!---
  Note: following css chunks are required for scrolling support beyond slide boundaries
-->

<style>
slides > slide {
  overflow-x: auto !important;
  overflow-y: auto !important;
}
</style>

<style type="text/css">
pre {
  max-height: 300px;
  overflow-y: auto;
}

pre[class] {
  max-height: 300px;
}
</style>

<style type="text/css">
.scroll-300 {
  max-height: 300px;
  overflow-y: auto;
  background-color: inherit;
}
</style>

## How to Navigate this Slide Show?

<br/>

- This __ioslides__ presentation contains scrollable slides. 
- Which slides are scrollable, is indicated by a tag at the bottom of the corresponding slides stating: 

<p style='text-align: center;'> __[ Scroll down to continue ]__ </p>

- The following single character keyboard shortcuts enable alternate display modes of __ioslides__:
    - `f`: enable fullscreen mode
    - `w`: toggle widescreen mode
    - `o`: enable overview mode
    - `h`: enable code highlight mode
- Pressing Esc exits all of these modes. Additional details can be found [here](https://bookdown.org/yihui/rmarkdown/ioslides-presentation.html).

# Outline

- <div class="white">__Overview__</div>
- Install
- File Import and Export
- Usage
- Chaining (Pipes)
- SQLite Databases
- References


## Overview

Modern object classes and methods for handling `data.frame` like structures
are provided by the `dplyr` (`tidyr`) and `data.table` packages. A related example is Bioconductor's 
`DataTable` object class. This tutorial provide a short introduction to the usage and 
functionalities of the `dplyr` and related packages.  

### Related documentation 

More detailed tutorials on this topic can be found here:

* [dplyr: A Grammar of Data Manipulation](https://rdrr.io/cran/dplyr/)
* [Introduction to `dplyr`](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)
* [Tutorial on `dplyr`](http://genomicsclass.github.io/book/pages/dplyr_tutorial.html)
* [Cheatsheet for Joins from Jenny Bryan](http://stat545.com/bit001_dplyr-cheatsheet.html)
* [Tibbles](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html)
* [Intro to `data.table` package](https://www.r-bloggers.com/intro-to-the-data-table-package/)
* [Big data with `dplyr` and `data.table`](https://www.r-bloggers.com/working-with-large-datasets-with-dplyr-and-data-table/)
* [Fast lookups with `dplyr` and `data.table`](https://www.r-bloggers.com/fast-data-lookups-in-r-dplyr-vs-data-table/)

# Outline

- Overview
- <div class="white">__Install__</div>
- File Import and Export
- Usage
- Chaining (Pipes)
- SQLite
- References

## Installation

The `dplyr` (`tidyr`) environment has evolved into an ecosystem of packages. To simplify
package management, one can install and load the entire collection via the
`tidyverse` package [@noauthor_undated-kc]. For more details on `tidyverse` see
[here](http://tidyverse.org/).



```r
install.packages("tidyverse")
```

Load `tidyverse` package environment


```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✔ ggplot2 3.3.3     ✔ purrr   0.3.4
## ✔ tibble  3.1.2     ✔ dplyr   1.0.6
## ✔ tidyr   1.1.3     ✔ stringr 1.4.0
## ✔ readr   1.4.0     ✔ forcats 0.5.1
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

# Outline

- Overview
- Install
- <div class="white">__File Import and Export__</div>
- Usage
- Chaining (Pipes)
- SQLite Databases
- References

## Construct objects

### Construct a `tibble`


```r
as_tibble(iris) # coerce data.frame to tibble tbl
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["fct"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Alternative function producing the same result include `tbl_df` and `as_data_frame` (latter
has been deprecated):


```r
tbl_df(iris) 
```
## Reading and writing tabular files

While the base R read/write utilities can be used for `data.frames`, best time
performance with the least amount of typing is achieved with the export/import
functions from the `readr` package. For very large files the `fread` function from 
the `data.table` package achieves the best time performance. 


### Import with `readr` 

Import functions provided by `readr` include:

* `read_csv()`: comma separated (CSV) files
* `read_tsv()`: tab separated files
* `read_delim()`: general delimited files
* `read_fwf()`: fixed width files
* `read_table()`: tabular files where colums are separated by white-space.
* `read_log()`: web log files

<p style='text-align: right;'> __[ scroll down to continue ]__ </p>
<br/><br/>
<br/><br/>


Create a sample tab delimited file for import


```r
write_tsv(iris, "iris.txt") # Creates sample file
```

Import with `read_tsv` 


```r
iris_df <- read_tsv("iris.txt") # Import with read_tbv from readr package
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   Sepal.Length = col_double(),
##   Sepal.Width = col_double(),
##   Petal.Length = col_double(),
##   Petal.Width = col_double(),
##   Species = col_character()
## )
```

```r
iris_df
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

To import Google Sheets directly into R, see [here](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/#reading-and-writing-external-data).

### Fast table import with `fread` 

The `fread` function from the `data.table` package provides the best time performance for reading large
tabular files into R.


```r
library(data.table)
```

```
## 
## Attaching package: 'data.table'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
```

```
## The following object is masked from 'package:purrr':
## 
##     transpose
```

```r
iris_df <- as_data_frame(fread("iris.txt")) # Import with fread and conversion to tibble
```

```
## Warning: `as_data_frame()` was deprecated in tibble 2.0.0.
## Please use `as_tibble()` instead.
## The signature and semantics have changed, see `?as_tibble`.
```

```r
iris_df
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Note: to ignore lines starting with comment signs, one can pass on to `fread` a shell
command for preprocessing the file. The following example illustrates this option.


```r
fread("grep -v '^#' iris.txt") 
```

### Export with `readr` 

Export function provided by `readr` inlcude

* `write_delim()`: general delimited files
* `write_csv()`: comma separated (CSV) files 
* `write_excel_csv()`: excel style CSV files
* `write_tsv()`: tab separated files

For instance, the `write_tsv` function writes a `data.frame` or `tibble` to a tab delimited file with much nicer
default settings than the base R `write.table` function. 


```r
write_tsv(iris_df, "iris.txt")
```

# Outline

- Overview
- Install
- File Import and Export
- <div class="white">__Usage__</div>
- Chaining (Pipes)
- SQLite Databases
- References

## How to use `tibbles`?

### Column and row binds

The equivalents to base R's `rbind` and `cbind` are `bind_rows` and `bind_cols`, respectively.


```r
bind_cols(iris_df, iris_df)[1:2,]
```

```
## New names:
## * Sepal.Length -> Sepal.Length...1
## * Sepal.Width -> Sepal.Width...2
## * Petal.Length -> Petal.Length...3
## * Petal.Width -> Petal.Width...4
## * Species -> Species...5
## * ...
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length...1"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width...2"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length...3"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width...4"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species...5"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Sepal.Length...6"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width...7"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["Petal.Length...8"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Petal.Width...9"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["Species...10"],"name":[10],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa","6":"5.1","7":"3.5","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa","6":"4.9","7":"3.0","8":"1.4","9":"0.2","10":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

<p style='text-align: right;'> __[ scroll down to continue ]__ </p>
<br/><br/>
<br/><br/>


```r
bind_rows(iris_df, iris_df)[1:2,]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Extract column as vector

The subsetting operators `[[` and `$`can be used to extract from a `tibble` single columns as vector.


```r
iris_df[[5]][1:12]
```

```
##  [1] "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa"
##  [9] "setosa" "setosa" "setosa" "setosa"
```

```r
iris_df$Species[1:12]
```

```
##  [1] "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa" "setosa"
##  [9] "setosa" "setosa" "setosa" "setosa"
```

### Important `dplyr` functions

1. `filter()` and `slice()`
2. `arrange()`
3. `select()` and `rename()`
4. `distinct()`
5. `mutate()` and `transmute()`
6. `summarise()`
7. `sample_n()` and `sample_frac()`


### Slice and filter functions 

#### Filter function


```r
filter(iris_df, Sepal.Length > 7.5, Species=="virginica")
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Base R code equivalent


```r
iris_df[iris_df[, "Sepal.Length"] > 7.5 & iris_df[, "Species"]=="virginica", ]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Including boolean operators


```r
filter(iris_df, Sepal.Length > 7.5 | Sepal.Length < 5.5, Species=="virginica")
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Subset rows by position

`dplyr` approach


```r
slice(iris_df, 1:2)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Base R code equivalent


```r
iris_df[1:2,]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Subset rows by names

Since `tibbles` do not contain row names, row wise subsetting via the `[,]` operator cannot be used.
However, the corresponding behavior can be achieved by passing to `select` a row position index 
obtained by basic R intersect utilities such as `match`.


Create a suitable test `tibble`


```r
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
```

```
## Warning: `data_frame()` was deprecated in tibble 1.1.0.
## Please use `tibble()` instead.
```

```r
df1
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g1","2":"1","3":"11","4":"21","5":"31"},{"1":"g2","2":"2","3":"12","4":"22","5":"32"},{"1":"g3","2":"3","3":"13","4":"23","5":"33"},{"1":"g4","2":"4","3":"14","4":"24","5":"34"},{"1":"g5","2":"5","3":"15","4":"25","5":"35"},{"1":"g6","2":"6","3":"16","4":"26","5":"36"},{"1":"g7","2":"7","3":"17","4":"27","5":"37"},{"1":"g8","2":"8","3":"18","4":"28","5":"38"},{"1":"g9","2":"9","3":"19","4":"29","5":"39"},{"1":"g10","2":"10","3":"20","4":"30","5":"40"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

`dplyr` approach


```r
slice(df1, match(c("g10", "g4", "g4"), ids1))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g10","2":"10","3":"20","4":"30","5":"40"},{"1":"g4","2":"4","3":"14","4":"24","5":"34"},{"1":"g4","2":"4","3":"14","4":"24","5":"34"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Base R equivalent


```r
df1_old <- as.data.frame(df1)
rownames(df1_old) <- df1_old[,1]
df1_old[c("g10", "g4", "g4"),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g10","2":"10","3":"20","4":"30","5":"40","_rn_":"g10"},{"1":"g4","2":"4","3":"14","4":"24","5":"34","_rn_":"g4"},{"1":"g4","2":"4","3":"14","4":"24","5":"34","_rn_":"g4.1"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Sorting with `arrange`

Row-wise ordering based on specific columns

`dplyr` approach


```r
arrange(iris_df, Species, Sepal.Length, Sepal.Width)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

For ordering descendingly use `desc()` function


```r
arrange(iris_df, desc(Species), Sepal.Length, Sepal.Width)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Base R code equivalent


```r
iris_df[order(iris_df$Species, iris_df$Sepal.Length, iris_df$Sepal.Width), ]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
iris_df[order(iris_df$Species, decreasing=TRUE), ] 
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Select columns with `select`

Select specific columns


```r
select(iris_df, Species, Petal.Length, Sepal.Length)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Species"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Petal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Sepal.Length"],"name":[3],"type":["dbl"],"align":["right"]}],"data":[{"1":"setosa","2":"1.4","3":"5.1"},{"1":"setosa","2":"1.4","3":"4.9"},{"1":"setosa","2":"1.3","3":"4.7"},{"1":"setosa","2":"1.5","3":"4.6"},{"1":"setosa","2":"1.4","3":"5.0"},{"1":"setosa","2":"1.7","3":"5.4"},{"1":"setosa","2":"1.4","3":"4.6"},{"1":"setosa","2":"1.5","3":"5.0"},{"1":"setosa","2":"1.4","3":"4.4"},{"1":"setosa","2":"1.5","3":"4.9"},{"1":"setosa","2":"1.5","3":"5.4"},{"1":"setosa","2":"1.6","3":"4.8"},{"1":"setosa","2":"1.4","3":"4.8"},{"1":"setosa","2":"1.1","3":"4.3"},{"1":"setosa","2":"1.2","3":"5.8"},{"1":"setosa","2":"1.5","3":"5.7"},{"1":"setosa","2":"1.3","3":"5.4"},{"1":"setosa","2":"1.4","3":"5.1"},{"1":"setosa","2":"1.7","3":"5.7"},{"1":"setosa","2":"1.5","3":"5.1"},{"1":"setosa","2":"1.7","3":"5.4"},{"1":"setosa","2":"1.5","3":"5.1"},{"1":"setosa","2":"1.0","3":"4.6"},{"1":"setosa","2":"1.7","3":"5.1"},{"1":"setosa","2":"1.9","3":"4.8"},{"1":"setosa","2":"1.6","3":"5.0"},{"1":"setosa","2":"1.6","3":"5.0"},{"1":"setosa","2":"1.5","3":"5.2"},{"1":"setosa","2":"1.4","3":"5.2"},{"1":"setosa","2":"1.6","3":"4.7"},{"1":"setosa","2":"1.6","3":"4.8"},{"1":"setosa","2":"1.5","3":"5.4"},{"1":"setosa","2":"1.5","3":"5.2"},{"1":"setosa","2":"1.4","3":"5.5"},{"1":"setosa","2":"1.5","3":"4.9"},{"1":"setosa","2":"1.2","3":"5.0"},{"1":"setosa","2":"1.3","3":"5.5"},{"1":"setosa","2":"1.4","3":"4.9"},{"1":"setosa","2":"1.3","3":"4.4"},{"1":"setosa","2":"1.5","3":"5.1"},{"1":"setosa","2":"1.3","3":"5.0"},{"1":"setosa","2":"1.3","3":"4.5"},{"1":"setosa","2":"1.3","3":"4.4"},{"1":"setosa","2":"1.6","3":"5.0"},{"1":"setosa","2":"1.9","3":"5.1"},{"1":"setosa","2":"1.4","3":"4.8"},{"1":"setosa","2":"1.6","3":"5.1"},{"1":"setosa","2":"1.4","3":"4.6"},{"1":"setosa","2":"1.5","3":"5.3"},{"1":"setosa","2":"1.4","3":"5.0"},{"1":"versicolor","2":"4.7","3":"7.0"},{"1":"versicolor","2":"4.5","3":"6.4"},{"1":"versicolor","2":"4.9","3":"6.9"},{"1":"versicolor","2":"4.0","3":"5.5"},{"1":"versicolor","2":"4.6","3":"6.5"},{"1":"versicolor","2":"4.5","3":"5.7"},{"1":"versicolor","2":"4.7","3":"6.3"},{"1":"versicolor","2":"3.3","3":"4.9"},{"1":"versicolor","2":"4.6","3":"6.6"},{"1":"versicolor","2":"3.9","3":"5.2"},{"1":"versicolor","2":"3.5","3":"5.0"},{"1":"versicolor","2":"4.2","3":"5.9"},{"1":"versicolor","2":"4.0","3":"6.0"},{"1":"versicolor","2":"4.7","3":"6.1"},{"1":"versicolor","2":"3.6","3":"5.6"},{"1":"versicolor","2":"4.4","3":"6.7"},{"1":"versicolor","2":"4.5","3":"5.6"},{"1":"versicolor","2":"4.1","3":"5.8"},{"1":"versicolor","2":"4.5","3":"6.2"},{"1":"versicolor","2":"3.9","3":"5.6"},{"1":"versicolor","2":"4.8","3":"5.9"},{"1":"versicolor","2":"4.0","3":"6.1"},{"1":"versicolor","2":"4.9","3":"6.3"},{"1":"versicolor","2":"4.7","3":"6.1"},{"1":"versicolor","2":"4.3","3":"6.4"},{"1":"versicolor","2":"4.4","3":"6.6"},{"1":"versicolor","2":"4.8","3":"6.8"},{"1":"versicolor","2":"5.0","3":"6.7"},{"1":"versicolor","2":"4.5","3":"6.0"},{"1":"versicolor","2":"3.5","3":"5.7"},{"1":"versicolor","2":"3.8","3":"5.5"},{"1":"versicolor","2":"3.7","3":"5.5"},{"1":"versicolor","2":"3.9","3":"5.8"},{"1":"versicolor","2":"5.1","3":"6.0"},{"1":"versicolor","2":"4.5","3":"5.4"},{"1":"versicolor","2":"4.5","3":"6.0"},{"1":"versicolor","2":"4.7","3":"6.7"},{"1":"versicolor","2":"4.4","3":"6.3"},{"1":"versicolor","2":"4.1","3":"5.6"},{"1":"versicolor","2":"4.0","3":"5.5"},{"1":"versicolor","2":"4.4","3":"5.5"},{"1":"versicolor","2":"4.6","3":"6.1"},{"1":"versicolor","2":"4.0","3":"5.8"},{"1":"versicolor","2":"3.3","3":"5.0"},{"1":"versicolor","2":"4.2","3":"5.6"},{"1":"versicolor","2":"4.2","3":"5.7"},{"1":"versicolor","2":"4.2","3":"5.7"},{"1":"versicolor","2":"4.3","3":"6.2"},{"1":"versicolor","2":"3.0","3":"5.1"},{"1":"versicolor","2":"4.1","3":"5.7"},{"1":"virginica","2":"6.0","3":"6.3"},{"1":"virginica","2":"5.1","3":"5.8"},{"1":"virginica","2":"5.9","3":"7.1"},{"1":"virginica","2":"5.6","3":"6.3"},{"1":"virginica","2":"5.8","3":"6.5"},{"1":"virginica","2":"6.6","3":"7.6"},{"1":"virginica","2":"4.5","3":"4.9"},{"1":"virginica","2":"6.3","3":"7.3"},{"1":"virginica","2":"5.8","3":"6.7"},{"1":"virginica","2":"6.1","3":"7.2"},{"1":"virginica","2":"5.1","3":"6.5"},{"1":"virginica","2":"5.3","3":"6.4"},{"1":"virginica","2":"5.5","3":"6.8"},{"1":"virginica","2":"5.0","3":"5.7"},{"1":"virginica","2":"5.1","3":"5.8"},{"1":"virginica","2":"5.3","3":"6.4"},{"1":"virginica","2":"5.5","3":"6.5"},{"1":"virginica","2":"6.7","3":"7.7"},{"1":"virginica","2":"6.9","3":"7.7"},{"1":"virginica","2":"5.0","3":"6.0"},{"1":"virginica","2":"5.7","3":"6.9"},{"1":"virginica","2":"4.9","3":"5.6"},{"1":"virginica","2":"6.7","3":"7.7"},{"1":"virginica","2":"4.9","3":"6.3"},{"1":"virginica","2":"5.7","3":"6.7"},{"1":"virginica","2":"6.0","3":"7.2"},{"1":"virginica","2":"4.8","3":"6.2"},{"1":"virginica","2":"4.9","3":"6.1"},{"1":"virginica","2":"5.6","3":"6.4"},{"1":"virginica","2":"5.8","3":"7.2"},{"1":"virginica","2":"6.1","3":"7.4"},{"1":"virginica","2":"6.4","3":"7.9"},{"1":"virginica","2":"5.6","3":"6.4"},{"1":"virginica","2":"5.1","3":"6.3"},{"1":"virginica","2":"5.6","3":"6.1"},{"1":"virginica","2":"6.1","3":"7.7"},{"1":"virginica","2":"5.6","3":"6.3"},{"1":"virginica","2":"5.5","3":"6.4"},{"1":"virginica","2":"4.8","3":"6.0"},{"1":"virginica","2":"5.4","3":"6.9"},{"1":"virginica","2":"5.6","3":"6.7"},{"1":"virginica","2":"5.1","3":"6.9"},{"1":"virginica","2":"5.1","3":"5.8"},{"1":"virginica","2":"5.9","3":"6.8"},{"1":"virginica","2":"5.7","3":"6.7"},{"1":"virginica","2":"5.2","3":"6.7"},{"1":"virginica","2":"5.0","3":"6.3"},{"1":"virginica","2":"5.2","3":"6.5"},{"1":"virginica","2":"5.4","3":"6.2"},{"1":"virginica","2":"5.1","3":"5.9"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Select range of columns by name


```r
select(iris_df, Sepal.Length : Petal.Width)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Drop specific columns (here range)


```r
select(iris_df, -(Sepal.Length : Petal.Width))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Species"],"name":[1],"type":["chr"],"align":["left"]}],"data":[{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"setosa"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"versicolor"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"},{"1":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Renaming columns with `rename`


`dplyr` approach


```r
rename(iris_df, new_col_name = Species)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["new_col_name"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Base R code approach


```r
colnames(iris_df)[colnames(iris_df)=="Species"] <- "new_col_names"
```

### Obtain unique rows with `distinct`

`dplyr` approach


```r
distinct(iris_df, Species, .keep_all=TRUE)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Base R code approach


```r
iris_df[!duplicated(iris_df$Species),]
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Add columns

#### `mutate`

The `mutate` function allows to append columns to existing ones.


```r
mutate(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Ratio"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Sum"],"name":[7],"type":["dbl"],"align":["right"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa","6":"1.457143","7":"8.6"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa","6":"1.633333","7":"7.9"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"1.468750","7":"7.9"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"1.483871","7":"7.7"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa","6":"1.388889","7":"8.6"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa","6":"1.384615","7":"9.3"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa","6":"1.352941","7":"8.0"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"1.470588","7":"8.4"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa","6":"1.517241","7":"7.3"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa","6":"1.580645","7":"8.0"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"1.459459","7":"9.1"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa","6":"1.411765","7":"8.2"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa","6":"1.600000","7":"7.8"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa","6":"1.433333","7":"7.3"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa","6":"1.450000","7":"9.8"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa","6":"1.295455","7":"10.1"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa","6":"1.384615","7":"9.3"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa","6":"1.457143","7":"8.6"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa","6":"1.500000","7":"9.5"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa","6":"1.342105","7":"8.9"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa","6":"1.588235","7":"8.8"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa","6":"1.378378","7":"8.8"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa","6":"1.277778","7":"8.2"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa","6":"1.545455","7":"8.4"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa","6":"1.411765","7":"8.2"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa","6":"1.666667","7":"8.0"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa","6":"1.470588","7":"8.4"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa","6":"1.485714","7":"8.7"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa","6":"1.529412","7":"8.6"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa","6":"1.468750","7":"7.9"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa","6":"1.548387","7":"7.9"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa","6":"1.588235","7":"8.8"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa","6":"1.268293","7":"9.3"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa","6":"1.309524","7":"9.7"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"1.580645","7":"8.0"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa","6":"1.562500","7":"8.2"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa","6":"1.571429","7":"9.0"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa","6":"1.361111","7":"8.5"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa","6":"1.466667","7":"7.4"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"1.500000","7":"8.5"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa","6":"1.428571","7":"8.5"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa","6":"1.956522","7":"6.8"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"1.375000","7":"7.6"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa","6":"1.428571","7":"8.5"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa","6":"1.342105","7":"8.9"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa","6":"1.600000","7":"7.8"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa","6":"1.342105","7":"8.9"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa","6":"1.437500","7":"7.8"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"1.432432","7":"9.0"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa","6":"1.515152","7":"8.3"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor","6":"2.187500","7":"10.2"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor","6":"2.000000","7":"9.6"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor","6":"2.225806","7":"10.0"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor","6":"2.391304","7":"7.8"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor","6":"2.321429","7":"9.3"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor","6":"2.035714","7":"8.5"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor","6":"1.909091","7":"9.6"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor","6":"2.041667","7":"7.3"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor","6":"2.275862","7":"9.5"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor","6":"1.925926","7":"7.9"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor","6":"2.500000","7":"7.0"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor","6":"1.966667","7":"8.9"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor","6":"2.727273","7":"8.2"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor","6":"2.103448","7":"9.0"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor","6":"1.931034","7":"8.5"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor","6":"2.161290","7":"9.8"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor","6":"1.866667","7":"8.6"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor","6":"2.148148","7":"8.5"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor","6":"2.818182","7":"8.4"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor","6":"2.240000","7":"8.1"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor","6":"1.843750","7":"9.1"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor","6":"2.178571","7":"8.9"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor","6":"2.520000","7":"8.8"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor","6":"2.178571","7":"8.9"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor","6":"2.206897","7":"9.3"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor","6":"2.200000","7":"9.6"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor","6":"2.428571","7":"9.6"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor","6":"2.233333","7":"9.7"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor","6":"2.068966","7":"8.9"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor","6":"2.192308","7":"8.3"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor","6":"2.291667","7":"7.9"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor","6":"2.291667","7":"7.9"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor","6":"2.148148","7":"8.5"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor","6":"2.222222","7":"8.7"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor","6":"1.800000","7":"8.4"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor","6":"1.764706","7":"9.4"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor","6":"2.161290","7":"9.8"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor","6":"2.739130","7":"8.6"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor","6":"1.866667","7":"8.6"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor","6":"2.200000","7":"8.0"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor","6":"2.115385","7":"8.1"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor","6":"2.033333","7":"9.1"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor","6":"2.230769","7":"8.4"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor","6":"2.173913","7":"7.3"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor","6":"2.074074","7":"8.3"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor","6":"1.900000","7":"8.7"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor","6":"1.965517","7":"8.6"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor","6":"2.137931","7":"9.1"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor","6":"2.040000","7":"7.6"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor","6":"2.035714","7":"8.5"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica","6":"1.909091","7":"9.6"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica","6":"2.148148","7":"8.5"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica","6":"2.366667","7":"10.1"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica","6":"2.172414","7":"9.2"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica","6":"2.166667","7":"9.5"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica","6":"2.533333","7":"10.6"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica","6":"1.960000","7":"7.4"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica","6":"2.517241","7":"10.2"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica","6":"2.680000","7":"9.2"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica","6":"2.000000","7":"10.8"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica","6":"2.031250","7":"9.7"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica","6":"2.370370","7":"9.1"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica","6":"2.266667","7":"9.8"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica","6":"2.280000","7":"8.2"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica","6":"2.071429","7":"8.6"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica","6":"2.000000","7":"9.6"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica","6":"2.166667","7":"9.5"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica","6":"2.026316","7":"11.5"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica","6":"2.961538","7":"10.3"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica","6":"2.727273","7":"8.2"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica","6":"2.156250","7":"10.1"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica","6":"2.000000","7":"8.4"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica","6":"2.750000","7":"10.5"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica","6":"2.333333","7":"9.0"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica","6":"2.030303","7":"10.0"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica","6":"2.250000","7":"10.4"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica","6":"2.214286","7":"9.0"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica","6":"2.033333","7":"9.1"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica","6":"2.285714","7":"9.2"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica","6":"2.400000","7":"10.2"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica","6":"2.642857","7":"10.2"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica","6":"2.078947","7":"11.7"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica","6":"2.285714","7":"9.2"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica","6":"2.250000","7":"9.1"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica","6":"2.346154","7":"8.7"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica","6":"2.566667","7":"10.7"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica","6":"1.852941","7":"9.7"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica","6":"2.064516","7":"9.5"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica","6":"2.000000","7":"9.0"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica","6":"2.225806","7":"10.0"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica","6":"2.161290","7":"9.8"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica","6":"2.225806","7":"10.0"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica","6":"2.148148","7":"8.5"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica","6":"2.125000","7":"10.0"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica","6":"2.030303","7":"10.0"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica","6":"2.233333","7":"9.7"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica","6":"2.520000","7":"8.8"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica","6":"2.166667","7":"9.5"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica","6":"1.823529","7":"9.6"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica","6":"1.966667","7":"8.9"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### `transmute`

The `transmute` function does the same as `mutate` but drops existing columns


```r
transmute(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Ratio"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sum"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"1.457143","2":"8.6"},{"1":"1.633333","2":"7.9"},{"1":"1.468750","2":"7.9"},{"1":"1.483871","2":"7.7"},{"1":"1.388889","2":"8.6"},{"1":"1.384615","2":"9.3"},{"1":"1.352941","2":"8.0"},{"1":"1.470588","2":"8.4"},{"1":"1.517241","2":"7.3"},{"1":"1.580645","2":"8.0"},{"1":"1.459459","2":"9.1"},{"1":"1.411765","2":"8.2"},{"1":"1.600000","2":"7.8"},{"1":"1.433333","2":"7.3"},{"1":"1.450000","2":"9.8"},{"1":"1.295455","2":"10.1"},{"1":"1.384615","2":"9.3"},{"1":"1.457143","2":"8.6"},{"1":"1.500000","2":"9.5"},{"1":"1.342105","2":"8.9"},{"1":"1.588235","2":"8.8"},{"1":"1.378378","2":"8.8"},{"1":"1.277778","2":"8.2"},{"1":"1.545455","2":"8.4"},{"1":"1.411765","2":"8.2"},{"1":"1.666667","2":"8.0"},{"1":"1.470588","2":"8.4"},{"1":"1.485714","2":"8.7"},{"1":"1.529412","2":"8.6"},{"1":"1.468750","2":"7.9"},{"1":"1.548387","2":"7.9"},{"1":"1.588235","2":"8.8"},{"1":"1.268293","2":"9.3"},{"1":"1.309524","2":"9.7"},{"1":"1.580645","2":"8.0"},{"1":"1.562500","2":"8.2"},{"1":"1.571429","2":"9.0"},{"1":"1.361111","2":"8.5"},{"1":"1.466667","2":"7.4"},{"1":"1.500000","2":"8.5"},{"1":"1.428571","2":"8.5"},{"1":"1.956522","2":"6.8"},{"1":"1.375000","2":"7.6"},{"1":"1.428571","2":"8.5"},{"1":"1.342105","2":"8.9"},{"1":"1.600000","2":"7.8"},{"1":"1.342105","2":"8.9"},{"1":"1.437500","2":"7.8"},{"1":"1.432432","2":"9.0"},{"1":"1.515152","2":"8.3"},{"1":"2.187500","2":"10.2"},{"1":"2.000000","2":"9.6"},{"1":"2.225806","2":"10.0"},{"1":"2.391304","2":"7.8"},{"1":"2.321429","2":"9.3"},{"1":"2.035714","2":"8.5"},{"1":"1.909091","2":"9.6"},{"1":"2.041667","2":"7.3"},{"1":"2.275862","2":"9.5"},{"1":"1.925926","2":"7.9"},{"1":"2.500000","2":"7.0"},{"1":"1.966667","2":"8.9"},{"1":"2.727273","2":"8.2"},{"1":"2.103448","2":"9.0"},{"1":"1.931034","2":"8.5"},{"1":"2.161290","2":"9.8"},{"1":"1.866667","2":"8.6"},{"1":"2.148148","2":"8.5"},{"1":"2.818182","2":"8.4"},{"1":"2.240000","2":"8.1"},{"1":"1.843750","2":"9.1"},{"1":"2.178571","2":"8.9"},{"1":"2.520000","2":"8.8"},{"1":"2.178571","2":"8.9"},{"1":"2.206897","2":"9.3"},{"1":"2.200000","2":"9.6"},{"1":"2.428571","2":"9.6"},{"1":"2.233333","2":"9.7"},{"1":"2.068966","2":"8.9"},{"1":"2.192308","2":"8.3"},{"1":"2.291667","2":"7.9"},{"1":"2.291667","2":"7.9"},{"1":"2.148148","2":"8.5"},{"1":"2.222222","2":"8.7"},{"1":"1.800000","2":"8.4"},{"1":"1.764706","2":"9.4"},{"1":"2.161290","2":"9.8"},{"1":"2.739130","2":"8.6"},{"1":"1.866667","2":"8.6"},{"1":"2.200000","2":"8.0"},{"1":"2.115385","2":"8.1"},{"1":"2.033333","2":"9.1"},{"1":"2.230769","2":"8.4"},{"1":"2.173913","2":"7.3"},{"1":"2.074074","2":"8.3"},{"1":"1.900000","2":"8.7"},{"1":"1.965517","2":"8.6"},{"1":"2.137931","2":"9.1"},{"1":"2.040000","2":"7.6"},{"1":"2.035714","2":"8.5"},{"1":"1.909091","2":"9.6"},{"1":"2.148148","2":"8.5"},{"1":"2.366667","2":"10.1"},{"1":"2.172414","2":"9.2"},{"1":"2.166667","2":"9.5"},{"1":"2.533333","2":"10.6"},{"1":"1.960000","2":"7.4"},{"1":"2.517241","2":"10.2"},{"1":"2.680000","2":"9.2"},{"1":"2.000000","2":"10.8"},{"1":"2.031250","2":"9.7"},{"1":"2.370370","2":"9.1"},{"1":"2.266667","2":"9.8"},{"1":"2.280000","2":"8.2"},{"1":"2.071429","2":"8.6"},{"1":"2.000000","2":"9.6"},{"1":"2.166667","2":"9.5"},{"1":"2.026316","2":"11.5"},{"1":"2.961538","2":"10.3"},{"1":"2.727273","2":"8.2"},{"1":"2.156250","2":"10.1"},{"1":"2.000000","2":"8.4"},{"1":"2.750000","2":"10.5"},{"1":"2.333333","2":"9.0"},{"1":"2.030303","2":"10.0"},{"1":"2.250000","2":"10.4"},{"1":"2.214286","2":"9.0"},{"1":"2.033333","2":"9.1"},{"1":"2.285714","2":"9.2"},{"1":"2.400000","2":"10.2"},{"1":"2.642857","2":"10.2"},{"1":"2.078947","2":"11.7"},{"1":"2.285714","2":"9.2"},{"1":"2.250000","2":"9.1"},{"1":"2.346154","2":"8.7"},{"1":"2.566667","2":"10.7"},{"1":"1.852941","2":"9.7"},{"1":"2.064516","2":"9.5"},{"1":"2.000000","2":"9.0"},{"1":"2.225806","2":"10.0"},{"1":"2.161290","2":"9.8"},{"1":"2.225806","2":"10.0"},{"1":"2.148148","2":"8.5"},{"1":"2.125000","2":"10.0"},{"1":"2.030303","2":"10.0"},{"1":"2.233333","2":"9.7"},{"1":"2.520000","2":"8.8"},{"1":"2.166667","2":"9.5"},{"1":"1.823529","2":"9.6"},{"1":"1.966667","2":"8.9"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### `bind_cols`

The `bind_cols` function is the equivalent of `cbind` in base R. To add rows, use the corresponding 
`bind_rows` function.


```r
bind_cols(iris_df, iris_df)
```

```
## New names:
## * Sepal.Length -> Sepal.Length...1
## * Sepal.Width -> Sepal.Width...2
## * Petal.Length -> Petal.Length...3
## * Petal.Width -> Petal.Width...4
## * Species -> Species...5
## * ...
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length...1"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width...2"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length...3"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width...4"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species...5"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Sepal.Length...6"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width...7"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["Petal.Length...8"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Petal.Width...9"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["Species...10"],"name":[10],"type":["chr"],"align":["left"]}],"data":[{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa","6":"5.1","7":"3.5","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa","6":"4.9","7":"3.0","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"4.7","7":"3.2","8":"1.3","9":"0.2","10":"setosa"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"4.6","7":"3.1","8":"1.5","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa","6":"5.0","7":"3.6","8":"1.4","9":"0.2","10":"setosa"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa","6":"5.4","7":"3.9","8":"1.7","9":"0.4","10":"setosa"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa","6":"4.6","7":"3.4","8":"1.4","9":"0.3","10":"setosa"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"5.0","7":"3.4","8":"1.5","9":"0.2","10":"setosa"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa","6":"4.4","7":"2.9","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa","6":"4.9","7":"3.1","8":"1.5","9":"0.1","10":"setosa"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"5.4","7":"3.7","8":"1.5","9":"0.2","10":"setosa"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa","6":"4.8","7":"3.4","8":"1.6","9":"0.2","10":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa","6":"4.8","7":"3.0","8":"1.4","9":"0.1","10":"setosa"},{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa","6":"4.3","7":"3.0","8":"1.1","9":"0.1","10":"setosa"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa","6":"5.8","7":"4.0","8":"1.2","9":"0.2","10":"setosa"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa","6":"5.7","7":"4.4","8":"1.5","9":"0.4","10":"setosa"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa","6":"5.4","7":"3.9","8":"1.3","9":"0.4","10":"setosa"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa","6":"5.1","7":"3.5","8":"1.4","9":"0.3","10":"setosa"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa","6":"5.7","7":"3.8","8":"1.7","9":"0.3","10":"setosa"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa","6":"5.1","7":"3.8","8":"1.5","9":"0.3","10":"setosa"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa","6":"5.4","7":"3.4","8":"1.7","9":"0.2","10":"setosa"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa","6":"5.1","7":"3.7","8":"1.5","9":"0.4","10":"setosa"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa","6":"4.6","7":"3.6","8":"1.0","9":"0.2","10":"setosa"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa","6":"5.1","7":"3.3","8":"1.7","9":"0.5","10":"setosa"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa","6":"4.8","7":"3.4","8":"1.9","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa","6":"5.0","7":"3.0","8":"1.6","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa","6":"5.0","7":"3.4","8":"1.6","9":"0.4","10":"setosa"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa","6":"5.2","7":"3.5","8":"1.5","9":"0.2","10":"setosa"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa","6":"5.2","7":"3.4","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa","6":"4.7","7":"3.2","8":"1.6","9":"0.2","10":"setosa"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa","6":"4.8","7":"3.1","8":"1.6","9":"0.2","10":"setosa"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa","6":"5.4","7":"3.4","8":"1.5","9":"0.4","10":"setosa"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa","6":"5.2","7":"4.1","8":"1.5","9":"0.1","10":"setosa"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa","6":"5.5","7":"4.2","8":"1.4","9":"0.2","10":"setosa"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"4.9","7":"3.1","8":"1.5","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa","6":"5.0","7":"3.2","8":"1.2","9":"0.2","10":"setosa"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa","6":"5.5","7":"3.5","8":"1.3","9":"0.2","10":"setosa"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa","6":"4.9","7":"3.6","8":"1.4","9":"0.1","10":"setosa"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa","6":"4.4","7":"3.0","8":"1.3","9":"0.2","10":"setosa"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"5.1","7":"3.4","8":"1.5","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa","6":"5.0","7":"3.5","8":"1.3","9":"0.3","10":"setosa"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa","6":"4.5","7":"2.3","8":"1.3","9":"0.3","10":"setosa"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"4.4","7":"3.2","8":"1.3","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa","6":"5.0","7":"3.5","8":"1.6","9":"0.6","10":"setosa"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa","6":"5.1","7":"3.8","8":"1.9","9":"0.4","10":"setosa"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa","6":"4.8","7":"3.0","8":"1.4","9":"0.3","10":"setosa"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa","6":"5.1","7":"3.8","8":"1.6","9":"0.2","10":"setosa"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa","6":"4.6","7":"3.2","8":"1.4","9":"0.2","10":"setosa"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"5.3","7":"3.7","8":"1.5","9":"0.2","10":"setosa"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa","6":"5.0","7":"3.3","8":"1.4","9":"0.2","10":"setosa"},{"1":"7.0","2":"3.2","3":"4.7","4":"1.4","5":"versicolor","6":"7.0","7":"3.2","8":"4.7","9":"1.4","10":"versicolor"},{"1":"6.4","2":"3.2","3":"4.5","4":"1.5","5":"versicolor","6":"6.4","7":"3.2","8":"4.5","9":"1.5","10":"versicolor"},{"1":"6.9","2":"3.1","3":"4.9","4":"1.5","5":"versicolor","6":"6.9","7":"3.1","8":"4.9","9":"1.5","10":"versicolor"},{"1":"5.5","2":"2.3","3":"4.0","4":"1.3","5":"versicolor","6":"5.5","7":"2.3","8":"4.0","9":"1.3","10":"versicolor"},{"1":"6.5","2":"2.8","3":"4.6","4":"1.5","5":"versicolor","6":"6.5","7":"2.8","8":"4.6","9":"1.5","10":"versicolor"},{"1":"5.7","2":"2.8","3":"4.5","4":"1.3","5":"versicolor","6":"5.7","7":"2.8","8":"4.5","9":"1.3","10":"versicolor"},{"1":"6.3","2":"3.3","3":"4.7","4":"1.6","5":"versicolor","6":"6.3","7":"3.3","8":"4.7","9":"1.6","10":"versicolor"},{"1":"4.9","2":"2.4","3":"3.3","4":"1.0","5":"versicolor","6":"4.9","7":"2.4","8":"3.3","9":"1.0","10":"versicolor"},{"1":"6.6","2":"2.9","3":"4.6","4":"1.3","5":"versicolor","6":"6.6","7":"2.9","8":"4.6","9":"1.3","10":"versicolor"},{"1":"5.2","2":"2.7","3":"3.9","4":"1.4","5":"versicolor","6":"5.2","7":"2.7","8":"3.9","9":"1.4","10":"versicolor"},{"1":"5.0","2":"2.0","3":"3.5","4":"1.0","5":"versicolor","6":"5.0","7":"2.0","8":"3.5","9":"1.0","10":"versicolor"},{"1":"5.9","2":"3.0","3":"4.2","4":"1.5","5":"versicolor","6":"5.9","7":"3.0","8":"4.2","9":"1.5","10":"versicolor"},{"1":"6.0","2":"2.2","3":"4.0","4":"1.0","5":"versicolor","6":"6.0","7":"2.2","8":"4.0","9":"1.0","10":"versicolor"},{"1":"6.1","2":"2.9","3":"4.7","4":"1.4","5":"versicolor","6":"6.1","7":"2.9","8":"4.7","9":"1.4","10":"versicolor"},{"1":"5.6","2":"2.9","3":"3.6","4":"1.3","5":"versicolor","6":"5.6","7":"2.9","8":"3.6","9":"1.3","10":"versicolor"},{"1":"6.7","2":"3.1","3":"4.4","4":"1.4","5":"versicolor","6":"6.7","7":"3.1","8":"4.4","9":"1.4","10":"versicolor"},{"1":"5.6","2":"3.0","3":"4.5","4":"1.5","5":"versicolor","6":"5.6","7":"3.0","8":"4.5","9":"1.5","10":"versicolor"},{"1":"5.8","2":"2.7","3":"4.1","4":"1.0","5":"versicolor","6":"5.8","7":"2.7","8":"4.1","9":"1.0","10":"versicolor"},{"1":"6.2","2":"2.2","3":"4.5","4":"1.5","5":"versicolor","6":"6.2","7":"2.2","8":"4.5","9":"1.5","10":"versicolor"},{"1":"5.6","2":"2.5","3":"3.9","4":"1.1","5":"versicolor","6":"5.6","7":"2.5","8":"3.9","9":"1.1","10":"versicolor"},{"1":"5.9","2":"3.2","3":"4.8","4":"1.8","5":"versicolor","6":"5.9","7":"3.2","8":"4.8","9":"1.8","10":"versicolor"},{"1":"6.1","2":"2.8","3":"4.0","4":"1.3","5":"versicolor","6":"6.1","7":"2.8","8":"4.0","9":"1.3","10":"versicolor"},{"1":"6.3","2":"2.5","3":"4.9","4":"1.5","5":"versicolor","6":"6.3","7":"2.5","8":"4.9","9":"1.5","10":"versicolor"},{"1":"6.1","2":"2.8","3":"4.7","4":"1.2","5":"versicolor","6":"6.1","7":"2.8","8":"4.7","9":"1.2","10":"versicolor"},{"1":"6.4","2":"2.9","3":"4.3","4":"1.3","5":"versicolor","6":"6.4","7":"2.9","8":"4.3","9":"1.3","10":"versicolor"},{"1":"6.6","2":"3.0","3":"4.4","4":"1.4","5":"versicolor","6":"6.6","7":"3.0","8":"4.4","9":"1.4","10":"versicolor"},{"1":"6.8","2":"2.8","3":"4.8","4":"1.4","5":"versicolor","6":"6.8","7":"2.8","8":"4.8","9":"1.4","10":"versicolor"},{"1":"6.7","2":"3.0","3":"5.0","4":"1.7","5":"versicolor","6":"6.7","7":"3.0","8":"5.0","9":"1.7","10":"versicolor"},{"1":"6.0","2":"2.9","3":"4.5","4":"1.5","5":"versicolor","6":"6.0","7":"2.9","8":"4.5","9":"1.5","10":"versicolor"},{"1":"5.7","2":"2.6","3":"3.5","4":"1.0","5":"versicolor","6":"5.7","7":"2.6","8":"3.5","9":"1.0","10":"versicolor"},{"1":"5.5","2":"2.4","3":"3.8","4":"1.1","5":"versicolor","6":"5.5","7":"2.4","8":"3.8","9":"1.1","10":"versicolor"},{"1":"5.5","2":"2.4","3":"3.7","4":"1.0","5":"versicolor","6":"5.5","7":"2.4","8":"3.7","9":"1.0","10":"versicolor"},{"1":"5.8","2":"2.7","3":"3.9","4":"1.2","5":"versicolor","6":"5.8","7":"2.7","8":"3.9","9":"1.2","10":"versicolor"},{"1":"6.0","2":"2.7","3":"5.1","4":"1.6","5":"versicolor","6":"6.0","7":"2.7","8":"5.1","9":"1.6","10":"versicolor"},{"1":"5.4","2":"3.0","3":"4.5","4":"1.5","5":"versicolor","6":"5.4","7":"3.0","8":"4.5","9":"1.5","10":"versicolor"},{"1":"6.0","2":"3.4","3":"4.5","4":"1.6","5":"versicolor","6":"6.0","7":"3.4","8":"4.5","9":"1.6","10":"versicolor"},{"1":"6.7","2":"3.1","3":"4.7","4":"1.5","5":"versicolor","6":"6.7","7":"3.1","8":"4.7","9":"1.5","10":"versicolor"},{"1":"6.3","2":"2.3","3":"4.4","4":"1.3","5":"versicolor","6":"6.3","7":"2.3","8":"4.4","9":"1.3","10":"versicolor"},{"1":"5.6","2":"3.0","3":"4.1","4":"1.3","5":"versicolor","6":"5.6","7":"3.0","8":"4.1","9":"1.3","10":"versicolor"},{"1":"5.5","2":"2.5","3":"4.0","4":"1.3","5":"versicolor","6":"5.5","7":"2.5","8":"4.0","9":"1.3","10":"versicolor"},{"1":"5.5","2":"2.6","3":"4.4","4":"1.2","5":"versicolor","6":"5.5","7":"2.6","8":"4.4","9":"1.2","10":"versicolor"},{"1":"6.1","2":"3.0","3":"4.6","4":"1.4","5":"versicolor","6":"6.1","7":"3.0","8":"4.6","9":"1.4","10":"versicolor"},{"1":"5.8","2":"2.6","3":"4.0","4":"1.2","5":"versicolor","6":"5.8","7":"2.6","8":"4.0","9":"1.2","10":"versicolor"},{"1":"5.0","2":"2.3","3":"3.3","4":"1.0","5":"versicolor","6":"5.0","7":"2.3","8":"3.3","9":"1.0","10":"versicolor"},{"1":"5.6","2":"2.7","3":"4.2","4":"1.3","5":"versicolor","6":"5.6","7":"2.7","8":"4.2","9":"1.3","10":"versicolor"},{"1":"5.7","2":"3.0","3":"4.2","4":"1.2","5":"versicolor","6":"5.7","7":"3.0","8":"4.2","9":"1.2","10":"versicolor"},{"1":"5.7","2":"2.9","3":"4.2","4":"1.3","5":"versicolor","6":"5.7","7":"2.9","8":"4.2","9":"1.3","10":"versicolor"},{"1":"6.2","2":"2.9","3":"4.3","4":"1.3","5":"versicolor","6":"6.2","7":"2.9","8":"4.3","9":"1.3","10":"versicolor"},{"1":"5.1","2":"2.5","3":"3.0","4":"1.1","5":"versicolor","6":"5.1","7":"2.5","8":"3.0","9":"1.1","10":"versicolor"},{"1":"5.7","2":"2.8","3":"4.1","4":"1.3","5":"versicolor","6":"5.7","7":"2.8","8":"4.1","9":"1.3","10":"versicolor"},{"1":"6.3","2":"3.3","3":"6.0","4":"2.5","5":"virginica","6":"6.3","7":"3.3","8":"6.0","9":"2.5","10":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica","6":"5.8","7":"2.7","8":"5.1","9":"1.9","10":"virginica"},{"1":"7.1","2":"3.0","3":"5.9","4":"2.1","5":"virginica","6":"7.1","7":"3.0","8":"5.9","9":"2.1","10":"virginica"},{"1":"6.3","2":"2.9","3":"5.6","4":"1.8","5":"virginica","6":"6.3","7":"2.9","8":"5.6","9":"1.8","10":"virginica"},{"1":"6.5","2":"3.0","3":"5.8","4":"2.2","5":"virginica","6":"6.5","7":"3.0","8":"5.8","9":"2.2","10":"virginica"},{"1":"7.6","2":"3.0","3":"6.6","4":"2.1","5":"virginica","6":"7.6","7":"3.0","8":"6.6","9":"2.1","10":"virginica"},{"1":"4.9","2":"2.5","3":"4.5","4":"1.7","5":"virginica","6":"4.9","7":"2.5","8":"4.5","9":"1.7","10":"virginica"},{"1":"7.3","2":"2.9","3":"6.3","4":"1.8","5":"virginica","6":"7.3","7":"2.9","8":"6.3","9":"1.8","10":"virginica"},{"1":"6.7","2":"2.5","3":"5.8","4":"1.8","5":"virginica","6":"6.7","7":"2.5","8":"5.8","9":"1.8","10":"virginica"},{"1":"7.2","2":"3.6","3":"6.1","4":"2.5","5":"virginica","6":"7.2","7":"3.6","8":"6.1","9":"2.5","10":"virginica"},{"1":"6.5","2":"3.2","3":"5.1","4":"2.0","5":"virginica","6":"6.5","7":"3.2","8":"5.1","9":"2.0","10":"virginica"},{"1":"6.4","2":"2.7","3":"5.3","4":"1.9","5":"virginica","6":"6.4","7":"2.7","8":"5.3","9":"1.9","10":"virginica"},{"1":"6.8","2":"3.0","3":"5.5","4":"2.1","5":"virginica","6":"6.8","7":"3.0","8":"5.5","9":"2.1","10":"virginica"},{"1":"5.7","2":"2.5","3":"5.0","4":"2.0","5":"virginica","6":"5.7","7":"2.5","8":"5.0","9":"2.0","10":"virginica"},{"1":"5.8","2":"2.8","3":"5.1","4":"2.4","5":"virginica","6":"5.8","7":"2.8","8":"5.1","9":"2.4","10":"virginica"},{"1":"6.4","2":"3.2","3":"5.3","4":"2.3","5":"virginica","6":"6.4","7":"3.2","8":"5.3","9":"2.3","10":"virginica"},{"1":"6.5","2":"3.0","3":"5.5","4":"1.8","5":"virginica","6":"6.5","7":"3.0","8":"5.5","9":"1.8","10":"virginica"},{"1":"7.7","2":"3.8","3":"6.7","4":"2.2","5":"virginica","6":"7.7","7":"3.8","8":"6.7","9":"2.2","10":"virginica"},{"1":"7.7","2":"2.6","3":"6.9","4":"2.3","5":"virginica","6":"7.7","7":"2.6","8":"6.9","9":"2.3","10":"virginica"},{"1":"6.0","2":"2.2","3":"5.0","4":"1.5","5":"virginica","6":"6.0","7":"2.2","8":"5.0","9":"1.5","10":"virginica"},{"1":"6.9","2":"3.2","3":"5.7","4":"2.3","5":"virginica","6":"6.9","7":"3.2","8":"5.7","9":"2.3","10":"virginica"},{"1":"5.6","2":"2.8","3":"4.9","4":"2.0","5":"virginica","6":"5.6","7":"2.8","8":"4.9","9":"2.0","10":"virginica"},{"1":"7.7","2":"2.8","3":"6.7","4":"2.0","5":"virginica","6":"7.7","7":"2.8","8":"6.7","9":"2.0","10":"virginica"},{"1":"6.3","2":"2.7","3":"4.9","4":"1.8","5":"virginica","6":"6.3","7":"2.7","8":"4.9","9":"1.8","10":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.1","5":"virginica","6":"6.7","7":"3.3","8":"5.7","9":"2.1","10":"virginica"},{"1":"7.2","2":"3.2","3":"6.0","4":"1.8","5":"virginica","6":"7.2","7":"3.2","8":"6.0","9":"1.8","10":"virginica"},{"1":"6.2","2":"2.8","3":"4.8","4":"1.8","5":"virginica","6":"6.2","7":"2.8","8":"4.8","9":"1.8","10":"virginica"},{"1":"6.1","2":"3.0","3":"4.9","4":"1.8","5":"virginica","6":"6.1","7":"3.0","8":"4.9","9":"1.8","10":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.1","5":"virginica","6":"6.4","7":"2.8","8":"5.6","9":"2.1","10":"virginica"},{"1":"7.2","2":"3.0","3":"5.8","4":"1.6","5":"virginica","6":"7.2","7":"3.0","8":"5.8","9":"1.6","10":"virginica"},{"1":"7.4","2":"2.8","3":"6.1","4":"1.9","5":"virginica","6":"7.4","7":"2.8","8":"6.1","9":"1.9","10":"virginica"},{"1":"7.9","2":"3.8","3":"6.4","4":"2.0","5":"virginica","6":"7.9","7":"3.8","8":"6.4","9":"2.0","10":"virginica"},{"1":"6.4","2":"2.8","3":"5.6","4":"2.2","5":"virginica","6":"6.4","7":"2.8","8":"5.6","9":"2.2","10":"virginica"},{"1":"6.3","2":"2.8","3":"5.1","4":"1.5","5":"virginica","6":"6.3","7":"2.8","8":"5.1","9":"1.5","10":"virginica"},{"1":"6.1","2":"2.6","3":"5.6","4":"1.4","5":"virginica","6":"6.1","7":"2.6","8":"5.6","9":"1.4","10":"virginica"},{"1":"7.7","2":"3.0","3":"6.1","4":"2.3","5":"virginica","6":"7.7","7":"3.0","8":"6.1","9":"2.3","10":"virginica"},{"1":"6.3","2":"3.4","3":"5.6","4":"2.4","5":"virginica","6":"6.3","7":"3.4","8":"5.6","9":"2.4","10":"virginica"},{"1":"6.4","2":"3.1","3":"5.5","4":"1.8","5":"virginica","6":"6.4","7":"3.1","8":"5.5","9":"1.8","10":"virginica"},{"1":"6.0","2":"3.0","3":"4.8","4":"1.8","5":"virginica","6":"6.0","7":"3.0","8":"4.8","9":"1.8","10":"virginica"},{"1":"6.9","2":"3.1","3":"5.4","4":"2.1","5":"virginica","6":"6.9","7":"3.1","8":"5.4","9":"2.1","10":"virginica"},{"1":"6.7","2":"3.1","3":"5.6","4":"2.4","5":"virginica","6":"6.7","7":"3.1","8":"5.6","9":"2.4","10":"virginica"},{"1":"6.9","2":"3.1","3":"5.1","4":"2.3","5":"virginica","6":"6.9","7":"3.1","8":"5.1","9":"2.3","10":"virginica"},{"1":"5.8","2":"2.7","3":"5.1","4":"1.9","5":"virginica","6":"5.8","7":"2.7","8":"5.1","9":"1.9","10":"virginica"},{"1":"6.8","2":"3.2","3":"5.9","4":"2.3","5":"virginica","6":"6.8","7":"3.2","8":"5.9","9":"2.3","10":"virginica"},{"1":"6.7","2":"3.3","3":"5.7","4":"2.5","5":"virginica","6":"6.7","7":"3.3","8":"5.7","9":"2.5","10":"virginica"},{"1":"6.7","2":"3.0","3":"5.2","4":"2.3","5":"virginica","6":"6.7","7":"3.0","8":"5.2","9":"2.3","10":"virginica"},{"1":"6.3","2":"2.5","3":"5.0","4":"1.9","5":"virginica","6":"6.3","7":"2.5","8":"5.0","9":"1.9","10":"virginica"},{"1":"6.5","2":"3.0","3":"5.2","4":"2.0","5":"virginica","6":"6.5","7":"3.0","8":"5.2","9":"2.0","10":"virginica"},{"1":"6.2","2":"3.4","3":"5.4","4":"2.3","5":"virginica","6":"6.2","7":"3.4","8":"5.4","9":"2.3","10":"virginica"},{"1":"5.9","2":"3.0","3":"5.1","4":"1.8","5":"virginica","6":"5.9","7":"3.0","8":"5.1","9":"1.8","10":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Summarize data

Summary calculation on single column


```r
summarize(iris_df, mean(Petal.Length))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["mean(Petal.Length)"],"name":[1],"type":["dbl"],"align":["right"]}],"data":[{"1":"3.758"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Summary calculation on many columns


```r
summarize_all(iris_df[,1:4], mean)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"5.843333","2":"3.057333","3":"3.758","4":"1.199333"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Summarize by grouping column


```r
summarize(group_by(iris_df, Species), mean(Petal.Length))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Species"],"name":[1],"type":["chr"],"align":["left"]},{"label":["mean(Petal.Length)"],"name":[2],"type":["dbl"],"align":["right"]}],"data":[{"1":"setosa","2":"1.462"},{"1":"versicolor","2":"4.260"},{"1":"virginica","2":"5.552"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Aggregate summaries


```r
summarize_all(group_by(iris_df, Species), mean) 
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Species"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Sepal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"setosa","2":"5.006","3":"3.428","4":"1.462","5":"0.246"},{"1":"versicolor","2":"5.936","3":"2.770","4":"4.260","5":"1.326"},{"1":"virginica","2":"6.588","3":"2.974","4":"5.552","5":"2.026"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Note: `group_by` does the looping for the user similar to `aggregate` or `tapply`.


### Merging tibbles

The `dplyr` package provides several join functions for merging `tibbles` by a common key column
similar to the `merge` function in base R. These `*_join` functions include: 

* `inner_join()`: returns join only for rows matching among both `tibbles`
* `full_join()`: returns join for all (matching and non-matching) rows of two `tibbles` 
* `left_join()`: returns join for all rows in first `tibble` 
* `right_join()`: returns join for all rows in second `tibble`
* `anti_join()`: returns for first `tibble` only those rows that have no match in the second one

Sample `tibbles` to illustrate `*.join` functions.


```r
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g1","2":"1","3":"11","4":"21","5":"31"},{"1":"g2","2":"2","3":"12","4":"22","5":"32"},{"1":"g3","2":"3","3":"13","4":"23","5":"33"},{"1":"g4","2":"4","3":"14","4":"24","5":"34"},{"1":"g5","2":"5","3":"15","4":"25","5":"35"},{"1":"g6","2":"6","3":"16","4":"26","5":"36"},{"1":"g7","2":"7","3":"17","4":"27","5":"37"},{"1":"g8","2":"8","3":"18","4":"28","5":"38"},{"1":"g9","2":"9","3":"19","4":"29","5":"39"},{"1":"g10","2":"10","3":"20","4":"30","5":"40"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
df2 <- bind_cols(data_frame(ids2=paste0("g", c(2,5,11,12))), as_data_frame(matrix(1:16, 4, 4, dimnames=list(1:4, paste0("CB", 1:4)))))
df2
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids2"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CB1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CB2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CB3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CB4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g2","2":"1","3":"5","4":"9","5":"13"},{"1":"g5","2":"2","3":"6","4":"10","5":"14"},{"1":"g11","2":"3","3":"7","4":"11","5":"15"},{"1":"g12","2":"4","3":"8","4":"12","5":"16"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
#### Inner join


```r
inner_join(df1, df2, by=c("ids1"="ids2"))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]},{"label":["CB1"],"name":[6],"type":["int"],"align":["right"]},{"label":["CB2"],"name":[7],"type":["int"],"align":["right"]},{"label":["CB3"],"name":[8],"type":["int"],"align":["right"]},{"label":["CB4"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"g2","2":"2","3":"12","4":"22","5":"32","6":"1","7":"5","8":"9","9":"13"},{"1":"g5","2":"5","3":"15","4":"25","5":"35","6":"2","7":"6","8":"10","9":"14"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Left join


```r
left_join(df1, df2, by=c("ids1"="ids2"))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]},{"label":["CB1"],"name":[6],"type":["int"],"align":["right"]},{"label":["CB2"],"name":[7],"type":["int"],"align":["right"]},{"label":["CB3"],"name":[8],"type":["int"],"align":["right"]},{"label":["CB4"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"g1","2":"1","3":"11","4":"21","5":"31","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g2","2":"2","3":"12","4":"22","5":"32","6":"1","7":"5","8":"9","9":"13"},{"1":"g3","2":"3","3":"13","4":"23","5":"33","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g4","2":"4","3":"14","4":"24","5":"34","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g5","2":"5","3":"15","4":"25","5":"35","6":"2","7":"6","8":"10","9":"14"},{"1":"g6","2":"6","3":"16","4":"26","5":"36","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g7","2":"7","3":"17","4":"27","5":"37","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g8","2":"8","3":"18","4":"28","5":"38","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g9","2":"9","3":"19","4":"29","5":"39","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g10","2":"10","3":"20","4":"30","5":"40","6":"NA","7":"NA","8":"NA","9":"NA"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Right join


```r
right_join(df1, df2, by=c("ids1"="ids2"))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]},{"label":["CB1"],"name":[6],"type":["int"],"align":["right"]},{"label":["CB2"],"name":[7],"type":["int"],"align":["right"]},{"label":["CB3"],"name":[8],"type":["int"],"align":["right"]},{"label":["CB4"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"g2","2":"2","3":"12","4":"22","5":"32","6":"1","7":"5","8":"9","9":"13"},{"1":"g5","2":"5","3":"15","4":"25","5":"35","6":"2","7":"6","8":"10","9":"14"},{"1":"g11","2":"NA","3":"NA","4":"NA","5":"NA","6":"3","7":"7","8":"11","9":"15"},{"1":"g12","2":"NA","3":"NA","4":"NA","5":"NA","6":"4","7":"8","8":"12","9":"16"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Full join


```r
full_join(df1, df2, by=c("ids1"="ids2"))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]},{"label":["CB1"],"name":[6],"type":["int"],"align":["right"]},{"label":["CB2"],"name":[7],"type":["int"],"align":["right"]},{"label":["CB3"],"name":[8],"type":["int"],"align":["right"]},{"label":["CB4"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"g1","2":"1","3":"11","4":"21","5":"31","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g2","2":"2","3":"12","4":"22","5":"32","6":"1","7":"5","8":"9","9":"13"},{"1":"g3","2":"3","3":"13","4":"23","5":"33","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g4","2":"4","3":"14","4":"24","5":"34","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g5","2":"5","3":"15","4":"25","5":"35","6":"2","7":"6","8":"10","9":"14"},{"1":"g6","2":"6","3":"16","4":"26","5":"36","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g7","2":"7","3":"17","4":"27","5":"37","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g8","2":"8","3":"18","4":"28","5":"38","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g9","2":"9","3":"19","4":"29","5":"39","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g10","2":"10","3":"20","4":"30","5":"40","6":"NA","7":"NA","8":"NA","9":"NA"},{"1":"g11","2":"NA","3":"NA","4":"NA","5":"NA","6":"3","7":"7","8":"11","9":"15"},{"1":"g12","2":"NA","3":"NA","4":"NA","5":"NA","6":"4","7":"8","8":"12","9":"16"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

#### Anti join


```r
anti_join(df1, df2, by=c("ids1"="ids2"))
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids1"],"name":[1],"type":["chr"],"align":["left"]},{"label":["CA1"],"name":[2],"type":["int"],"align":["right"]},{"label":["CA2"],"name":[3],"type":["int"],"align":["right"]},{"label":["CA3"],"name":[4],"type":["int"],"align":["right"]},{"label":["CA4"],"name":[5],"type":["int"],"align":["right"]}],"data":[{"1":"g1","2":"1","3":"11","4":"21","5":"31"},{"1":"g3","2":"3","3":"13","4":"23","5":"33"},{"1":"g4","2":"4","3":"14","4":"24","5":"34"},{"1":"g6","2":"6","3":"16","4":"26","5":"36"},{"1":"g7","2":"7","3":"17","4":"27","5":"37"},{"1":"g8","2":"8","3":"18","4":"28","5":"38"},{"1":"g9","2":"9","3":"19","4":"29","5":"39"},{"1":"g10","2":"10","3":"20","4":"30","5":"40"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

For additional join options users want to cosult the `*_join` help pages.

# Outline

- Overview
- Install
- File Import and Export
- Usage
- <div class="white">__Chaining (Pipes)__</div>
- SQLite Databases
- References

## Chaining

<br/><br/>
<br/><br/>

To simplify chaining of serveral operations, `dplyr` provides the `%>%`
operator, where `x %>% f(y)` turns into `f(x, y)`. This way one can pipe
together multiple operations by writing them from left-to-right or
top-to-bottom. This makes for easy to type and readable code.

<br/><br/>
<br/><br/>
<br/><br/>
<br/><br/>

<p style='text-align: right;'> __[ scroll down to continue ]__ </p>
<br/><br/>
<br/><br/>


### Example 1

Series of data manipulations and export


```r
read_tsv("iris.txt") %>% # Import with read_tbv from readr package
    as_tibble() %>% # Declare to use tibble
    select(Sepal.Length:Species) %>% # Select columns
    filter(Species=="setosa") %>% # Filter rows by some value
    arrange(Sepal.Length) %>% # Sort by some column
    mutate(Subtract=Petal.Length - Petal.Width) # Calculate and append
```

```
## 
## ── Column specification ────────────────────────────────────────────────────────
## cols(
##   Sepal.Length = col_double(),
##   Sepal.Width = col_double(),
##   Petal.Length = col_double(),
##   Petal.Width = col_double(),
##   Species = col_character()
## )
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sepal.Length"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[5],"type":["chr"],"align":["left"]},{"label":["Subtract"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"4.3","2":"3.0","3":"1.1","4":"0.1","5":"setosa","6":"1.0"},{"1":"4.4","2":"2.9","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"4.4","2":"3.0","3":"1.3","4":"0.2","5":"setosa","6":"1.1"},{"1":"4.4","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"1.1"},{"1":"4.5","2":"2.3","3":"1.3","4":"0.3","5":"setosa","6":"1.0"},{"1":"4.6","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"4.6","2":"3.4","3":"1.4","4":"0.3","5":"setosa","6":"1.1"},{"1":"4.6","2":"3.6","3":"1.0","4":"0.2","5":"setosa","6":"0.8"},{"1":"4.6","2":"3.2","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"4.7","2":"3.2","3":"1.3","4":"0.2","5":"setosa","6":"1.1"},{"1":"4.7","2":"3.2","3":"1.6","4":"0.2","5":"setosa","6":"1.4"},{"1":"4.8","2":"3.4","3":"1.6","4":"0.2","5":"setosa","6":"1.4"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.1","5":"setosa","6":"1.3"},{"1":"4.8","2":"3.4","3":"1.9","4":"0.2","5":"setosa","6":"1.7"},{"1":"4.8","2":"3.1","3":"1.6","4":"0.2","5":"setosa","6":"1.4"},{"1":"4.8","2":"3.0","3":"1.4","4":"0.3","5":"setosa","6":"1.1"},{"1":"4.9","2":"3.0","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.1","5":"setosa","6":"1.4"},{"1":"4.9","2":"3.1","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"4.9","2":"3.6","3":"1.4","4":"0.1","5":"setosa","6":"1.3"},{"1":"5.0","2":"3.6","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"5.0","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"5.0","2":"3.0","3":"1.6","4":"0.2","5":"setosa","6":"1.4"},{"1":"5.0","2":"3.4","3":"1.6","4":"0.4","5":"setosa","6":"1.2"},{"1":"5.0","2":"3.2","3":"1.2","4":"0.2","5":"setosa","6":"1.0"},{"1":"5.0","2":"3.5","3":"1.3","4":"0.3","5":"setosa","6":"1.0"},{"1":"5.0","2":"3.5","3":"1.6","4":"0.6","5":"setosa","6":"1.0"},{"1":"5.0","2":"3.3","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"5.1","2":"3.5","3":"1.4","4":"0.3","5":"setosa","6":"1.1"},{"1":"5.1","2":"3.8","3":"1.5","4":"0.3","5":"setosa","6":"1.2"},{"1":"5.1","2":"3.7","3":"1.5","4":"0.4","5":"setosa","6":"1.1"},{"1":"5.1","2":"3.3","3":"1.7","4":"0.5","5":"setosa","6":"1.2"},{"1":"5.1","2":"3.4","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"5.1","2":"3.8","3":"1.9","4":"0.4","5":"setosa","6":"1.5"},{"1":"5.1","2":"3.8","3":"1.6","4":"0.2","5":"setosa","6":"1.4"},{"1":"5.2","2":"3.5","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"5.2","2":"3.4","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"5.2","2":"4.1","3":"1.5","4":"0.1","5":"setosa","6":"1.4"},{"1":"5.3","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"5.4","2":"3.9","3":"1.7","4":"0.4","5":"setosa","6":"1.3"},{"1":"5.4","2":"3.7","3":"1.5","4":"0.2","5":"setosa","6":"1.3"},{"1":"5.4","2":"3.9","3":"1.3","4":"0.4","5":"setosa","6":"0.9"},{"1":"5.4","2":"3.4","3":"1.7","4":"0.2","5":"setosa","6":"1.5"},{"1":"5.4","2":"3.4","3":"1.5","4":"0.4","5":"setosa","6":"1.1"},{"1":"5.5","2":"4.2","3":"1.4","4":"0.2","5":"setosa","6":"1.2"},{"1":"5.5","2":"3.5","3":"1.3","4":"0.2","5":"setosa","6":"1.1"},{"1":"5.7","2":"4.4","3":"1.5","4":"0.4","5":"setosa","6":"1.1"},{"1":"5.7","2":"3.8","3":"1.7","4":"0.3","5":"setosa","6":"1.4"},{"1":"5.8","2":"4.0","3":"1.2","4":"0.2","5":"setosa","6":"1.0"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
    # write_tsv("iris.txt") # Export to file, omitted here to show result 
```

### Example 2

Series of summary calculations for grouped data (`group_by`)


```r
iris_df %>% # Declare tibble to use 
    group_by(Species) %>% # Group by species
    summarize(Mean_Sepal.Length=mean(Sepal.Length), 
              Max_Sepal.Length=max(Sepal.Length),
              Min_Sepal.Length=min(Sepal.Length),
              SD_Sepal.Length=sd(Sepal.Length),
              Total=n()) 
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["Species"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Mean_Sepal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Max_Sepal.Length"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Min_Sepal.Length"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["SD_Sepal.Length"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Total"],"name":[6],"type":["int"],"align":["right"]}],"data":[{"1":"setosa","2":"5.006","3":"5.8","4":"4.3","5":"0.3524897","6":"50"},{"1":"versicolor","2":"5.936","3":"7.0","4":"4.9","5":"0.5161711","6":"50"},{"1":"virginica","2":"6.588","3":"7.9","4":"4.9","5":"0.6358796","6":"50"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Example 3

Combining `dplyr` chaining with `ggplot`


```r
iris_df %>% 
    group_by(Species) %>% 
    summarize_all(mean) %>% 
    reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%
    ggplot(aes(Samples, Values, fill = Species)) + 
           geom_bar(position="dodge", stat="identity")
```

![](dplyr_slides_files/figure-html/plyr_chaining3-1.png)<!-- -->

# Outline

- Overview
- Install
- File Import and Export
- Usage
- Chaining (Pipes)
- <div class="white">__SQLite Databases__</div>
- References

## SQLite Databases

<br/><br/>
<br/><br/>

`SQLite` is a lightweight relational database solution. The `RSQLite` package provides an easy to use interface to create, manage and query `SQLite` databases directly from R. Basic instructions
for using `SQLite` from the command-line are available [here](https://www.sqlite.org/cli.html). A short introduction to `RSQLite` is available [here](https://github.com/rstats-db/RSQLite/blob/master/vignettes/RSQLite.Rmd).

## Loading data into SQLite databases

The following loads two `data.frames` derived from the `iris` data set (here `mydf1` and `mydf2`) 
into an SQLite database (here `test.db`).


```r
library(RSQLite)
unlink("test.db") # Delete any existing test.db
mydb <- dbConnect(SQLite(), "test.db") # Creates database file test.db
mydf1 <- data.frame(ids=paste0("id", seq_along(iris[,1])), iris)
mydf2 <- mydf1[sample(seq_along(mydf1[,1]), 10),]
dbWriteTable(mydb, "mydf1", mydf1)
dbWriteTable(mydb, "mydf2", mydf2)
```


### List names of tables in database


```r
dbListTables(mydb)
```

```
## [1] "mydf1" "mydf2"
```

<p style='text-align: right;'> __[ scroll down to continue ]__ </p>
<br/><br/>
<br/><br/>


### Import table into `data.frame`


```r
dbGetQuery(mydb, 'SELECT * FROM mydf2')
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Sepal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"id75","2":"6.4","3":"2.9","4":"4.3","5":"1.3","6":"versicolor"},{"1":"id5","2":"5.0","3":"3.6","4":"1.4","5":"0.2","6":"setosa"},{"1":"id49","2":"5.3","3":"3.7","4":"1.5","5":"0.2","6":"setosa"},{"1":"id83","2":"5.8","3":"2.7","4":"3.9","5":"1.2","6":"versicolor"},{"1":"id48","2":"4.6","3":"3.2","4":"1.4","5":"0.2","6":"setosa"},{"1":"id133","2":"6.4","3":"2.8","4":"5.6","5":"2.2","6":"virginica"},{"1":"id43","2":"4.4","3":"3.2","4":"1.3","5":"0.2","6":"setosa"},{"1":"id31","2":"4.8","3":"3.1","4":"1.6","5":"0.2","6":"setosa"},{"1":"id63","2":"6.0","3":"2.2","4":"4.0","5":"1.0","6":"versicolor"},{"1":"id89","2":"5.6","3":"3.0","4":"4.1","5":"1.3","6":"versicolor"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Query database


```r
dbGetQuery(mydb, 'SELECT * FROM mydf1 WHERE "Sepal.Length" < 4.6')
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Sepal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[6],"type":["chr"],"align":["left"]}],"data":[{"1":"id9","2":"4.4","3":"2.9","4":"1.4","5":"0.2","6":"setosa"},{"1":"id14","2":"4.3","3":"3.0","4":"1.1","5":"0.1","6":"setosa"},{"1":"id39","2":"4.4","3":"3.0","4":"1.3","5":"0.2","6":"setosa"},{"1":"id42","2":"4.5","3":"2.3","4":"1.3","5":"0.3","6":"setosa"},{"1":"id43","2":"4.4","3":"3.2","4":"1.3","5":"0.2","6":"setosa"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

### Join tables

The two tables can be joined on the shared `ids` column as follows. 


```r
dbGetQuery(mydb, 'SELECT * FROM mydf1, mydf2 WHERE mydf1.ids = mydf2.ids')
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["ids"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Sepal.Length"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[6],"type":["chr"],"align":["left"]},{"label":["ids"],"name":[7],"type":["chr"],"align":["left"]},{"label":["Sepal.Length"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["Sepal.Width"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["Petal.Length"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["Petal.Width"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["Species"],"name":[12],"type":["chr"],"align":["left"]}],"data":[{"1":"id5","2":"5.0","3":"3.6","4":"1.4","5":"0.2","6":"setosa","7":"id5","8":"5.0","9":"3.6","10":"1.4","11":"0.2","12":"setosa"},{"1":"id31","2":"4.8","3":"3.1","4":"1.6","5":"0.2","6":"setosa","7":"id31","8":"4.8","9":"3.1","10":"1.6","11":"0.2","12":"setosa"},{"1":"id43","2":"4.4","3":"3.2","4":"1.3","5":"0.2","6":"setosa","7":"id43","8":"4.4","9":"3.2","10":"1.3","11":"0.2","12":"setosa"},{"1":"id48","2":"4.6","3":"3.2","4":"1.4","5":"0.2","6":"setosa","7":"id48","8":"4.6","9":"3.2","10":"1.4","11":"0.2","12":"setosa"},{"1":"id49","2":"5.3","3":"3.7","4":"1.5","5":"0.2","6":"setosa","7":"id49","8":"5.3","9":"3.7","10":"1.5","11":"0.2","12":"setosa"},{"1":"id63","2":"6.0","3":"2.2","4":"4.0","5":"1.0","6":"versicolor","7":"id63","8":"6.0","9":"2.2","10":"4.0","11":"1.0","12":"versicolor"},{"1":"id75","2":"6.4","3":"2.9","4":"4.3","5":"1.3","6":"versicolor","7":"id75","8":"6.4","9":"2.9","10":"4.3","11":"1.3","12":"versicolor"},{"1":"id83","2":"5.8","3":"2.7","4":"3.9","5":"1.2","6":"versicolor","7":"id83","8":"5.8","9":"2.7","10":"3.9","11":"1.2","12":"versicolor"},{"1":"id89","2":"5.6","3":"3.0","4":"4.1","5":"1.3","6":"versicolor","7":"id89","8":"5.6","9":"3.0","10":"4.1","11":"1.3","12":"versicolor"},{"1":"id133","2":"6.4","3":"2.8","4":"5.6","5":"2.2","6":"virginica","7":"id133","8":"6.4","9":"2.8","10":"5.6","11":"2.2","12":"virginica"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

## Session Info


```r
sessionInfo()
```

```
## R version 4.1.0 (2021-05-18)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 10 (buster)
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] RSQLite_2.2.7     data.table_1.14.0 forcats_0.5.1     stringr_1.4.0    
##  [5] dplyr_1.0.6       purrr_0.3.4       readr_1.4.0       tidyr_1.1.3      
##  [9] tibble_3.1.2      ggplot2_3.3.3     tidyverse_1.3.1  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.6        lubridate_1.7.10  ps_1.6.0          assertthat_0.2.1 
##  [5] digest_0.6.27     utf8_1.2.1        R6_2.5.0          cellranger_1.1.0 
##  [9] plyr_1.8.6        backports_1.2.1   reprex_2.0.0      evaluate_0.14    
## [13] httr_1.4.2        highr_0.9         pillar_1.6.1      rlang_0.4.11     
## [17] readxl_1.3.1      rstudioapi_0.13   blob_1.2.1        jquerylib_0.1.4  
## [21] rmarkdown_2.8     labeling_0.4.2    bit_4.0.4         munsell_0.5.0    
## [25] broom_0.7.6       compiler_4.1.0    modelr_0.1.8      xfun_0.23        
## [29] pkgconfig_2.0.3   htmltools_0.5.1.1 tidyselect_1.1.1  fansi_0.4.2      
## [33] crayon_1.4.1      dbplyr_2.1.1      withr_2.4.2       grid_4.1.0       
## [37] jsonlite_1.7.2    gtable_0.3.0      lifecycle_1.0.0   DBI_1.1.1        
## [41] magrittr_2.0.1    scales_1.1.1      cachem_1.0.5      cli_2.5.0        
## [45] stringi_1.6.2     farver_2.1.0      reshape2_1.4.4    fs_1.5.0         
## [49] xml2_1.3.2        bslib_0.2.5.1     ellipsis_0.3.2    generics_0.1.0   
## [53] vctrs_0.3.8       tools_4.1.0       bit64_4.0.5       glue_1.4.2       
## [57] hms_1.1.0         fastmap_1.1.0     yaml_2.2.1        colorspace_2.0-1 
## [61] rvest_1.0.0       memoise_2.0.0     knitr_1.33        haven_2.4.1      
## [65] sass_0.4.0
```

# Outline

- Overview
- Install
- File Import and Export
- Usage
- Chaining (Pipes)
- SQLite Databases
- <div class="white">__References__</div>

## References



