---
title: "Programming in R" 
author: "Author: Thomas Girke"
date: "Last update: 06 June, 2021" 
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

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('Programming_in_R.Rmd', c('html_document'), clean=F); knitr::knit('Programming_in_R.Rmd', tangle=TRUE)"; Rscript ../md2jekyll.R Programming_in_R.knit.md 9; Rscript -e "rmarkdown::render('Programming_in_R.Rmd', c('pdf_document'))"
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
[ [.Rmd](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rprogramming/rprogramming.Rmd) ] &nbsp; &nbsp; 
[ [.R](https://raw.githubusercontent.com/tgirke/GEN242//main/content/en/tutorials/rprogramming/rprogramming.R) ]
</div>


## Overview

One of the main attractions of using the R
([http://cran.at.r-project.org](http://cran.at.r-project.org)) environment is
the ease with which users can write their own programs and custom functions.
The R programming syntax is extremely easy to learn, even for users with no
previous programming experience. Once the basic R programming control
structures are understood, users can use the R language as a powerful
environment to perform complex custom analyses of almost any type of data [@Gentleman2008-xo].


## Why Programming in R?

* Powerful statistical environment and programming language
* Facilitates reproducible research
* Efficient data structures make programming very easy
* Ease of implementing custom functions
* Powerful graphics
* Access to fast growing number of analysis packages
* One of the most widely used languages in bioinformatics
* Is standard for data mining and biostatistical analysis
* Technical advantages: free, open-source, available for all OSs


## R Basics 

The previous [Rbasics](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/rbasics/rbasics/) tutorial provides a general introduction to the usage of the R environment and its basic command syntax.
More details can be found in the R & BioConductor manual [here](http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual).

## Code Editors for R

Several excellent code editors are available that provide functionalities like R syntax highlighting, auto code indenting and utilities to send code/functions to the R console.

* [RStudio](https://www.rstudio.com/products/rstudio/features/): GUI-based IDE for R
* [Vim-R-Tmux](https://girke.bioinformatics.ucr.edu/GEN242/tutorials/linux/linux/#nvim-r-tmux-essentials): R working environment based on vim and tmux
* [Emacs](http://www.xemacs.org/Download/index.html) ([ESS add-on package](http://ess.r-project.org/))
* [gedit](https://wiki.gnome.org/Apps/Gedit) and [Rgedit](https://wiki.gnome.org/Apps/Gedit)
* [RKWard](https://rkward.kde.org/)
* [Eclipse](http://www.walware.de/goto/statet)
* [Tinn-R](https://sourceforge.net/projects/tinn-r/files/Tinn-R%20portable/)
* [Notepad++ (NppToR)](https://sourceforge.net/projects/npptor/)

<center>Programming in R using RStudio</center>
<center><img title="R_Interfaces" src="../images/rstudio.png"/></center>

<center> Programming in R using Vim or Emacs</center>
<center><img title="vim-r" src="../images/vimR.png"/></center>

## Finding Help

Reference list on R programming (selection)

* [Advanced R](http://adv-r.had.co.nz/), by Hadley Wickham
* [R Programming for Bioinformatics](http://master.bioconductor.org/help/publications/books/r-programming-for-bioinformatics/), by Robert Gentleman
* [S Programming](http://www.stats.ox.ac.uk/pub/MASS3/Sprog/), by W. N. Venables and B. D. Ripley
* [Programming with Data](http://www.amazon.com/Programming-Data-Language-Lecture-Economics/dp/0387985034), by John M. Chambers
* [R Help](http://www1.maths.lth.se/help/R/) & [R Coding Conventions](http://www1.maths.lth.se/help/R/RCC/), Henrik Bengtsson, Lund University
* [Programming in R](http://zoonek2.free.fr/UNIX/48_R/02.html) (Vincent Zoonekynd)
* [Peter's R Programming Pages](http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r), University of Warwick
* [Rtips](http://pj.freefaculty.org/R/statsRus.html), Paul Johnsson, University of Kansas
* [R for Programmers](http://heather.cs.ucdavis.edu/~matloff/r.html), Norm Matloff, UC Davis
* [High-Performance R](http://www.statistik.uni-dortmund.de/useR-2008/tutorials/useR2008introhighperfR.pdf), Dirk Eddelbuettel tutorial presented at useR-2008
* [C/C++ level programming for R](http://www.stat.harvard.edu/ccr2005/index.html), Gopi Goswami


## Control Structures

### Important Operators

#### Comparison operators

* `==` (equal)
* `!=` (not equal)
* `>` (greater than)
* `>=` (greater than or equal)
* `<` (less than)
* `<=` (less than or equal)

#### Logical operators
		
* `&` (and) 
* `&&` (and) 
* `|` (or) 
* `||` (or) 
* `!` (not)

Note: `&` and `&&` indicate logical AND, while `|` and `||` indicate logical OR. The shorter form performs element-wise comparisons of same-length vectors. 
The longer form evaluates left to right examining only the first element of each vector (can be of different lengths). Evaluation proceeds only until the result 
is determined. The longer form is preferred for programming control-flow, _e.g._ via `if` clauses.

### Conditional Executions: `if` Statements

An `if` statement operates on length-one logical vectors.

__Syntax__

```r
if (TRUE) { 
    statements_1 
} else { 
    statements_2 
}
```

In the `else` component, avoid inserting newlines between `} else`. For details on how to best and consistently format R code, 
this [style guide](http://adv-r.had.co.nz/Style.html) is a good start. In addition, the [`formatR`](https://yihui.org/formatr/) package can be helpful.

__Example__

```r
if (1==0) { 
    print(1) 
} else { 
    print(2) 
}
```

```
## [1] 2
```

### Conditional Executions: `ifelse` Statements

The `ifelse` statement operates on vectors.

__Syntax__

```r
ifelse(test, true_value, false_value)
```
__Example__

```r
x <- 1:10 
ifelse(x<5, sqrt(x), 0)
```

```
##  [1] 1.000000 1.414214 1.732051 2.000000 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000
```

## Loops

### `for` loop

`for` loops iterate over elements of a looping vector.

__Syntax__

```r
for(variable in sequence) { 
	statements 
}
```
__Example__

```r
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

__Note:__ Inject into objecs is much faster than append approach with `c`, `cbind`, etc.

__Example__

```r
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]
```

```
## [1] 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

#### Conditional Stop of Loops

The `stop` function can be used to break out of a loop (or a function) when a condition becomes `TRUE`. In addition, an error message will be printed.

__Example__

```r
x <- 1:10
z <- NULL
for(i in seq(along=x)) { 
	if (x[i] < 5) { 
		z <- c(z, x[i]-1)  
	} else { 
		stop("values need to be < 5") 
	}
}
```

### `while` loop

Iterates as long as a condition is true.

__Syntax__

```r
while(condition) {
	statements
}
```

__Example__

```r
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}
```

```
## [1] 2
## [1] 4
## [1] 6
```

### The `apply` Function Family

#### `apply`

__Syntax__

```r
apply(X, MARGIN, FUN, ARGs)
```

__Arguments__

* `X`: `array`, `matrix` or `data.frame`
* `MARGIN`: `1` for rows, `2` for columns
* `FUN`: one or more functions
* `ARGs`: possible arguments for functions

__Example__

```r
apply(iris[1:8,1:3], 1, mean)
```

```
##        1        2        3        4        5        6        7        8 
## 3.333333 3.100000 3.066667 3.066667 3.333333 3.666667 3.133333 3.300000
```

#### `tapply`

Applies a function to vector components that are defined by a factor.

__Syntax__

```r
tapply(vector, factor, FUN)
```

__Example__

```r
iris[1:2,]
```

```
##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
## 1          5.1         3.5          1.4         0.2  setosa
## 2          4.9         3.0          1.4         0.2  setosa
```

```r
tapply(iris$Sepal.Length, iris$Species, mean)
```

```
##     setosa versicolor  virginica 
##      5.006      5.936      6.588
```

#### `sapply`, `lapply` and `vapply`

The iterator functions `sapply`, `lapply` and `vapply` apply a function to
vectors or lists. The `lapply` function always returns a list, while `sapply`
returns `vector` or `matrix` objects when possible. If not then a list is
returned.  The `vapply` function returns a vector or array of type matching the
`FUN.VALUE`. Compared to `sapply`, `vapply` is a safer choice with respect to
controlling specific output types to avoid exception handling problems.

__Examples__

```r
l <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(l, mean)
```

```
## $a
## [1] 5.5
## 
## $beta
## [1] 4.535125
## 
## $logic
## [1] 0.5
```

```r
sapply(l, mean)
```

```
##        a     beta    logic 
## 5.500000 4.535125 0.500000
```

```r
vapply(l, mean, FUN.VALUE=numeric(1))
```

```
##        a     beta    logic 
## 5.500000 4.535125 0.500000
```

Often used in combination with a function definition:

```r
lapply(names(l), function(x) mean(l[[x]]))
sapply(names(l), function(x) mean(l[[x]]))
vapply(names(l), function(x) mean(l[[x]]), FUN.VALUE=numeric(1))
```

### Improving Speed Performance of Loops

Looping over very large data sets can become slow in R. However, this
limitation can be overcome by eliminating certain operations in loops or
avoiding loops over the data intensive dimension in an object altogether. The
latter can be achieved by performing mainly vector-to-vecor or
matrix-to-matrix computations. These vectorized operations run in R often over
100 times faster than the corresponding `for()` or `apply()` loops. In
addition, one can make use of the existing speed-optimized C-level functions in
R, such as `rowSums`, `rowMeans`, `table`, and `tabulate`.  Moreover, one can
design custom functions that avoid expensive R loops by using vector- or
matrix-based approaches. Alternatively, one can write programs that will
perform all time consuming computations on the C-level.

The following code samples illustrate the time-performance differences among 
the different approaches of running iterative operations in R.

#### 1. `for` loop with append versus inject approach

The following runs a `for` loop where the result is appended in each iteration
with the `c()` function. The corresponding `cbind` and `rbind` for two dimensional 
data objects would have a similar performance impact as `c()`.


```r
myMA <- matrix(rnorm(1000000), 100000, 10, dimnames=list(1:100000, paste("C", 1:10, sep="")))
results <- NULL
system.time(for(i in seq(along=myMA[,1])) results <- c(results, mean(myMA[i,])))
   user  system elapsed
 39.156   6.369  45.559
```

Now the for loop is run with an inject approach for storing the results in each iteration.


```r
results <- numeric(length(myMA[,1]))
system.time(for(i in seq(along=myMA[,1])) results[i] <- mean(myMA[i,]))
   user  system elapsed
  1.550   0.005   1.556 
```

As one can see from the output of `system.time`, the inject approach is 20-50 times faster. 

#### 2. `apply` loop versus `rowMeans` 

The following performs a row-wise mean calculation on a large matrix first with an `apply` 
loop and then with the `rowMeans` function.


```r
system.time(myMAmean <- apply(myMA, 1, mean))
  user  system elapsed
 1.452   0.005   1.456

system.time(myMAmean <- rowMeans(myMA))
   user  system elapsed
  0.005   0.001   0.006
```

Based on the results from `system.time`, the `rowMeans` approach is over 200 times faster 
than the `apply` loop.

#### 3. `apply` loop versus vectorized approach

In this example row-wise standard deviations are computed with an `apply` loop and then 
in a vectorized manner.


```r
system.time(myMAsd <- apply(myMA, 1, sd))
   user  system elapsed
  3.707   0.014   3.721
myMAsd[1:4]
        1         2         3         4
0.8505795 1.3419460 1.3768646 1.3005428

system.time(myMAsd <- sqrt((rowSums((myMA-rowMeans(myMA))^2)) / (length(myMA[1,])-1)))
   user  system elapsed
  0.020   0.009   0.028

myMAsd[1:4]
        1         2         3         4
0.8505795 1.3419460 1.3768646 1.3005428
```

The vector-based approach in the last step is over 200 times faster than the apply loop.


## Functions

### Function Overview

A very useful feature of the R environment is the possibility to expand existing functions and to easily write custom functions. In fact, most of the R software can be viewed as a series of R functions.

__Syntax__ to define function

```r
myfct <- function(arg1, arg2, ...) { 
	function_body 
}
```
__Syntax__ to call functions

```r
myfct(arg1=..., arg2=...)
```
The value returned by a function is the value of the function body, which is usually an unassigned final expression, _e.g._: `return()`

### Function Syntax Rules
	
__General__

* Functions are defined by 
    1. The assignment with the keyword `function`
    2. The declaration of arguments/variables (`arg1, arg2, ...`) 
    3. The definition of operations (`function_body`) that perform computations on the provided arguments. A function name needs to be assigned to call the function.

__Naming__ 

* Function names can be almost anything. However, the usage of names of existing functions should be avoided.
	
__Arguments__ 

* It is often useful to provide default values for arguments (_e.g._: `arg1=1:10`). This way they don't need to be provided in a function call. The argument list can also be left empty (`myfct <- function() { fct_body }`) if a function is expected to return always the same value(s). The argument `...` can be used to allow one function to pass on argument settings to another.

__Body__

* The actual expressions (commands/operations) are defined in the function body which should be enclosed by braces. The individual commands are separated by semicolons or new lines (preferred).

__Usage__ 

* Functions are called by their name followed by parentheses containing possible argument names. Empty parenthesis after the function name will result in an error message when a function requires certain arguments to be provided by the user. The function name alone will print the definition of a function.

__Scope__

* Variables created inside a function exist only for the life time of a function. Thus, they are not accessible outside of the function. To force variables in functions to exist globally, one can use the double assignment operator: `<<-` 

### Examples

__Define sample function__


```r
myfct <- function(x1, x2=5) { 
	z1 <- x1 / x1
	z2 <- x2 * x2
        myvec <- c(z1, z2) 
        return(myvec)
} 
```

__Function usage__


Apply function to values `2` and `5`

```r
myfct(x1=2, x2=5) 
```

```
## [1]  1 25
```

Run without argument names

```r
myfct(2, 5) 
```

```
## [1]  1 25
```

Makes use of default value `5`

```r
myfct(x1=2) 
```

```
## [1]  1 25
```
Print function definition (often unintended) 

```r
myfct 
```

```
## function(x1, x2=5) { 
## 	z1 <- x1 / x1
## 	z2 <- x2 * x2
##         myvec <- c(z1, z2) 
##         return(myvec)
## }
## <bytecode: 0x5865ccde4df0>
```

## Useful Utilities

### Debugging Utilities

Several debugging utilities are available for R. They include:

* `traceback`
* `browser`
* `options(error=recover)`
* `options(error=NULL)`
* `debug`

The [Debugging in R page](https://adv-r.hadley.nz/debugging.html) provides an overview of the available resources.

### Regular Expressions

R's regular expression utilities work similar as in other languages. To learn how to use them in R, one can consult the main help page on this topic with `?regexp`.

#### String matching with `grep`

The grep function can be used for finding patterns in strings, here letter `A` in vector `month.name`.

```r
month.name[grep("A", month.name)] 
```

```
## [1] "April"  "August"
```

#### String substitution with `gsub`

Example for using regular expressions to substitute a pattern by another one using a back reference. Remember: single escapes `\` need to be double escaped `\\` in R.


```r
gsub('(i.*a)', 'xxx_\\1', "virginica", perl = TRUE) 
```

```
## [1] "vxxx_irginica"
```

### Interpreting a Character String as Expression

Some useful examples

Generates vector of object names in session

```r
myfct <- function(x) x^2
mylist <- ls()
n <- which(mylist %in% "myfct")
mylist[n] 
```

```
## [1] "myfct"
```

Executes entry in position `n` as expression


```r
get(mylist[n])
```

```
## function(x) x^2
```

```r
get(mylist[n])(2)
```

```
## [1] 4
```

Alternative approach 

```r
eval(parse(text=mylist[n])) 
```

```
## function(x) x^2
```

### Replacement, Split and Paste Functions for Strings

__Selected examples__

Substitution with back reference which inserts in this example `_` character

```r
x <- gsub("(a)","\\1_", month.name[1], perl=T) 
x
```

```
## [1] "Ja_nua_ry"
```

Split string on inserted character from above

```r
strsplit(x,"_")
```

```
## [[1]]
## [1] "Ja"  "nua" "ry"
```

Reverse a character string by splitting first all characters into vector fields


```r
paste(rev(unlist(strsplit(x, NULL))), collapse="") 
```

```
## [1] "yr_aun_aJ"
```

### Time, Date and Sleep

__Selected examples__

Return CPU (and other) times that an expression used (here ls)

```r
system.time(ls()) 
```

```
##    user  system elapsed 
##       0       0       0
```

Return the current system date and time

```r
date() 
```

```
## [1] "Sat Feb 13 18:30:39 2021"
```

Pause execution of R expressions for a given number of seconds (e.g. in loop)

```r
Sys.sleep(1) 
```

#### Example

##### Import of Specific File Lines with Regular Expression

The following example demonstrates the retrieval of specific lines from an external file with a regular expression. First, an external file is created with the `cat` function, all lines of this file are imported into a vector with `readLines`, the specific elements (lines) are then retieved with the `grep` function, and the resulting lines are split into vector fields with `strsplit`.


```r
cat(month.name, file="zzz.txt", sep="\n")
x <- readLines("zzz.txt")
x[1:6] 
```

```
## [1] "January"  "February" "March"    "April"    "May"      "June"
```

```r
x <- x[c(grep("^J", as.character(x), perl = TRUE))]
t(as.data.frame(strsplit(x, "u")))
```

```
##                 [,1]  [,2] 
## c..Jan....ary.. "Jan" "ary"
## c..J....ne..    "J"   "ne" 
## c..J....ly..    "J"   "ly"
```
## Calling External Software

External command-line software can be called with `system`. The following example calls `blastall` from R

```r
system("blastall -p blastp -i seq.fasta -d uniprot -o seq.blastp")
```

## Running R Scripts

### Possibilities for Executing R Scripts

#### R console

```r
source("my_script.R")
```

#### Command-line


```sh
Rscript my_script.R # or just ./myscript.R after making it executable
R CMD BATCH my_script.R # Alternative way 1 
R --slave < my_script.R # Alternative way 2
```
#### Passing arguments from command-line to R

Create an R script named `test.R` with the following content:


```sh
myarg <- commandArgs()
print(iris[1:myarg[6], ])
```

Then run it from the command-line like this:

```sh
Rscript test.R 10
```

In the given example the number `10` is passed on from the command-line as an argument to the R script which is used to return to `STDOUT` the first 10 rows of the `iris` sample data. If several arguments are provided, they will be interpreted as one string and need to be split in R with the strsplit function. A more detailed example can be found [here](http://manuals.bioinformatics.ucr.edu/home/ht-seq\#TOC-Quality-Reports-of-FASTQ-Files-).


## Object-Oriented Programming (OOP)

R supports several systems for object-oriented programming (OOP). An older S3
system, and the more recently introduced R6 and S4 systems. The latter is the
most formal version that supports multiple inheritance, multiple dispatch and
introspection.  Many of these features are not available in the older S3
system. In general, the OOP approach taken by R is to separate the class
specifications from the specifications of generic functions (function-centric
system). The following introduction is restricted to the S4 system since it is
nowadays the preferred OOP method for package development in Bioconductor. More
information about OOP in R can be found in the following introductions: 

+ [Vincent Zoonekynd's introduction to S3 Classes](http://zoonek2.free.fr/UNIX/48_R/02.html#4)
+ [Christophe Genolini's S4 Intro](https://cran.r-project.org/doc/contrib/Genolini-S4tutorialV0-5en.pdf)
+ [Advanced Bioconductor Courses](http://master.bioconductor.org/help/course-materials/2008/advanced_R/)
+ [Programming with R by John Chambers](https://www.springer.com/gp/book/9780387759357) 
+ [R Programming for Bioinformatics by Robert Gentleman](http://www.bioconductor.org/help/publications/books/r-programming-for-bioinformatics/)
+ [Advanced R online book by Hadley Wichham](https://adv-r.hadley.nz/r6.html)

### Define S4 Classes

#### 1. Define S4 Classes with `setClass()` and `new()`


```r
y <- matrix(1:10, 2, 5) # Sample data set
setClass(Class="myclass",
    representation=representation(a="ANY"),
    prototype=prototype(a=y[1:2,]), # Defines default value (optional)
    validity=function(object) { # Can be defined in a separate step using setValidity
        if(class(object@a)[1]!="matrix") {
            return(paste("expected matrix, but obtained", class(object@a)))
        } else {
            return(TRUE)
        }
    }
)
```

The setClass function defines classes. Its most important arguments are

+ `Class`: the name of the class
+ `representation`: the slots that the new class should have and/or other classes that this class extends.
+ `prototype`: an object providing default data for the slots.
+ `contains`: the classes that this class extends.
+ `validity`, `access`, `version`: control arguments included for compatibility with S-Plus.
+ `where`: the environment to use to store or remove the definition as meta data.

#### 2. Create new class instance

The function `new` creates an instance of a class (here `myclass`).


```r
myobj <- new("myclass", a=y)
myobj
```

```
## An object of class "myclass"
## Slot "a":
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    3    5    7    9
## [2,]    2    4    6    8   10
```

If evaluated the following would return an error due to wrong input type (`data.frame` instead of `matrix`).


```r
new("myclass", a=iris) # Returns error due to wrong input  
```

The arguments of `new` are:

+ `Class`: the name of the class
+ `...`: data to include in the new object with arguments according to slots in class definition

#### 3. Initialization method

A more generic way of creating class instances is to define an initialization
method (more details below).


```r
setMethod("initialize", "myclass", function(.Object, a) {
    .Object@a <- a/a
    .Object
})
new("myclass", a = y)
```

```
## An object of class "myclass"
## Slot "a":
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    1    1    1    1
## [2,]    1    1    1    1    1
```

#### 4. Usage and helper functions

The '@' operator extracts the contents of a slot. Its usage should be limited to internal 
functions.


```r
myobj@a 
```

```
##      [,1] [,2] [,3] [,4] [,5]
## [1,]    1    3    5    7    9
## [2,]    2    4    6    8   10
```

Create a new S4 object from an old one.


```r
initialize(.Object=myobj, a=as.matrix(cars[1:2,])) 
```

```
## An object of class "myclass"
## Slot "a":
##   speed dist
## 1     1    1
## 2     1    1
```

If evaluated the `removeClass` function removes an object from the current session.
This does not apply to associated methods.


```r
removeClass("myclass") 
```

#### 5. Inheritance

Inheritance allows to define new classes that inherit all properties (e.g. data slots, methods) 
from their existing parent classes. The `contains` argument used below allows to extend 
existing classes. This propagates all slots of parent classes.


```r
setClass("myclass1", representation(a = "character", b = "character"))
setClass("myclass2", representation(c = "numeric", d = "numeric"))
setClass("myclass3", contains=c("myclass1", "myclass2"))
new("myclass3", a=letters[1:4], b=letters[1:4], c=1:4, d=4:1)
```

```
## An object of class "myclass3"
## Slot "a":
## [1] "a" "b" "c" "d"
## 
## Slot "b":
## [1] "a" "b" "c" "d"
## 
## Slot "c":
## [1] 1 2 3 4
## 
## Slot "d":
## [1] 4 3 2 1
```

```r
getClass("myclass1")
```

```
## Class "myclass1" [in ".GlobalEnv"]
## 
## Slots:
##                           
## Name:          a         b
## Class: character character
## 
## Known Subclasses: "myclass3"
```

```r
getClass("myclass2")
```

```
## Class "myclass2" [in ".GlobalEnv"]
## 
## Slots:
##                       
## Name:        c       d
## Class: numeric numeric
## 
## Known Subclasses: "myclass3"
```

```r
getClass("myclass3")
```

```
## Class "myclass3" [in ".GlobalEnv"]
## 
## Slots:
##                                               
## Name:          a         b         c         d
## Class: character character   numeric   numeric
## 
## Extends: "myclass1", "myclass2"
```

#### 6. Coerce objects to another class

The following defines a coerce method. After this the standard `as(..., "...")`
syntax can be used to coerce the new class to another one.














































