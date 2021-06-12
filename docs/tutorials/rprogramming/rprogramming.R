## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))


## ----if_statement, eval=FALSE---------------------------------------------------------------------
## if (TRUE) {
##     statements_1
## } else {
##     statements_2
## }


## ----if_statement_example, eval=TRUE--------------------------------------------------------------
if (1==0) { 
    print(1) 
} else { 
    print(2) 
}


## ----ifelse_statement, eval=FALSE-----------------------------------------------------------------
## ifelse(test, true_value, false_value)


## ----ifelse_statement_example, eval=TRUE----------------------------------------------------------
x <- 1:10 
ifelse(x<5, sqrt(x), 0)


## ----for_loop, eval=FALSE-------------------------------------------------------------------------
## for(variable in sequence) {
## 	statements
## }


## ----for_loop_example, eval=TRUE------------------------------------------------------------------
mydf <- iris
myve <- NULL
for(i in seq(along=mydf[,1])) {
	myve <- c(myve, mean(as.numeric(mydf[i,1:3])))
}
myve[1:8]


## ----for_loop_inject_example, eval=TRUE-----------------------------------------------------------
myve <- numeric(length(mydf[,1]))
for(i in seq(along=myve)) {
	myve[i] <- mean(as.numeric(mydf[i,1:3]))
}
myve[1:8]


## ----for_loop_stop_example, eval=FALSE------------------------------------------------------------
## x <- 1:10
## z <- NULL
## for(i in seq(along=x)) {
## 	if (x[i] < 5) {
## 		z <- c(z, x[i]-1)
## 	} else {
## 		stop("values need to be < 5")
## 	}
## }


## ----while_loop, eval=FALSE-----------------------------------------------------------------------
## while(condition) {
## 	statements
## }


## ----while_loop_example, eval=TRUE----------------------------------------------------------------
z <- 0
while(z<5) { 
	z <- z + 2
	print(z)  
}


## ----apply_loop, eval=FALSE-----------------------------------------------------------------------
## apply(X, MARGIN, FUN, ARGs)


## ----apply_loop_example, eval=TRUE----------------------------------------------------------------
apply(iris[1:8,1:3], 1, mean)


## ----tapply_loop, eval=FALSE----------------------------------------------------------------------
## tapply(vector, factor, FUN)


## ----tapply_loop_example, eval=TRUE---------------------------------------------------------------
iris[1:2,]
tapply(iris$Sepal.Length, iris$Species, mean)


## ----lapply_loop_example, eval=TRUE---------------------------------------------------------------
l <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
lapply(l, mean)
sapply(l, mean)
vapply(l, mean, FUN.VALUE=numeric(1))


## ----lapply_loop_fct_example, eval=FALSE----------------------------------------------------------
## lapply(names(l), function(x) mean(l[[x]]))
## sapply(names(l), function(x) mean(l[[x]]))
## vapply(names(l), function(x) mean(l[[x]]), FUN.VALUE=numeric(1))


## ----for_loop_with_c_append, eval=FALSE-----------------------------------------------------------
## myMA <- matrix(rnorm(1000000), 100000, 10, dimnames=list(1:100000, paste("C", 1:10, sep="")))
## results <- NULL
## system.time(for(i in seq(along=myMA[,1])) results <- c(results, mean(myMA[i,])))
##    user  system elapsed
##  39.156   6.369  45.559


## ----for_loop_with_inject, eval=FALSE-------------------------------------------------------------
## results <- numeric(length(myMA[,1]))
## system.time(for(i in seq(along=myMA[,1])) results[i] <- mean(myMA[i,]))
##    user  system elapsed
##   1.550   0.005   1.556


## ----apply_loop_mean, eval=FALSE------------------------------------------------------------------
## system.time(myMAmean <- apply(myMA, 1, mean))
##   user  system elapsed
##  1.452   0.005   1.456
## 
## system.time(myMAmean <- rowMeans(myMA))
##    user  system elapsed
##   0.005   0.001   0.006


## ----apply_loop_vs_vectorized, eval=FALSE---------------------------------------------------------
## system.time(myMAsd <- apply(myMA, 1, sd))
##    user  system elapsed
##   3.707   0.014   3.721
## myMAsd[1:4]
##         1         2         3         4
## 0.8505795 1.3419460 1.3768646 1.3005428
## 
## system.time(myMAsd <- sqrt((rowSums((myMA-rowMeans(myMA))^2)) / (length(myMA[1,])-1)))
##    user  system elapsed
##   0.020   0.009   0.028
## 
## myMAsd[1:4]
##         1         2         3         4
## 0.8505795 1.3419460 1.3768646 1.3005428


## ----create_stats_ma, eval=TRUE-------------------------------------------------------------------
lfcPvalMA <- function(Nrow=200, Ncol=4, stats_labels=c("lfc", "pval")) {
    set.seed(1410)
    assign(stats_labels[1], runif(n = Nrow * Ncol, min = -4, max = 4))
    assign(stats_labels[2], runif(n = Nrow * Ncol, min = 0, max = 1))
    lfc_ma <- matrix(lfc, Nrow, Ncol, dimnames=list(paste("g", 1:Nrow, sep=""), paste("t", 1:Ncol, "_", stats_labels[1], sep=""))) 
    pval_ma <- matrix(pval, Nrow, Ncol, dimnames=list(paste("g", 1:Nrow, sep=""), paste("t", 1:Ncol, "_", stats_labels[2], sep=""))) 
    statsMA <- cbind(lfc_ma, pval_ma)
    return(statsMA[, order(colnames(statsMA))])
}
degMA <- lfcPvalMA(Nrow=200, Ncol=4, stats_labels=c("lfc", "pval"))
dim(degMA) 
degMA[1:4,] # Prints first 4 rows of DEG matrix generated as a test data set 


## ----stats_ma_to_list, eval=TRUE------------------------------------------------------------------
degList <- list(lfc=degMA[ , grepl("lfc", colnames(degMA))], pval=degMA[ , grepl("pval", colnames(degMA))])
names(degList)
sapply(degList, dim)


## ----filter_stats, eval=TRUE----------------------------------------------------------------------
queryResult <- (degList$lfc <= 1 | degList$lfc <= -1) & degList$pval <= 0.5 
colnames(queryResult) <- gsub("_.*", "", colnames(queryResult)) # Adjust column names 
queryResult[1:4,]


## ----id_result_list1, eval=TRUE-------------------------------------------------------------------
matchingIDlist <- sapply(colnames(queryResult), function(x) names(queryResult[queryResult[ , x] , x]), simplify=FALSE)
matchingIDlist


## ----id_result_list2, eval=TRUE-------------------------------------------------------------------
matchingID <- rowSums(queryResult) > 2 
queryResult[matchingID, , drop=FALSE]
names(matchingID[matchingID])


## ----function_def_syntax, eval=FALSE--------------------------------------------------------------
## myfct <- function(arg1, arg2, ...) {
## 	function_body
## }


## ----function_call_syntax, eval=FALSE-------------------------------------------------------------
## myfct(arg1=..., arg2=...)


## ----define_function_example, eval=TRUE-----------------------------------------------------------
myfct <- function(x1, x2=5) { 
	z1 <- x1 / x1
	z2 <- x2 * x2
        myvec <- c(z1, z2) 
        return(myvec)
} 


## ----usage_function_example1, eval=TRUE-----------------------------------------------------------
myfct(x1=2, x2=5) 


## ----usage_function_example2, eval=TRUE-----------------------------------------------------------
myfct(2, 5) 


## ----usage_function_example3, eval=TRUE-----------------------------------------------------------
myfct(x1=2) 


## ----usage_function_example4, eval=TRUE-----------------------------------------------------------
myfct 


## ----grep_fct, eval=TRUE--------------------------------------------------------------------------
month.name[grep("A", month.name)] 


## ----gsub_fct, eval=TRUE--------------------------------------------------------------------------
gsub('(i.*a)', 'xxx_\\1', "virginica", perl = TRUE) 


## ----ls_fct, eval=TRUE----------------------------------------------------------------------------
myfct <- function(x) x^2
mylist <- ls()
n <- which(mylist %in% "myfct")
mylist[n] 


## ----eval_expr, eval=TRUE-------------------------------------------------------------------------
get(mylist[n])
get(mylist[n])(2)


## ----eval_expr2, eval=TRUE------------------------------------------------------------------------
eval(parse(text=mylist[n])) 


## ----back_ref, eval=TRUE--------------------------------------------------------------------------
x <- gsub("(a)","\\1_", month.name[1], perl=T) 
x


## ----split_string, eval=TRUE----------------------------------------------------------------------
strsplit(x,"_")


## ----reverse_string, eval=TRUE--------------------------------------------------------------------
paste(rev(unlist(strsplit(x, NULL))), collapse="") 


## ----sys_time, eval=TRUE--------------------------------------------------------------------------
system.time(ls()) 


## ----sys_date, eval=TRUE--------------------------------------------------------------------------
date() 


## ----sys_sleep, eval=TRUE-------------------------------------------------------------------------
Sys.sleep(1) 


## ----read_lines, eval=TRUE------------------------------------------------------------------------
cat(month.name, file="zzz.txt", sep="\n")
x <- readLines("zzz.txt")
x[1:6] 
x <- x[c(grep("^J", as.character(x), perl = TRUE))]
t(as.data.frame(strsplit(x, "u")))


## ----system_blast, eval=FALSE---------------------------------------------------------------------
## system("blastall -p blastp -i seq.fasta -d uniprot -o seq.blastp")


## ----r_script1, eval=FALSE------------------------------------------------------------------------
## source("my_script.R")


## Rscript my_script.R # or just ./myscript.R after making it executable

## R CMD BATCH my_script.R # Alternative way 1

## R --slave < my_script.R # Alternative way 2


## myarg <- commandArgs()

## print(iris[1:myarg[6], ])


## Rscript test.R 10


## ----define_s4, eval=TRUE-------------------------------------------------------------------------
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


## ----new_s4, eval=TRUE----------------------------------------------------------------------------
myobj <- new("myclass", a=y)
myobj


## ----new_s4_error, eval=FALSE---------------------------------------------------------------------
## new("myclass", a=iris) # Returns error due to wrong input


## ----s4_init_method, eval=TRUE--------------------------------------------------------------------
setMethod("initialize", "myclass", function(.Object, a) {
    .Object@a <- a/a
    .Object
})
new("myclass", a = y)


## ----s4_helper_fct1, eval=TRUE--------------------------------------------------------------------
myobj@a 


## ----s4_helper_fct2, eval=TRUE--------------------------------------------------------------------
initialize(.Object=myobj, a=as.matrix(cars[1:2,])) 


## ----s4_helper_fct3, eval=FALSE-------------------------------------------------------------------
## removeClass("myclass")


## ----s4_inheritance1, eval=TRUE-------------------------------------------------------------------
setClass("myclass1", representation(a = "character", b = "character"))
setClass("myclass2", representation(c = "numeric", d = "numeric"))
setClass("myclass3", contains=c("myclass1", "myclass2"))
new("myclass3", a=letters[1:4], b=letters[1:4], c=1:4, d=4:1)
getClass("myclass1")
getClass("myclass2")
getClass("myclass3")


## ----s4_coerce, eval=TRUE-------------------------------------------------------------------------
setAs(from="myclass", to="character", def=function(from) as.character(as.matrix(from@a)))
as(myobj, "character")


## ----s4_virtual, eval=TRUE------------------------------------------------------------------------
setClass("myVclass")
setClass("myVclass", representation(a = "character", "VIRTUAL"))


## ----s4_setgeneric, eval=TRUE---------------------------------------------------------------------
setGeneric(name="acc", def=function(x) standardGeneric("acc"))
setMethod(f="acc", signature="myclass", definition=function(x) {
	return(x@a)
})
acc(myobj)


## ----s4_replace_acc, eval=TRUE--------------------------------------------------------------------
setGeneric(name="acc<-", def=function(x, value) standardGeneric("acc<-"))
setReplaceMethod(f="acc", signature="myclass", definition=function(x, value) {
				 x@a <- value
                 return(x)
})
## After this the following replace operations with 'acc' work on new object class
acc(myobj)[1,1] <- 999 # Replaces first value
colnames(acc(myobj)) <- letters[1:5] # Assigns new column names
rownames(acc(myobj)) <- letters[1:2] # Assigns new row names
myobj


## ----s4_replace_bracket, eval=TRUE----------------------------------------------------------------
setReplaceMethod(f="[", signature="myclass", definition=function(x, i, j, value) {
				 x@a[i,j] <- value
                 return(x)
})
myobj[1,2] <- 999
myobj


## ----s4_bracket_subsetting, eval=TRUE-------------------------------------------------------------
setMethod(f="[", signature="myclass",
		  definition=function(x, i, j, ..., drop) {
          x@a <- x@a[i,j]
          return(x)
})
myobj[1:2, 1:3] # Standard subsetting works now on new class


## ----s4_printing, eval=TRUE-----------------------------------------------------------------------
setMethod(f="show", signature="myclass", definition=function(object) {
		  cat("An instance of ", "\"", class(object), "\"", " with ", length(acc(object)[,1]), " elements", "\n", sep="")
		  if(length(acc(object)[,1])>=5) {
				print(as.data.frame(rbind(acc(object)[1:2,], ...=rep("...", length(acc(object)[1,])), acc(object)[(length(acc(object)[,1])-1):length(acc(object)[,1]),])))
		  } else {
				print(acc(object))
		  }
})
myobj # Prints object with custom method


## ----s4_custom_methods, eval=TRUE-----------------------------------------------------------------
setGeneric(name="randomize", def=function(x) standardGeneric("randomize"))
setMethod(f="randomize", signature="myclass", definition=function(x) {
		  acc(x)[sample(1:length(acc(x)[,1]), length(acc(x)[,1])), ]
})
randomize(myobj)


## ----s4_plot_methods, eval=TRUE-------------------------------------------------------------------
setMethod(f="plot", signature="myclass", definition=function(x, ...) {
		  barplot(as.matrix(acc(x)), ...)
})
plot(myobj)


## ----exercise1_for, eval=TRUE---------------------------------------------------------------------
myMA <- matrix(rnorm(500), 100, 5, dimnames=list(1:100, paste("C", 1:5, sep="")))
myve_for <- NULL
for(i in seq(along=myMA[,1])) {
	myve_for <- c(myve_for, mean(as.numeric(myMA[i, ])))
}
myResult <- cbind(myMA, mean_for=myve_for)
myResult[1:4, ]


## ----exercise1_while, eval=TRUE-------------------------------------------------------------------
z <- 1
myve_while <- NULL
while(z <= length(myMA[,1])) {
	myve_while <- c(myve_while, mean(as.numeric(myMA[z, ])))
	z <- z + 1
}
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while)
myResult[1:4, -c(1,2)]


## ----exercise1_confirm, eval=TRUE-----------------------------------------------------------------
all(myResult[,6] == myResult[,7])


## ----exercise1_apply, eval=TRUE-------------------------------------------------------------------
myve_apply <- apply(myMA, 1, mean)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply)
myResult[1:4, -c(1,2)]


## ----exercise1_noloops, eval=TRUE-----------------------------------------------------------------
mymean <- rowMeans(myMA)
myResult <- cbind(myMA, mean_for=myve_for, mean_while=myve_while, mean_apply=myve_apply, mean_int=mymean)
myResult[1:4, -c(1,2,3)]


## ----exercise2_fct, eval=TRUE---------------------------------------------------------------------
myMA <- matrix(rnorm(100000), 10000, 10, dimnames=list(1:10000, paste("C", 1:10, sep="")))
myMA[1:2,]
myList <- tapply(colnames(myMA), c(1,1,1,2,2,2,3,3,4,4), list) 
names(myList) <- sapply(myList, paste, collapse="_")
myMAmean <- sapply(myList, function(x) apply(myMA[,x], 1, mean))
myMAmean[1:4,] 


## ----exercise2_fct_solution, eval=FALSE, echo=FALSE, keep.source=TRUE-----------------------------
## myMAcomp <- function(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean) {
## 	myList <- tapply(colnames(myMA), group, list)
## 	names(myList) <- sapply(myList, paste, collapse="_")
## 	myMAmean <- sapply(myList, function(x) apply(myMA[, x, drop=FALSE], 1, myfct))
## 	return(myMAmean)
## }
## myMAcomp(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean)[1:2,]


## ----nested_loops1, eval=TRUE---------------------------------------------------------------------
setlist <- lapply(11:30, function(x) sample(letters, x, replace=TRUE))
names(setlist) <- paste("S", seq(along=setlist), sep="") 
setlist[1:6]


## ----nested_loops2, eval=TRUE---------------------------------------------------------------------
setlist <- sapply(setlist, unique)
olMA <- sapply(names(setlist), function(x) sapply(names(setlist), 
               function(y) sum(setlist[[x]] %in% setlist[[y]])))
olMA[1:12,] 


## ----nested_loops3, eval=TRUE---------------------------------------------------------------------
library(pheatmap); library("RColorBrewer")
pheatmap(olMA, color=brewer.pal(9,"Blues"), cluster_rows=FALSE, cluster_cols=FALSE, display_numbers=TRUE, number_format="%.0f", fontsize_number=10)
# image(olMA) 


## ----package_skeleton2, eval=FALSE----------------------------------------------------------------
## package.skeleton(name="mypackage", code_files=c("script1.R"))


## ----build_package_tar, eval=FALSE----------------------------------------------------------------
## system("R CMD build mypackage")


## ----install_package_tar, eval=FALSE--------------------------------------------------------------
## install.packages("mypackage_1.0.tar.gz", repos=NULL, type="source")
## library(mypackage)
## ?myMAcomp # Opens help for function defined by mypackage


## ----pattern_matching, eval=FALSE-----------------------------------------------------------------
## source("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rprogramming/scripts/patternSearch.R")


## ----word_finder, eval=FALSE----------------------------------------------------------------------
## source("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rprogramming/scripts/wordFinder.R")


## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

