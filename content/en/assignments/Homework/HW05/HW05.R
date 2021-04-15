########################
## HW5: R Programming ##
########################

## A. Reverse and complement of DNA

## Task 1: Write a `RevComp` function that returns the reverse and complement
## of a DNA sequence string. Include an argument that will allow to return only
## the reversed sequence, the complemented sequence or the reversed and
## complemented sequence. The following R functions will be useful for the
## implementation: 

## Sample sequence
x <- c("ATGCATTGGACGTTAG")  
x

## Vectorize sequence
x <- substring(x, 1:nchar(x), 1:nchar(x)) 
x

## Form reverse of sequence
x <- rev(x) 
x

## Convert back to character string
x <- paste(x, collapse="")
x

## Form complement of sequence
chartr("ATGC", "TACG", x) 


## Task 2: Write a function that applies the `RevComp` function to many sequences stored in a vector.

## B. Translate DNA into Protein

## Task 3: Write a function that will translate one or many DNA sequences in
## all three reading frames into proteins. The following commands will simplify
## this task:

## Import lookup table for genetic code
AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=TRUE, sep="\t") 
AAdf[1:4,]

## Generate named character vector from relevant columns (translation vector)
AAv <- as.character(AAdf[,2]) 
names(AAv) <- AAdf[,1] 
AAv

## Tripletize DNA sequence
y <- gsub("(...)", "\\1_", x) 
y <- unlist(strsplit(y, "_")) 
y <- y[grep("^...$", y)] 

## Translate to protein by passing on triplets to translation vector AAv
AAv[y] 


## Homework submission

## The tasks 1-3 of this homework can be summarized as follows: submit the 3
## functions in one well structured and annotated R script. The script should
## include instructions on how to use the functions.






