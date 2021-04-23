#####################################################
## Homework: Basic R Programming for Sequence Data ##
#####################################################
## Author: Thomas Girke
## Date: April 3, 2021
## Utility: Compute for one or many DNA sequences the 
##                (1) Reverse & Complement
##                (2) Translation into Protein

#################################################
## (A) REVERSE AND COMPLEMENT OF DNA SEQUENCES ##
#################################################

## (Task 1) Write a revComp() function that returns the reverse and complement of a
## DNA sequence string. Include an argument that will allow to return only the
## reversed sequence, the complemented sequence or the reversed and complemented
## sequence. The following R functions will be useful for the implementation:
##	   x <- c("ATGCATTGGACGTTAG")
##	   substring(x, 1:nchar(x), 1:nchar(x))
##	   rev(x)
##	   paste(x, collapse="")
##	   chartr("ATGC", "TACG", x)

## Define function
revComp <- function(x, type="RC") {
    ## Validity check of input
    ## Check length
    if(length(x) != 1) {
        x <- x[1]
        warning("Only the first sequence will be used.")
    }
    ## Enforces upper case in input
    x <- toupper(x) 
    
    ## Check for correct DNA alphabet 
    if(any(grepl("[^ATCG]", x))) stop("Sequence contains non-DNA characters.")

    ## Reverse
	if(type=="R") {
		x <- substring(x, 1:nchar(x), 1:nchar(x))
		x <- rev(x)
		x <- paste(x, collapse="")  	
    ## Complement
	} else if(type=="C") { 
                x <- chartr("ATGC", "TACG", x) 
    ## Reverse and Complement
	} else if(type=="RC") {
		x <- chartr("ATGC", "TACG", x) 
		x <- substring(x, 1:nchar(x), 1:nchar(x))
		x <- rev(x)
		x <- paste(x, collapse="")  	
	} else {
        stop("Argument type can only be assigned one of: C, R or RC.")
    }
	return(x)
}

## Usage
## (A.1) Creates random sample DNA sequences 
myDNA <- sapply(letters[1:10], function(x) paste(sample(c("A","T","G","C"), 21, replace=TRUE), collapse=""))

## (A.2) How to run revComp function 
revComp(x=myDNA[1], type="RC")
	# Option R: reverse 
	# Option C: complement
	# Option RC: reverse and complement

## (Task 2) Apply function to many sequences.
sapply(myDNA, function(i) revComp(x=i, type="RC"))

##############################################
## (B) TRANSLATION OF DNA INTO AA SEQUENCES ##
##############################################

## (Task 3) Write a function that will translate one or many DNA sequences in all three 
## reading frames into proteins. The following commands will simplify this task:
##          AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=TRUE, sep="\t")
##          AAv <- AAdf[,2]; names(AAv) <- AAdf[,1]
##          y <- gsub("(...)", "\\1_", y) # Inserts "_" after each triplet.
##          y <- unlist(strsplit(y, "_")) # Splits on "_" and returns vectorized triplets.
##          y <- y[grep("^...$", y)] # Removes incomplete triplets.
##          AAv[y] # Translation into protein by name-based subsetting.

## Define translate function
translateDNA <- function(x, frame=1) {
        ## Validity checks of inputs
        if(!any(1:3 %in% frame)) stop("Argument type can be assigned only one of these vaules: 1, 2 or 3.")
    
        # Enforce upper case in input
        x <- toupper(x) 
        
        ## Check for correct DNA alphabet 
        if(any(grepl("[^ATCG]", x))) stop("Sequence contains non-DNA characters.")
        
        ## Import genetic code and store it in named vector AAv
        AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=TRUE, sep="\t")
        AAv <- AAdf[,2]; names(AAv) <- AAdf[,1]

        ## Translate function
        translateOne <- function(x, frame) {
                ## Internal function for converting DNA string to vector of triplets 
        	tripSplit <- function(x=x) { 
	        	x <- gsub("(...)", "\\1_", x) 
	        	x <- unlist(strsplit(x, "_")) 
	        	x <- x[grep("^...$", x)] 
	        	return(x)	
	        }
                ## Perform translation on different frames
            if(frame==1) {
                        x <- tripSplit(x) 
	        }
	        else if(frame==2) {
		        x <- gsub("^.", "", x)
		        x <- tripSplit(x) 
	        }
	        else if(frame==3) {
		        x <- gsub("^..", "", x)
		        x <- tripSplit(x) 
	        } else {
                stop("The argument 'frame' can only be assined one of: 1, 2 or 3.")
            }
	        myAA <- AAv[x]
	        myAA <- paste(myAA, collapse="")
	        return(myAA)
        }
        if(length(x) == 1) return(translateOne(x, frame))
        if(length(x) > 1) return(sapply(x, function(i) translateOne(x=i, frame=frame)))
}

## Usage

## (B.1) How to run translateDNA function
mypep <- translateDNA(x=myDNA, frame=1)
		# The three possible reading frames (1, 2 or 3) can be selected under the frame option.

## (B.2) Return result in data frame
data.frame(DNA=names(mypep), PEP=as.vector(mypep))  

## (B.3) Return result as AAStringSet object
library(Biostrings)
AAStringSet(mypep)

## (B.4) Perform same operation on DNAStringSet object from Biostrings
translateDNA(DNAStringSet(myDNA))

## (B.5) Compare results
all(as.character(translateDNA(DNAStringSet(myDNA))) == mypep)


