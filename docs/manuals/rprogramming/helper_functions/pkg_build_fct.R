#' Talk to me function
#'
#' This function triggers response from your computer.
#' @param talk Are you in the mood to talk? Defaults to TRUE.
#' @keywords talk
#' @export
#' @examples
#' talkToMe()

talkToMe <- function(talk=TRUE){
    if(talk==TRUE){
        print("I am not in the mood to talk.")
    }
    else {
        print("I don't understand humans.")
    }
}

#' @title Row-wise summary stats
#'
#' @description Computes row-wise summary statistics for any combination of
#' columns of a `matrix`. The first argument (`myMA`)
#' specifies the input `matrix`; the second one (`group`) selects the columns
#' based on a grouping vector, and the third one the summary statistics (e.g. mean, sd, max).  
#' @param myMA PARAM_DESCRIPTION, Default: myMA
#' @param group PARAM_DESCRIPTION, Default: c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4)
#' @param myfct PARAM_DESCRIPTION, Default: mean
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'    myMA <- matrix(rnorm(100000), 10000, 10, dimnames=list(1:10000, paste("C", 1:10, sep="")))
#'    myMAcomp(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean)[1:2,] 
#'  }
#' }
#' @rdname myMAcomp
#' @export 

myMAcomp <- function(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean) {
	myList <- tapply(colnames(myMA), group, list)
	names(myList) <- sapply(myList, paste, collapse="_")
	myMAmean <- sapply(myList, function(x) apply(myMA[, x, drop=FALSE], 1, myfct))
	return(myMAmean)
}
## Usage:
# myMAcomp(myMA=myMA, group=c(1,1,1,2,2,2,3,3,4,4), myfct=mean)[1:2,] 
