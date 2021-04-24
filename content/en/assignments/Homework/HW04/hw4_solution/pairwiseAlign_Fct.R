#############################################################
## R code for creating dynamic programming matrix for HW4B ##
#############################################################
## Author: Thomas Girke
## Last update: Apr 3, 2019

###########################################
## Generate Dynamic Programming Matrices ##
###########################################
dynProgMatrix <- function(S1, S2, align_method="global", gap_penalty=8, substitutionMA="BLOSUM50") { 
    ## Validity checks
    if(!align_method %in% c("local", "global")) stop("'local' or 'global' expected in 'align_method'")
    # More input validity checks should be added ...
    
    ## Vectorize input sequences
    x <- substring(S1, 1:nchar(S1), 1:nchar(S1))
    y <- substring(S2, 1:nchar(S2), 1:nchar(S2))
    
    ## Load substitution matrix from Biostrings package
    library(Biostrings)
    ssMA <- data(list=substitutionMA)

    ## BLOSUM50 matrix subsetted to sequences
    s_ma <- get(ssMA)[y, x]

    ## Gap penality
    gp <- gap_penalty 

    ## Create global dynamic programming matrix based on x, y and gp
    if(align_method=="global") {
        glob_ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
        glob_ma[1,] <- seq(0, -(length(glob_ma[1,])-1) * gp, -gp)
        glob_ma[,1] <- seq(0, -(length(glob_ma[,1])-1) * gp, -gp)
        glob_ma

        ## Store path
        glob_ma_path <- glob_ma

        ## Global dyn programming matrix
        for(j in 2:length(glob_ma[1,])) {
            for(i in 2:length(glob_ma[,1])) {
                diagonal <- s_ma[i-1,j-1] + glob_ma[i-1,j-1]
                left <- glob_ma[i,j-1]-gp
                up <- glob_ma[i-1,j]-gp
                glob_ma[i, j] <- max(diagonal, up, left)
                if(max(diagonal, up, left) == diagonal) {
                    glob_ma_path[i,j] <- 1
                } else { 
                    if(max(diagonal, up, left) == up) { 
                        glob_ma_path[i,j] <- 2
                } else {
                    glob_ma_path[i,j] <- 3}
                }
                i <- i+1
            }  
            j <- j+1
        }
        ## Assign corresponding paths to boundaries
        glob_ma_path[1,-1] <- 3
        glob_ma_path[-1,1] <- 2
        
        ## Return results in list
        return(list(glob_ma=glob_ma, glob_ma_path=glob_ma_path))
    }

    ## Create local dynamic programming matrix based on x, y and gp
    if(align_method=="local") {
        loc_ma <- matrix(NA, length(y)+1, length(x)+1, dimnames=list(c("gp", y), c("gp", x)))
        loc_ma[1,] <- 0; loc_ma[,1] <- 0
        loc_ma

        ## Store path
        loc_ma_path <- loc_ma

        ## Local dyn programming matrix
        for(j in 2:length(loc_ma[1,])) {
            for(i in 2:length(loc_ma[,1])) {
                diagonal = s_ma[i-1,j-1] + loc_ma[i-1,j-1]
                left = loc_ma[i,j-1]-gp
                up = loc_ma[i-1,j]-gp
                loc_ma[i,j] = max(diagonal, up, left, 0)
                if(max(diagonal, up, left) == diagonal) {
                    loc_ma_path[i,j] = 1
                } else {
                    if(max(diagonal, up, left) == up) {
                        loc_ma_path[i,j] = 2
                    } else {
                        loc_ma_path[i,j] = 3
                    }
                }
                i=i+1
            }  
            j=j+1
        }

        ## Return results in list
        return(list(loc_ma=loc_ma, loc_ma_path=loc_ma_path))
    }
}

## Usage:

## Define input sequences
# S1 <- "PFGFGKRSCMGRRLA"
# S2 <- "FIPFSAGPRNCIGQK"

# S1 <- "HEAGAWGHEE"
# S2 <- "PAWHEAE"

# dynMA <- dynProgMatrix(S1, S2, align_method="local", gap_penalty=8, substitutionMA="BLOSUM50")

##################################
## Traceback Alignment Function ##
##################################
alignmentTraceback <- function(ma, ma_path, align_method) {
    ## Validity checks
    if(!align_method %in% c("local", "global")) stop("'local' or 'global' expected in 'align_method'")
    # More input validity checks should be added ...

    ## Save original input matrices for return object
    orig_ma <- ma
    orig_ma_path <- ma_path

    ## For local alignment subset matrix to max value coordinate 
    if(align_method=="local") {
        traceback_coor <- which(ma == max(ma), arr.ind=TRUE)
        ma <- ma[1:traceback_coor[1,1], 1:traceback_coor[1,2]]
        ma_path <- ma_path[1:traceback_coor[1,1], 1:traceback_coor[1,2]]
    } 
    s1 <- colnames(ma_path)[-1]
    as1 <- NULL
    s2 <- rownames(ma_path)[-1]
    as2 <- NULL
    xc <- length(ma_path[1,]) 
    yc <- length(ma_path[,2])
    path_coor <- matrix(NA, nrow=length(c(s1, s2)), ncol = 2, dimnames=list(NULL, c("x", "y")))
    counter <- 1
    path_coor[counter, ] <- c(xc, yc) # handles intial bottom right corner coordinates 
    counter <- counter + 1 # now counter needs to be increased by 1 
    ## The following if statatemt is to use same while loop condition for global and local
    if(align_method=="global") {
        ma[ma==0] <- 0.9; ma[1,1] <- 0 
    }
    ## Store score of alignment
    my_score <- ma[nrow(ma), ncol(ma)]
    
    ## Creat alignment with traceback
    while(ma[yc, xc] != 0) {
        ## Align residues, or gap in s1 or s2
        if(ma_path[yc, xc] == 1) { 
            as1 <- c(as1, s1[length(s1)])
            s1 <- s1[-length(s1)]
            as2 <- c(as2, s2[length(s2)])
            s2 <- s2[-length(s2)]
            direction <- 1 
        } else if(ma_path[yc, xc] == 2) {
            as1 <- c(as1, "-")
            as2 <- c(as2, s2[length(s2)])
            s2 <- s2[-length(s2)]
            direction <- 2 
        } else if(ma_path[yc, xc] == 3) {
            as1 <- c(as1, s1[length(s1)])
            s1 <- s1[-length(s1)]
            as2 <- c(as2, "-")
            direction <- 3 
        } else {
            stop("Exception in ma_path")
        }
        ## Update coordinates
        if(direction == 1) { 
            xc <- xc - 1
            yc <- yc - 1
        } else if(direction == 2) {
            yc <- yc - 1
        } else if(direction == 3) {
            xc <- xc - 1
        } else {
            stop("Exception in coordinates")
        }
        path_coor[counter, ] <- c(xc, yc)
        counter <- counter + 1
    }
    ## Remove NA rows in path_coor
    path_coor <- na.omit(path_coor)
    attributes(path_coor)$na.action <- NULL
    ## Return alignment
    alignList <- list(
                      ma=orig_ma,
                      ma_path=orig_ma_path,
                      path_coor=path_coor,
                      as1=paste(rev(as1), collapse=""),
                      consensus=paste(rev(ifelse(as1==as2, "|", " ")), collapse=""),
                      as2=paste(rev(as2), collapse=""),
                      score=my_score
    )
    return(alignList)
}
## Usage: 
# alignList <- alignmentTraceback(ma=dynMA[[1]], ma_path=dynMA[[2]], align_method="local") 

## Print alignment and score
printAlign <- function(x) {
    x <- cat("\n", "S1: ", alignList$as1, "\n", 
        "    ", alignList$consensus, "\n", 
        "S2: ", alignList$as2, "\n\n",
        "Score of alignment:", alignList$score, "\n")
}
## Usage:
# printAlign(x=alignList)

###################################################################
## Print Dynamic programming matrix with traceback path in color ##
###################################################################
printColMa <- function(alignList) {
    library(kableExtra) # See vignette here: https://bit.ly/3sauMY8
    ma <- alignList$ma
    pc <- alignList$path_coor
    ima <- ma; ima[,] <- 0
    for(i in seq_along(pc[,1])) ima[pc[i,2], pc[i,1]] <- 1 
    ma[ima==1] <- cell_spec(ma[ima==1], color = "red")
    kbl(ma, escape = FALSE) %>% 
        kable_paper("striped", full_width = FALSE) 
}
## Usage:
# printColMa(alignList)



