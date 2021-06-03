##############################################################
## Run baySeq on sample comparisons defined in targets file ##
##############################################################

## Load libraries, and targets and read counts table obtained from RNA-Seq workflow
library(baySeq); library(systemPipeR)
targets <- read.delim("targetsPE.txt", comment.char="#") # Imports targets file used here: https://bit.ly/3fMB5hm
samples <- as.character(targets$Factor)
names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
group <- as.character(samples)
counts <- read.delim("results/countDFeByg.xls", row.names=1, comment.char="#") # Imports read count table generated here: https://bit.ly/3fMB5hm
counts <- counts[, names(samples)]
counts[is.na(counts)] <- 0

## Run baySeq for individual comparisons defined in header of targets file
cmp <- systemPipeR::readComp(file="targetsPE.txt", format="matrix", delim="-")[[1]]
rownames(cmp) <- paste(cmp[,1], cmp[,2], sep="-")
run_baySeq <- function(counts, cmp=cmp, samplesize=1000){
    group <- samples[samples %in% cmp]
    counts_sub <- counts[, names(group)]
    group_list <- list(NDE=rep(1, length(group)), DE=rep(1:2, table(group)))
    cd_ob <- new("countData", data=as.matrix(counts_sub), replicates=group, groups=group_list)
    densityFunction(cd_ob) <- nbinomDensity
    libsizes(cd_ob) <- getLibsizes(cd_ob) 
    cd_ob@annotation <- data.frame(rownames(counts_sub))   
    cd_ob <- getPriors.NB(cd_ob, samplesize=samplesize, consensus=FALSE, cl=NULL)
    cd_ob <- getLikelihoods(cd_ob, pET="BIC", cl=NULL)
    resultDF <- topCounts(cd_ob, group="DE", number=nrow(counts))  
    rownames(resultDF) <- resultDF[,1]
    resultDF <- resultDF[,-1]
    ## baySeq lacks LFC output. To allow filtering, as with other DEG methods they are computed here in a naive manner
    source("https://raw.githubusercontent.com/tgirke/GEN242/main/content/en/tutorials/rpackages/helper_functions/pkg_build_fct.R")
    myMA <- resultDF[,names(group)]
    meanma <- myMAcomp(myMA=myMA, group=group, myfct=mean)
    lfc <- log2(meanma[,1]) - log2(meanma[,2])
    resultDF <- cbind(resultDF, logFC=lfc, FDR=resultDF$'FDR.DE') 
    return(resultDF)
}

## Run for single comparison
# single_comparison <- run_baySeq(counts, cmp=cmp["M1-A1",], samplesize=1000)
## Run for all comparisons and store results in list where each result is named by comparison
## For debugging/testing reduce count matrix and sample size 
all_comparisons <- sapply(rownames(cmp), function(x) run_baySeq(counts=counts[1:1000,], cmp=cmp[x,], samplesize=10), simplify=FALSE)
## For full run use this line instead!!!
# all_comparisons <- sapply(rownames(cmp), function(x) run_baySeq(counts=counts, cmp=cmp[x,], samplesize=1000), simplify=FALSE)
names(all_comparisons)
sapply(all_comparisons, dim)

## Function to assemble list in master DEG data.frame
assembleFromList <- function(x=all_comparisons, counts) {
    masterDEG <- data.frame(row.names = rownames(counts))
    for(i in seq_along(x)) {
        singleDF <- x[[i]]
        fix_colnames <- colnames(singleDF)
        # fix_colnames[fix_colnames %in% c("logFC", "adj.P.Val")] <- c("logFC", "FDR")
        fix_colnames <- paste(paste(names(x)[i], collapse = "-"), fix_colnames, sep = "_")
        colnames(singleDF) <- fix_colnames
        masterDEG <- cbind(masterDEG, singleDF[rownames(masterDEG),])
    }
    masterDEGdf <- na.omit(masterDEG)
    return(masterDEG)
}
masterDEGdf <- assembleFromList(x=all_comparisons, counts=counts)
write.table(masterDEGdf, "./results/bayseq_allcomp.xls", quote = FALSE, sep = "\t", col.names = NA)

## Plot DEG results
bayseqDF <- read.delim("results/bayseq_allcomp.xls", row.names = 1, check.names = FALSE)
DEG_list <- filterDEGs(degDF=bayseqDF, filter = c(Fold = 2, FDR = 10))

