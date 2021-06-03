##################################################################
## Run limma voom on sample comparisons defined in targets file ##
##################################################################

## (A) Run sample comparisons independently by subsetting count matrix to replicates
## that are part of a comparison 

## Load libraries, and targets and read counts table obtained from RNA-Seq workflow
## Load packages and data
library(limma); library(edgeR); library(systemPipeR)
targets <- read.delim("targetsPE.txt", comment.char="#")
samples <- as.character(targets$Factor)
names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
group <- as.character(samples)
counts <- read.delim("results/countDFeByg.xls", row.names=1, comment.char="#")
counts <- counts[, names(samples)]
counts[is.na(counts)] <- 0

## Run limma voom for individual comparisons defined in header of targets file
cmp <- systemPipeR::readComp(file="targetsPE.txt", format="matrix", delim="-")[[1]]
rownames(cmp) <- paste(cmp[,1], cmp[,2], sep="-")
run_voom <- function(counts, cmp) {
	group <- samples[samples %in% cmp]
	design <- model.matrix(~group)
	counts_sub <- counts[,names(group)]
	dge <- DGEList(counts=counts_sub, group=group)
	keep <- filterByExpr(dge, design)
	dge <- dge[keep,,keep.lib.sizes=FALSE]
	dge <- calcNormFactors(dge)
	v <- voom(dge, design, plot=FALSE)
	fit <- lmFit(v, design)
	fit <- eBayes(fit)
	topTable(fit, coef=ncol(design), nrow(counts))
}
## Run for single comparison
single_comparison <- run_voom(counts=counts, cmp=cmp["M1-A1",])
## Run for all comparisons and store results in list where each result is named by comparison
all_comparisons <- sapply(rownames(cmp), function(x) run_voom(counts=counts, cmp=cmp[x,]), simplify=FALSE)
names(all_comparisons)
sapply(all_comparisons, dim)

## Function to assemble list in master DEG data.frame
assembleFromList <- function(x=all_comparisons, counts) {
    masterDEG <- data.frame(row.names = rownames(counts))
    for(i in seq_along(x)) {
        singleDF <- x[[i]]
        fix_colnames <- colnames(singleDF)
        fix_colnames[fix_colnames %in% c("logFC", "adj.P.Val")] <- c("logFC", "FDR")
        fix_colnames <- paste(paste(names(x)[i], collapse = "-"), fix_colnames, sep = "_")
        colnames(singleDF) <- fix_colnames
        masterDEG <- cbind(masterDEG, singleDF[rownames(masterDEG),])
    }
    masterDEGdf <- na.omit(masterDEG)
    return(masterDEG)
}
masterDEGdf <- assembleFromList(x=all_comparisons, counts=counts)
write.table(masterDEGdf, "./results/limmavoom_allcomp.xls", quote = FALSE, sep = "\t", col.names = NA)

## Plot DEG results
voomDF <- read.delim("results/limmavoom_allcomp.xls", row.names = 1, check.names = FALSE)
DEG_list <- filterDEGs(degDF = voomDF, filter = c(Fold = 2, FDR = 50))

## (B) Run sample comparisons with all replicates included. In this case
## both a design and contrast matrix are used.

## Load libraries, and targets and read counts table obtained from RNA-Seq workflow
## Load packages and data
library(limma); library(edgeR); library(systemPipeR)
targets <- read.delim("targetsPE.txt", comment.char="#")
counts <- read.delim("results/countDFeByg.xls", row.names=1, comment.char="#")
samples <- as.character(targets$Factor)
names(samples) <- paste(as.character(targets$SampleName), "", sep = "")
counts <- counts[, names(samples)]
counts[is.na(counts)] <- 0

## Create desing and contrast matrix
targets$Factor <- factor(targets$Factor, levels=unique(targets$Factor), ordered=TRUE) # Order factor as in given column
design <- model.matrix(~0 + factor(targets$Factor)) 
colnames(design) <- levels(targets$Factor)
cmp <- systemPipeR::readComp(file="targetsPE.txt", format="vector", delim="-")[[1]]
contrast.matrix <- makeContrasts(contrasts=cmp, levels=colnames(design))

## Run limma voom
dge <- DGEList(counts=counts, group=as.character(targets$Factor))
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, contrast.matrix) # Computes estimated coefficients and standard errors for a given set of contrasts.
fit2 <- eBayes(fit2)
all_comparisons <- sapply(cmp, function(x) topTable(fit2, coef=x, number=nrow(counts)), simplify=FALSE)
names(all_comparisons)
sapply(all_comparisons, dim)

## Function to assemble list in master DEG data.frame. Note it is identical to function under (A)
assembleFromList <- function(x=all_comparisons, counts) {
    masterDEG <- data.frame(row.names = rownames(counts))
    for(i in seq_along(x)) {
        singleDF <- x[[i]]
        fix_colnames <- colnames(singleDF)
        fix_colnames[fix_colnames %in% c("logFC", "adj.P.Val")] <- c("logFC", "FDR")
        fix_colnames <- paste(paste(names(x)[i], collapse = "-"), fix_colnames, sep = "_")
        colnames(singleDF) <- fix_colnames
        masterDEG <- cbind(masterDEG, singleDF[rownames(masterDEG),])
    }
    masterDEGdf <- na.omit(masterDEG)
    return(masterDEG)
}
masterDEGdf <- assembleFromList(x=all_comparisons, counts=counts) # Function defined under (A)
write.table(masterDEGdf, "./results/limmavoom_allcomp.xls", quote = FALSE, sep = "\t", col.names = NA)

## Plot DEG results
voomDF <- read.delim("results/limmavoom_allcomp.xls", row.names = 1, check.names = FALSE)
DEG_list <- filterDEGs(degDF = voomDF, filter = c(Fold = 2, FDR = 10))

