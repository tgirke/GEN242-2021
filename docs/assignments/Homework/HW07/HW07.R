## A. Unstranded and strand-specific read counting for features of interest

## Task 1: Rerun the RNA-Seq workflow
## with the toy data sets up to the read quantification step.  In the read
## quantification step with `summarizeOverlaps` generate count tables for exons by
## genes (`eByg`) of the following three strand modes:

##   1. Unstranded 
##   2. Strand-specific for positive (sense) strand
##   3. Strand-specific for negative (antisense) strand
   
## The solution for generating the unstranded read counts is given below. Note,
## the upstream steps 1-4 in the RNA-Seq workflow only need to be rerun to
## generate the proper inputs for the read counting. Thus, they are not
## required to be included in the homework results (see `HW7.R` below).

summarizeOverlaps(eByg, bfl, mode="Union", 
                    ignore.strand=TRUE, 
                    # preprocess.reads=invertStrand,
                    inter.feature=FALSE, 
                    singleEnd=TRUE)

## Before attempting to solve this homework task please read the vignette
## Counting reads with `summarizeOverlaps`
## from the `GenomicAlignments` package that defines the `summarizeOverlap`
## function.

## Task 2: Provide R code that demonstrates that the two strand-specific count
## tables sum up to the values of the unstranded count table. 

## Task 3: Explain the utility (biological relevance) of the different strand
## counting modes used under Task 1. Include your explanation as comment text
## in your homework script (see `HW7.R` below). 

## B. Read counting for different feature types
## Task 4: Compute strand-specific count tables for the positive (sense) strand
## of the following feature types. The help files of `?exonsBy` and `?transcripts'
## provide useful information to solve these tasks.

##   1. Genes
##   2. Exons
##   3. Exons by genes 
##   4. Introns by transcripts
##   5. 5'-UTRs by transcripts

## C. DEG analysis

## Task 5: Perform the DEG analysis with `edgeR` as outlined under section 6 of
## the RNA-Seq workflow
## Use in one case for the DEG analysis the unstranded count table as input
## (from Task 1.1) and in another the sense strand count table (from Task 1.2).
## Compare the DEG result of the two methods in two separate 4-way Venn
## diagrams for the same sample comparisons used in the workflow example.

##   1. 4-way Venn diagram for unstranded count table
##   2. 4-way Venn diagram for sense strand count table


