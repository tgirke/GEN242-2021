---
title: NGS Analysis Basics
author: "Author: Thomas Girke"
date: "Last update: 13 February, 2021" 
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
weith: 7
type: docs
---

<!--
- Compile from command-line
Rscript -e "rmarkdown::render('Rsequences.Rmd', c('html_document'), clean=F); knitr::knit('Rsequences.Rmd', tangle=TRUE)"; Rscript ../md2jekyll.R Rsequences.knit.md 10; Rscript -e "rmarkdown::render('Rsequences.Rmd', c('pdf_document'))"
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



# Overview

## Sequence Analysis in R and Bioconductor

__R Base__

* Some basic string handling utilities. Wide spectrum of numeric data analysis tools.

__Bioconductor__

Bioconductor packages provide much more sophisticated string handling utilities for sequence analysis [@Lawrence2013-kt, @Huber2015-ag].

* [Biostrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html): general sequence analysis environment
* [ShortRead](http://bioconductor.org/packages/release/bioc/html/ShortRead.html): pipeline for short read data
* [IRanges](http://bioconductor.org/packages/release/bioc/html/IRanges.html): low-level infrastructure for range data
* [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html): high-level infrastructure for range data
* [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html): managing transcript centric annotations
* [GenomicAlignments](http://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html): handling short genomic alignments
* [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html): interface to  `samtools`, `bcftools` and `tabix` 
* [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html): genome annotation data
* [biomaRt](http://bioconductor.org/packages/release/bioc/html/biomaRt.html): interface to BioMart annotations
* [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html): Annotation imports, interface to online genome browsers
* [HelloRanges](http://bioconductor.org/packages/release/bioc/html/HelloRanges.html): Bedtools semantics in Bioc's Ranges infrastructure


# Package Requirements

Several Bioconductor packages are required for this tutorial. To install them, execute
the following lines in the R console. Please also make sure that you have a recent R version
installed on your system. R versions `3.3.x` or higher are recommended.


```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings", "GenomicRanges", "GenomicRanges", "rtracklayer", "systemPipeR", "seqLogo", "ShortRead"))
```

# Strings in R Base

## Basic String Matching and Parsing

### String matching

Generate sample sequence data set


```r
myseq <- c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC")
```

String searching with regular expression support

```r
myseq[grep("ATG", myseq)] 
```

```
## [1] "ATGCAGACATAGTG" "ATGAACATAGATCC"
```

Searches `myseq` for first match of pattern "AT"

```r
pos1 <- regexpr("AT", myseq) 
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
```

```
## [1] 1 1 7
```

```
## [1] 2 2 2
```

Searches `myseq` for all matches of pattern "AT"

```r
pos2 <- gregexpr("AT", myseq) 
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Returns positions of matches in first sequence
```

```
## [1] 1 9
```

```
## [1] 2 2
```

String substitution with regular expression support

```r
gsub("^ATG", "atg", myseq) 
```

```
## [1] "atgCAGACATAGTG" "atgAACATAGATCC" "GTACAGATCAC"
```

### Positional parsing

```r
nchar(myseq) # Computes length of strings
```

```
## [1] 14 14 11
```

```r
substring(myseq[1], c(1,3), c(2,5)) # Positional parsing of several fragments from one string
```

```
## [1] "AT"  "GCA"
```

```r
substring(myseq, c(1,4,7), c(2,6,10)) # Positional parsing of many strings
```

```
## [1] "AT"   "AAC"  "ATCA"
```

## Random Sequence Generation

### Random DNA sequences of any length


```r
rand <- sapply(1:100, function(x) paste(sample(c("A","T","G","C"), sample(10:20), replace=T), collapse=""))
rand[1:3]
```

```
## [1] "CACGATCACGTCCGCTA" "AGTATTCTAGGCGATGG" "GAGAACTTCTCGAG"
```

### Count identical sequences


```r
table(c(rand[1:4], rand[1]))
```

```
## 
##    AGTATTCTAGGCGATGG    CACGATCACGTCCGCTA       GAGAACTTCTCGAG GGTCCAATAGTGTGACGACC 
##                    1                    2                    1                    1
```

### Extract reads from reference

Note: this requires `Biostrings` package.


```r
library(Biostrings)
ref <- DNAString(paste(sample(c("A","T","G","C"), 100000, replace=T), collapse=""))
randstart <- sample(1:(length(ref)-15), 1000)
randreads <- Views(ref, randstart, width=15)
rand_set <- DNAStringSet(randreads)
unlist(rand_set)
```

```
## 15000-letter DNAString object
## seq: TTCGTACGGCAGACGCCAGGGAACGGATCCCCCCGTAGGAAGGCGG...ACCTGCCGTAGCTGGCCTGAGACCTCCTCGAACGACTATGGCATTA
```

# Sequences in Bioconductor

## Important Data Objects of Biostrings

### `XString` for single sequence

* `DNAString`: for DNA
* `RNAString`: for RNA
* `AAString`: for amino acid 
* `BString`: for any string

### `XStringSet` for many sequences
        
* `DNAStringSet``: for DNA
* `RNAStringSet`: for RNA
* `AAStringSet`: for amino acid 
* `BStringSet`: for any string

### `QualityScaleXStringSet` for sequences with quality data

* `QualityScaledDNAStringSet`: for DNA
* `QualityScaledRNAStringSet`: for RNA
* `QualityScaledAAStringSet`: for amino acid 
* `QualityScaledBStringSet`: for any string

## Sequence Import and Export

Download the following sequences to your current working directory and then import them into R: 
[ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn](ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn)


```r
dir.create("data", showWarnings = FALSE)
# system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn")
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.ffn", "data/AE004437.ffn")
```

Import FASTA file with `readDNAStringSet`

```r
myseq <- readDNAStringSet("data/AE004437.ffn")
myseq[1:3]
```

```
## DNAStringSet object of length 3:
##     width seq                                                                   names               
## [1]  1206 ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTC...GTCGTCGTTGTTCGACGCTGGCGGAACCCATGA gi|12057215|gb|AE...
## [2]   666 ATGAGCATCATCGAACTCGAAGGCGTGGTCAAA...GTCAACCTCGTCGATGGGGTGTTACACACGTGA gi|12057215|gb|AE...
## [3]  1110 ATGGCGTGGCGGAACCTCGGGCGGAACCGCGTG...AACGATCCGCCCGTCGAGGCGCTCGGCGAATGA gi|12057215|gb|AE...
```

Subset sequences with regular expression on sequence name field

```r
sub <- myseq[grep("99.*", names(myseq))]
length(sub)
```

```
## [1] 170
```

Export subsetted sequences to FASTA file

```r
writeXStringSet(sub, file="./data/AE004437sub.ffn", width=80)
```

Now inspect exported sequence file `AE004437sub.ffn` in a text editor
	
    
## Working with `XString` Containers

The `XString` stores the different types of biosequences in dedicated containers

```r
library(Biostrings)
d <- DNAString("GCATAT-TAC")
d
```

```
## 10-letter DNAString object
## seq: GCATAT-TAC
```

```r
d[1:4]
```

```
## 4-letter DNAString object
## seq: GCAT
```

RNA sequences

```r
r <- RNAString("GCAUAU-UAC") 
r <- RNAString(d) # Converts d to RNAString object
r
```

```
## 10-letter RNAString object
## seq: GCAUAU-UAC
```
Protein sequences

```r
p <- AAString("HCWYHH")
p
```

```
## 6-letter AAString object
## seq: HCWYHH
```

Any type of character strings

```r
b <- BString("I store any set of characters. Other XString objects store only the IUPAC characters.")
b
```

```
## 85-letter BString object
## seq: I store any set of characters. Other XString objects store only the IUPAC characters.
```

## Working with `XStringSet` Containers

`XStringSet` containers allow to store many biosequences in one object

```r
dset <- DNAStringSet(c("GCATATTAC", "AATCGATCC", "GCATATTAC")) 
names(dset) <- c("seq1", "seq2", "seq3") # Assigns names
dset[1:2]
```

```
## DNAStringSet object of length 2:
##     width seq                                                                   names               
## [1]     9 GCATATTAC                                                             seq1
## [2]     9 AATCGATCC                                                             seq2
```

Important utilities for `XStringSet` containers

```r
width(dset) # Returns the length of each sequences
```

```
## [1] 9 9 9
```

```r
d <- dset[[1]] # The [[ subsetting operator returns a single entry as XString object
dset2 <- c(dset, dset) # Appends/concatenates two XStringSet objects
dsetchar <- as.character(dset) # Converts XStringSet to named vector 
dsetone <- unlist(dset) # Collapses many sequences to a single one stored in a DNAString container
```

Sequence subsetting by positions:

```r
DNAStringSet(dset, start=c(1,2,3), end=c(4,8,5)) 
```

```
## DNAStringSet object of length 3:
##     width seq                                                                   names               
## [1]     4 GCAT                                                                  seq1
## [2]     7 ATCGATC                                                               seq2
## [3]     3 ATA                                                                   seq3
```

## Multiple Alignment Class

The `XMultipleAlignment` class stores the different types of multiple sequence alignments:


```r
origMAlign <- readDNAMultipleAlignment(filepath = system.file("extdata",
              "msx2_mRNA.aln", package = "Biostrings"), format = "clustal")
origMAlign
```

```
## DNAMultipleAlignment with 8 rows and 2343 columns
##      aln                                                                        names               
## [1] -----TCCCGTCTCCGCAGCAAAAAAGTTTGAGTCG...TTGTCCAAACTCACAATTAAAAAAAAAAAAAAAAA gi|84452153|ref|N...
## [2] ------------------------------------...----------------------------------- gi|208431713|ref|...
## [3] ------------------------------------...----------------------------------- gi|118601823|ref|...
## [4] ----------------------AAAAGTTGGAGTCT...----------------------------------- gi|114326503|ref|...
## [5] ------------------------------------...----------------------------------- gi|119220589|ref|...
## [6] ------------------------------------...----------------------------------- gi|148540149|ref|...
## [7] --------------CGGCTCCGCAGCGCCTCACTCG...----------------------------------- gi|45383056|ref|N...
## [8] GGGGGAGACTTCAGAAGTTGTTGTCCTCTCCGCTGA...----------------------------------- gi|213515133|ref|...
```

## Basic Sequence Manipulations

### Reverse and Complement


```r
randset <- DNAStringSet(rand)
complement(randset[1:2])
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    17 GTGCTAGTGCAGGCGAT
## [2]    17 TCATAAGATCCGCTACC
```

```r
reverse(randset[1:2])
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    17 ATCGCCTGCACTAGCAC
## [2]    17 GGTAGCGGATCTTATGA
```

```r
reverseComplement(randset[1:2])
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    17 TAGCGGACGTGATCGTG
## [2]    17 CCATCGCCTAGAATACT
```

## Translate DNA into Protein

```r
translate(randset[1:2])
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[1]]':
## last 2 bases were ignored
```

```
## Warning in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet], : in 'x[[2]]':
## last 2 bases were ignored
```

```
## AAStringSet object of length 2:
##     width seq
## [1]     5 HDHVR
## [2]     5 SILGD
```

## Pattern Matching

### Pattern matching with mismatches

Find pattern matches in reference 

```r
myseq1 <- readDNAStringSet("./data/AE004437.ffn") 
mypos <- matchPattern("ATGGTG", myseq1[[1]], max.mismatch=1) 
```

Count only the corresponding matches

```r
countPattern("ATGGCT", myseq1[[1]], max.mismatch=1) 
```

```
## [1] 3
```

Count matches in many sequences

```r
vcountPattern("ATGGCT", myseq1, max.mismatch=1)[1:20]
```

```
##  [1] 3 0 5 4 1 2 2 1 4 3 0 0 1 2 0 1 4 0 0 1
```

Results shown in DNAStringSet object

```r
tmp <- c(DNAStringSet("ATGGTG"), DNAStringSet(mypos)) 
```

Return a consensus  matrix for query and hits

```r
consensusMatrix(tmp)[1:4,] 
```

```
##   [,1] [,2] [,3] [,4] [,5] [,6]
## A    3    0    0    0    0    0
## C    1    1    0    0    0    0
## G    0    0    4    4    1    4
## T    0    3    0    0    3    0
```

Find all pattern matches in reference

```r
myvpos <- vmatchPattern("ATGGCT", myseq1, max.mismatch=1) 
myvpos # The results are stored as MIndex object.
```

```
## MIndex object of length 2058
## $`gi|12057215|gb|AE004437.1|:248-1453 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 3 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         6         6
##   [2]       383       388         6
##   [3]       928       933         6
## 
## $`gi|12057215|gb|AE004437.1|:1450-2115 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 0 ranges and 0 metadata columns:
##        start       end     width
##    <integer> <integer> <integer>
## 
## $`gi|12057215|gb|AE004437.1|:2145-3254 Halobacterium sp. NRC-1, complete genome`
## IRanges object with 5 ranges and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]         1         6         6
##   [2]        94        99         6
##   [3]       221       226         6
##   [4]       535       540         6
##   [5]       601       606         6
## 
## ...
## <2055 more elements>
```

```r
Views(myseq1[[1]], start(myvpos[[1]]), end(myvpos[[1]])) # Retrieves the result for single entry
```

```
## Views on a 1206-letter DNAString subject
## subject: ATGACTCGGCGGTCTCGTGTCGGTGCCGGCCTCGCAGCCATTGT...TTGCGATCGTCGTCGTCGTTGTTCGACGCTGGCGGAACCCATGA
## views:
##       start end width
##   [1]     1   6     6 [ATGACT]
##   [2]   383 388     6 [ATGGCA]
##   [3]   928 933     6 [ATGACT]
```

Return all matches

```r
sapply(seq(along=myseq1), function(x) 
       as.character(Views(myseq1[[x]], start(myvpos[[x]]), end(myvpos[[x]]))))[1:4] 
```

### Pattern matching with regular expression support


```r
myseq <- DNAStringSet(c("ATGCAGACATAGTG", "ATGAACATAGATCC", "GTACAGATCAC"))
myseq[grep("^ATG", myseq, perl=TRUE)] # String searching with regular expression support
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    14 ATGCAGACATAGTG
## [2]    14 ATGAACATAGATCC
```

```r
pos1 <- regexpr("AT", myseq) # Searches 'myseq' for first match of pattern "AT"
as.numeric(pos1); attributes(pos1)$match.length # Returns position information of matches
```

```
## [1] 1 1 7
```

```
## [1] 2 2 2
```

```r
pos2 <- gregexpr("AT", myseq) # Searches 'myseq' for all matches of pattern "AT"
as.numeric(pos2[[1]]); attributes(pos2[[1]])$match.length # Match positions in first sequence
```

```
## [1] 1 9
```

```
## [1] 2 2
```

```r
DNAStringSet(gsub("^ATG", "NNN", myseq)) # String substitution with regular expression support
```

```
## DNAStringSet object of length 3:
##     width seq
## [1]    14 NNNCAGACATAGTG
## [2]    14 NNNAACATAGATCC
## [3]    11 GTACAGATCAC
```

## PWM Viewing and Searching

### Plot with `seqLogo`


```r
library(seqLogo) 
```

```
## Loading required package: grid
```

```r
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
pwm
```

```
##        [,1]      [,2]      [,3]
## A 0.0000000 0.0000000 0.2312611
## C 0.0000000 0.3157205 0.0000000
## G 0.3685591 0.2312611 0.0000000
## T 0.0000000 0.0000000 0.3157205
```

```r
seqLogo(t(t(pwm) * 1/colSums(pwm)))
```

<img src="Rsequences_files/figure-html/pwm_logo-1.png" width="672" />

### Plot with `ggseqlogo`

The `ggseqlogo` package ([manual](https://omarwagih.github.io/ggseqlogo/))
provides many customization options for plotting sequence logos. It also supports
various alphabets including sequence logos for amino acid sequences.



```r
library(ggplot2); library(ggseqlogo)
pwm <- PWM(DNAStringSet(c("GCT", "GGT", "GCA"))) 
ggseqlogo(pwm)
```

<img src="Rsequences_files/figure-html/pwm_logo2-1.png" width="672" />

Search sequence for PWM matches with score better than `min.score`

```r
chr <- DNAString("AAAGCTAAAGGTAAAGCAAAA") 
matchPWM(pwm, chr, min.score=0.9) 
```

```
## Views on a 21-letter DNAString subject
## subject: AAAGCTAAAGGTAAAGCAAAA
## views:
##       start end width
##   [1]     4   6     3 [GCT]
##   [2]    10  12     3 [GGT]
##   [3]    16  18     3 [GCA]
```

# NGS Sequences

## Sequence and Quality Data: FASTQ Format

Four lines per sequence:

1. ID
2. Sequence
3. ID
4. Base call qualities (Phred scores) as ASCII characters

The following gives an example of 3 Illumina reads in a FASTQ file. The numbers
at the beginning of each line are not part of the FASTQ format. They have been added 
solely for illustration purposes.

```
1. @SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
2. CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
3. +SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
4. BA@7>B=>:>>7@7@>>9=BAA?;>52;>:9=8.=A
1. @SRR038845.41 HWI-EAS038:6:1:0:1474 length=36
2. CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
3. +SRR038845.41 HWI-EAS038:6:1:0:1474 length=36
4. BCCBA@BB@BBBBAB@B9B@=BABA@A:@693:@B=
1. @SRR038845.53 HWI-EAS038:6:1:1:360 length=36
2. GTTCAAAAAGAACTAAATTGTGTCAATAGAAAACTC
3. +SRR038845.53 HWI-EAS038:6:1:1:360 length=36
4. BBCBBBBBB@@BAB?BBBBCBC>BBBAA8>BBBAA@
```

## Sequence and Quality Data: `QualityScaleXStringSet`

Phred quality scores are integers from 0-50 that are
stored as ASCII characters after adding 33. The basic R functions `rawToChar` and
`charToRaw` can be used to interconvert among their representations.

Phred score interconversion

```r
phred <- 1:9
phreda <- paste(sapply(as.raw((phred)+33), rawToChar), collapse="")
phreda
```

```
## [1] "\"#$%&'()*"
```

```r
as.integer(charToRaw(phreda))-33 
```

```
## [1] 1 2 3 4 5 6 7 8 9
```

Construct `QualityScaledDNAStringSet` from scratch

```r
dset <- DNAStringSet(sapply(1:100, function(x) paste(sample(c("A","T","G","C"), 20, replace=T), collapse=""))) # Creates random sample sequence.
myqlist <- lapply(1:100, function(x) sample(1:40, 20, replace=T)) # Creates random Phred score list.
myqual <- sapply(myqlist, function(x) toString(PhredQuality(x))) # Converts integer scores into ASCII characters.
myqual <- PhredQuality(myqual) # Converts to a PhredQuality object.
dsetq1 <- QualityScaledDNAStringSet(dset, myqual) # Combines DNAStringSet and quality data in QualityScaledDNAStringSet object.
dsetq1[1:2]
```

```
##   A QualityScaledDNAStringSet instance containing:
## 
## DNAStringSet object of length 2:
##     width seq
## [1]    20 ACCAGCAGGGGTGTCGGTCA
## [2]    20 GGACCTGCACAGAGCAGCGT
## 
## PhredQuality object of length 2:
##     width seq
## [1]    20 G;7H5GA7-(B6@0/GB'67
## [2]    20 6D:31,,<A&%3%:9H2I'E
```

## Processing FASTQ Files with ShortRead

The following expains the basic usage of `ShortReadQ` objects. To make the sample code work, 
download and unzip this [file](http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_6_10_2012/Rsequences/data.zip) to your current working directory.
The following code performs the download for you.


```r
library(ShortRead)
download.file("http://faculty.ucr.edu/~tgirke/HTML_Presentations/Manuals/Workshop_Dec_6_10_2012/Rsequences/data.zip", "data.zip")
unzip("data.zip")
```

Important utilities for accessing FASTQ files

```r
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 
(fq <- readFastq(fastq[1])) # Imports first FASTQ file
```

```
## class: ShortReadQ
## length: 1000 reads; width: 36 cycles
```

```r
countLines(dirPath="./data", pattern=".fastq$")/4 # Counts numbers of reads in FASTQ files
```

```
## SRR038845.fastq SRR038846.fastq SRR038848.fastq SRR038850.fastq 
##            1000            1000            1000            1000
```

```r
id(fq)[1] # Returns ID field
```

```
## BStringSet object of length 1:
##     width seq
## [1]    43 SRR038845.3 HWI-EAS038:6:1:0:1938 length=36
```

```r
sread(fq)[1] # Returns sequence
```

```
## DNAStringSet object of length 1:
##     width seq
## [1]    36 CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
```

```r
quality(fq)[1] # Returns Phred scores 
```

```
## class: FastqQuality
## quality:
## BStringSet object of length 1:
##     width seq
## [1]    36 BA@7>B=>:>>7@7@>>9=BAA?;>52;>:9=8.=A
```

```r
as(quality(fq), "matrix")[1:4,1:12] # Coerces Phred scores to numeric matrix
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
## [1,]   33   32   31   22   29   33   28   29   25    29    29    22
## [2,]   33   34   34   33   32   31   33   33   31    33    33    33
## [3,]   33   33   34   33   33   33   33   33   33    31    31    33
## [4,]   33   33   33   33   31   33   28   31   28    32    33    33
```

```r
ShortReadQ(sread=sread(fq), quality=quality(fq), id=id(fq)) # Constructs a ShortReadQ from components
```

```
## class: ShortReadQ
## length: 1000 reads; width: 36 cycles
```

## FASTQ Quality Reports 

### Using `systemPipeR`

The following `seeFastq` and `seeFastqPlot` functions generate and plot a series of useful quality statistics for a set of FASTQ files.


```r
library(systemPipeR)
fqlist <- seeFastq(fastq=fastq, batchsize=800, klength=8) # For real data set batchsize to at least 10^5 
seeFastqPlot(fqlist)
```

<img src="Rsequences_files/figure-html/see_fastq-1.png" width="1344" />

Handles many samples in one PDF file. For more details see [here](http://bioconductor.org/packages/devel/bioc/vignettes/systemPipeR/inst/doc/systemPipeR.html) 

### Using `ShortRead`

The `ShortRead` package contains several FASTQ quality reporting functions.

```r
sp <- SolexaPath(system.file('extdata', package='ShortRead'))
fl <- file.path(analysisPath(sp), "s_1_sequence.txt") 
fls <- c(fl, fl) 
coll <- QACollate(QAFastqSource(fls), QAReadQuality(), QAAdapterContamination(), 
	    QANucleotideUse(), QAQualityUse(), QASequenceUse(), QAFrequentSequence(n=10), 
		QANucleotideByCycle(), QAQualityByCycle())
x <- qa2(coll, verbose=TRUE)
res <- report(x)
if(interactive())
browseURL(res) 
```

## Filtering and Trimming FASTQ Files with ShortRead 

### Adaptor trimming


```r
fqtrim <- trimLRPatterns(Rpattern="GCCCGGGTAA", subject=fq)
sread(fq)[1:2] # Before trimming
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    36 CAACGAGTTCACACCTTGGCCGACAGGCCCGGGTAA
## [2]    36 CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
```

```r
sread(fqtrim)[1:2] # After trimming
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]    26 CAACGAGTTCACACCTTGGCCGACAG
## [2]    36 CCAATGATTTTTTTCCGTGTTTCAGAATACGGTTAA
```
### Read counting and duplicate removal


```r
tables(fq)$distribution # Counts read occurences
```

```
##   nOccurrences nReads
## 1            1    948
## 2            2     26
```

```r
sum(srduplicated(fq)) # Identifies duplicated reads
```

```
## [1] 26
```

```r
fq[!srduplicated(fq)]
```

```
## class: ShortReadQ
## length: 974 reads; width: 36 cycles
```

### Trimming low quality tails


```r
cutoff <- 30
cutoff <- rawToChar(as.raw(cutoff+33))
sread(trimTails(fq, k=2, a=cutoff, successive=FALSE))[1:2]
```

```
## DNAStringSet object of length 2:
##     width seq
## [1]     4 CAAC
## [2]    20 CCAATGATTTTTTTCCGTGT
```

### Removal of reads with Phred scores below a threshold value


```r
cutoff <- 30
qcount <- rowSums(as(quality(fq), "matrix") <= 20) 
fq[qcount == 0] # Number of reads where all Phred scores >= 20
```

```
## class: ShortReadQ
## length: 349 reads; width: 36 cycles
```

### Removal of reads with x Ns and/or low complexity segments


```r
filter1 <- nFilter(threshold=1) # Keeps only reads without Ns
filter2 <- polynFilter(threshold=20, nuc=c("A","T","G","C")) # Removes reads with >=20 of one nucleotide
filter <- compose(filter1, filter2)
fq[filter(fq)]
```

```
## class: ShortReadQ
## length: 989 reads; width: 36 cycles
```

## Memory Efficient FASTQ Processing

Streaming through FASTQ files with `FastqStreamer` and random sampling reads with `FastqSampler`


```r
fq <- yield(FastqStreamer(fastq[1], 50)) # Imports first 50 reads 
fq <- yield(FastqSampler(fastq[1], 50)) # Random samples 50 reads 
```

Streaming through a FASTQ file while applying filtering/trimming functions and writing the results to a new file
 here `SRR038845.fastq_sub` in `data` directory.


```r
f <- FastqStreamer(fastq[1], 50) 
while(length(fq <- yield(f))) {
	fqsub <- fq[grepl("^TT", sread(fq))] 
	writeFastq(fqsub, paste(fastq[1], "sub", sep="_"), mode="a", compress=FALSE)
}
close(f)
```

# Range Operations  

## Important Data Objects for Range Operations

* `IRanges`: stores range data only (IRanges library)
* `GRanges`: stores ranges and annotations (GenomicRanges library)
* `GRangesList`: list version of GRanges container (GenomicRanges library)

## Range Data Are Stored in `IRanges` and `GRanges` Containers

### Construct `GRanges` Object 


```r
library(GenomicRanges); library(rtracklayer)
gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)), score = 1:10, GC = seq(1, 0, length = 10)) # Example of creating a GRanges object with its constructor function.
```

### Import GFF into `GRanges` Object

```r
gff <- import.gff("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/gff3.gff") # Imports a simplified GFF3 genome annotation file.
seqlengths(gff) <- end(ranges(gff[which(values(gff)[,"type"]=="chromosome"),])) 
names(gff) <- 1:length(gff) # Assigns names to corresponding slot
gff[1:4,]
```

```
## GRanges object with 4 ranges and 10 metadata columns:
##     seqnames     ranges strand |   source       type     score     phase                  ID
##        <Rle>  <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer>         <character>
##   1     Chr1 1-30427671      + |   TAIR10 chromosome        NA      <NA>                Chr1
##   2     Chr1  3631-5899      + |   TAIR10 gene              NA      <NA>           AT1G01010
##   3     Chr1  3631-5899      + |   TAIR10 mRNA              NA      <NA>         AT1G01010.1
##   4     Chr1  3760-5630      + |   TAIR10 protein           NA      <NA> AT1G01010.1-Protein
##            Name                Note          Parent       Index Derives_from
##     <character>     <CharacterList> <CharacterList> <character>  <character>
##   1        Chr1                                            <NA>         <NA>
##   2   AT1G01010 protein_coding_gene                        <NA>         <NA>
##   3 AT1G01010.1                           AT1G01010           1         <NA>
##   4 AT1G01010.1                                            <NA>  AT1G01010.1
##   -------
##   seqinfo: 7 sequences from an unspecified genome
```

### Coerce `GRanges` object to `data.frame`

```r
as.data.frame(gff)[1:4, 1:7]
```

```
##   seqnames start      end    width strand source       type
## 1     Chr1     1 30427671 30427671      + TAIR10 chromosome
## 2     Chr1  3631     5899     2269      + TAIR10       gene
## 3     Chr1  3631     5899     2269      + TAIR10       mRNA
## 4     Chr1  3760     5630     1871      + TAIR10    protein
```

### Coerce `GRanges` to `RangedData` object and vice versa
















































