##########################################################
### Identify Over-Represented Strings in Sequence Sets ###
##########################################################
## Authors: Ugoeze Nwokedi, Kevin Horan, Thomas Girke (thomas.girke@ucr.edu, UC Riverside)
## Last update: July 14, 2007
## Utility: 
##	(1) Generate a database of all possible DNA/RNA/PEP strings of a given length
##	(2) Sequence import into R 
##	(3) Determine the frequency of the strings in reference and test sequences 
##	(4) Calculate over-representation p-value
##	(5) Large-scale meta analysis by looping over entries in string database

## How to run the script:
## 	source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/wordFinder.R")

##############################################################################################
## (1) Function to Create a Datatabes of All Possible DNA/RNA/PEP Strings of a Given Length ##
##############################################################################################
## Create X'mer Database
makeXmerDB <- function(type="DNA", xmer=4) {
	## Define sequence type
	if(type=="DNA") {
		N <- c("A", "T", "G", "C")
	} 
	if(type=="RNA") {
		N <- c("A", "U", "G", "C")
	}
	if(type=="PEP") {
		N <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	}
	
	## Create sequence strings
	myseq <- N
	seqlength <- xmer
	if(xmer>=2) {
		for(i in 1:(seqlength-1)) {
			myseq <- sapply(as.vector(myseq), function(x) paste(x, N, sep=""))
		}
	}
	myseq <- as.vector(myseq)
	return(myseq)
}
cat("\n# (1.1) Execute the following command, to generate a 6-mer database:\n\t myDB <- makeXmerDB(type=\"DNA\", xmer=6)\n")

## Remove reverse and complement strings in pattern database. They create duplicated results when searching both DNA strands.
## Note that the number of removed sequences will be less than half of the number of elements in myDB, because the reverse and 
## complement of palindromic sequences is identical to their original sequences.
removeDups <- function(myDBF=myDB) {
	myDB <- myDB
	temp <- paste(myDB, as.vector(sapply(myDB, rev_comp_fct, rev=T, comp=T)), sep="_")
	temp <- strsplit(temp, "_")
	temp <- lapply(temp, sort)
	temp <- unlist(lapply(temp, paste, collapse="_"))
	myDB <- myDB[!duplicated(temp)]
	return(myDB)
}
cat("\n# (1.2) Only for searches of both DNA strands: remove reverse and complement patterns which create duplicated results: \n\t myDB <- removeDups(myDB=myDB) \n")

##################################
## (2) Sequence Import Function ##
##################################
seq_imp_fct <- function(fileloc) {
	my_fasta <- readLines(fileloc) # reads file line-wise into vector 
	y <- regexpr("^[^>]", my_fasta, perl=T) # identifies all fields that do not start with a '>' sign 
	y <- as.vector(y);  y[y==-1] <- 0
	index <- which(y==0)
	distance <- data.frame(start=index[1:(length(index)-1)], end=index[2:length(index)])
	distance <- rbind(distance, c(distance[length(distance[,1]),2], length(y)+1)) # gets data for last entry 
	distance <- data.frame(distance, dist=distance[,2]-distance[,1])
	seq_no <- 1:length(y[y==0])
	index <- rep(seq_no, as.vector(distance[,3]))
	my_fasta <- data.frame(index, y, my_fasta)
	my_fasta[my_fasta[,2]==0,1] <- 0
	seq <- tapply(as.vector(my_fasta[,3]), factor(my_fasta[,1]), paste, collapse="", simplify=F)
	seq <- as.vector(seq[2:length(seq)])
	Desc <- as.vector(my_fasta[c(grep("^>", as.character(my_fasta[,3]), perl = TRUE)),3])
	ID <- gsub("^>| .*", "", as.character(Desc), perl=T)
	Desc <- gsub("^.*? ", "", as.character(Desc), perl=T)
	my_fasta <- data.frame(ID, Desc, Length=nchar(seq), seq)
	my_fasta
}	
cat("\n# (2.1) Import sample sequences into R like this: \n\t myseq <- seq_imp_fct(\"http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/Samples/sampleDNA.txt\"); myseq <- myseq[1:8,]\n")
cat("\n# (2.2) Run this command to set all sequences to upper case since the function gregexpr() is case sensitive! \n\t myseq$seq <- gsub(\"([A-Z])\", \"\\\\U\\\\1\", as.character(myseq$seq), perl=T, ignore.case=T)\n")
cat("\n# (2.3) To speed up the sequence loading process, save and load them like this: \n\t save(myseq, file=\"myseq\"); load(file=\"myseq\") \n")   

## Reverse and/or Complement Function
rev_comp_fct <- function(seq=mystr, rev=T, comp=T) {
	if(rev==T) {
		seq <- as.vector(unlist(strsplit(seq, split="")))
		seq <- rev(seq)
		seq <- paste(seq, collapse="")
	}	
	if(comp==T) {
		seq <- gsub("A", "1", seq, ignore.case = T)
		seq <- gsub("T", "2", seq, ignore.case = T)
		seq <- gsub("C", "3", seq, ignore.case = T)
		seq <- gsub("G", "4", seq, ignore.case = T)		
		seq <- gsub("1", "T", seq, ignore.case = T)
		seq <- gsub("2", "A", seq, ignore.case = T)
		seq <- gsub("3", "G", seq, ignore.case = T)
		seq <- gsub("4", "C", seq, ignore.case = T)
	}
	seq
}

#################################
## (3) Match Counting Function ##
#################################
## The matchCounter function can be used to calculate how often every word in the counts 

## Define matchCounter function
matchCounter <- function(seq, ids, myDB, strands=2, pos=T, mm=0) {
	testseq <- seq
	myids <- ids
	
	## Search only positive strand
	if(strands==1) {
		mymatch <- lapply(myDB, function(x) gregexpr(x, testseq, perl=T))
	}

	## Search both strands
	if(strands==2) {
		myDB2 <- paste(myDB, "|", rev_comp_fct(myDB, rev=T, comp=T), sep="")
		mymatch <- lapply(myDB2, function(x) gregexpr(x, testseq, perl=T))
	}
	names(mymatch) <- myDB
	motifs <- myDB[sort(rep(1:length(myDB), length(myids)))] 
	ids <- myids[rep(1:length(myids), length(myDB))]
	count <- unlist(lapply(names(mymatch), function(x) lapply(mymatch[[x]], function(y)  sum(y!=-1) )))
	if(pos==F) {
		mymatchDF <- data.frame(motifs, ids, count)
		return(mymatchDF)
	}
	if(pos==T) {
		pos <- unlist(lapply(names(mymatch), function(x) lapply(mymatch[[x]], function(y)  paste(y, collapse=", ") )))
		pos[pos=="-1"] <- "0"
		mymatchDF <- data.frame(motifs, ids, count, pos)
		if(mm>0) {
			mm <- lapply(myDB, function(x) agrep(x, testseq, max=mm, value=F))
			fixzero <- sapply(mm, length)
			mm[fixzero==0] <- 0
			mm <- lapply(mm, function(x) 1:length(testseq) %in% x)
			names(mm) <- myDB
			mm <- unlist(mm)
			mymatchDF <- data.frame(motifs, ids, count, pos, mm)
		}
		return(mymatchDF)
	}
}
cat("\n# (3.1) Run matchCounter function: \n\t mymatchDF <- matchCounter(seq=as.character(myseq$seq), ids=as.character(myseq$ID), myDB=myDB, strands=2, pos=T, mm=0); mymatchDF[1:12,] \n \t\t # To speed up this process, set \"pos=F\".\n")


################################################
## (4) Calculate Over-Representation P-Values ##
################################################
## Formula and variables:
	## phyper(k, D, n-D, N, lower.tail=F)
		## k: number of motif matches in sample sequence (number of good items in sequence sample)  
		## D: number of motif matches in sequence database (number of good items in sequence database)
		## n: number of words in sequence database minus D (number of bad items in sequence database)
		## N: number of words in sample sequence (number of drawn items) 

enrichTest <- function(myseq=myseq, mymatchDF=mymatchDF, xmer=6) {
	
	## Obtain values for required variables
	motifMatches <- tapply(mymatchDF$count, mymatchDF$motif, sum)
	k <- as.vector(mymatchDF$count)
	motifMatches <- tapply(mymatchDF$count, mymatchDF$motif, sum)
	D <- as.vector(motifMatches[mymatchDF$motifs])
	n <- sum(nchar(myseq$seq)-(xmer-1))
	N <- nchar(myseq$seq)-(xmer-1)
	names(N) <- myseq$ID
	N <- as.vector(N[mymatchDF$ids])
	
	## Calculate hypergeometric p-values
	pval <- phyper(k-1, D, n-D, N, lower.tail=F)
	
	## Bonferroni p-value adjustment. Multiplyer is the number of sequences that show at least one match for a given motif.
	N_Match_Seq <- tapply(mymatchDF$count, mymatchDF$motif, function(x) sum(x>=1))
	N_Match_Seq <- as.vector(N_Match_Seq[mymatchDF$motifs])
	N_Match_Seq[N_Match_Seq==0] <- 1
	adj_pval <- pval * N_Match_Seq
	adj_pval[adj_pval>1] <- 1

	## Assemble results in data frame
	mymatchDF <- data.frame(mymatchDF, D=D, n=n, N=N, pval=pval, adj_pval)
	return(mymatchDF)
}
cat("\n# (4.1) Calculate over-representation P-values: \n\t mymatchDF <- enrichTest(myseq=myseq, mymatchDF=mymatchDF, xmer=6); mymatchDF[order(mymatchDF$pval),][1:12,]\n")

##############################################################################
## (5) Large-Scale Meta Analysis by Looping Over Entries in String Database ##
##############################################################################
## This meta function is for genome-wide analyses. It runs the above steps in a memory 
## efficient manner by looping over the entries in the string database and storing 
## only the results for enriched patterns.

metaAnalysis <- function(seq=as.character(myseq$seq), myDB=myDB, cutoff=0.01) {
	cont <- NULL
	for(i in 1:length(myDB)) {
		mymatchDF <- matchCounter(seq=as.character(myseq$seq), ids=as.character(myseq$ID), myDB=myDB[i], strands=2, pos=T, mm=0)
		mymatchDF <- enrichTest(myseq=myseq, mymatchDF=mymatchDF, xmer=6)
		mymatchDF_filter <- mymatchDF[mymatchDF$adj_pval <= cutoff, ]
		cont <- rbind(cont, mymatchDF_filter)
		cat("Iteration: ", i, ", N rows: ", length(cont[,1]), "\n", sep="")
	}
	return(cont)
}
cat("\n# (5.1) Run Genome-Wide Analysis: \n\t finalDF <- metaAnalysis(seq=as.character(myseq$seq), myDB=myDB, cutoff=0.05) \n\n")

## How to run meta analysis:
	## Generate motif database
	## Import genome sequences: ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR7_blastsets/TAIR7_upstream_1000_20070405
	## Run metaAnalysis function





