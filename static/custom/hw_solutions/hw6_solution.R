###########################################
## Solutions for HW6 - Sequence Analysis ##
###########################################


## Environment and Input Data
## Note; if not available, please install R.utils package for unzip() step
library(ShortRead)
download.file("http://cluster.hpcc.ucr.edu/~tgirke/HTML_Presentations/Manuals/testdata/samplefastq/data.zip", "data.zip")
unzip("data.zip")
fastq <- list.files("data", "*.fastq$"); fastq <- paste("data/", fastq, sep="")
names(fastq) <- paste("flowcell6_lane", 1:length(fastq), sep="_") 


## A. Demultiplexing 
library(ShortRead); library(Biostrings); library(GenomicRanges); 
demultiplex <- function(x, barcode, nreads, cutoff) {
	f <- FastqStreamer(x, nreads) 
	while(length(fq <- yield(f))) {
		for(i in barcode) {
			pattern <- paste("^", i, sep="")
			fqsub <- fq[grepl(pattern, sread(fq))] 
			if(length(fqsub) > 0) {
			    fqsub <- trimTails(fqsub, k=2, a=rawToChar(as.raw(cutoff+33)), successive=FALSE)
				writeFastq(fqsub, paste(x, i, sep="_"), mode="a", compress=FALSE)
			}
		}
	}
	close(f)
}
demultiplex(x=fastq[1], barcode=c("TT", "AA", "GG"), nreads=50, cutoff=20)

## B. Sequence Parsing 

## Task 1
library(Biostrings); library(rtracklayer)
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.gff", "data/AE004437.gff")
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Halobacterium_sp_uid217/AE004437.fna", "data/AE004437.fna")
chr <- readDNAStringSet("data/AE004437.fna")
gff <- import("data/AE004437.gff")
gffgene <- gff[values(gff)[,"type"]=="gene"]
gene <- DNAStringSet(Views(chr[[1]], IRanges(start(gffgene), end(gffgene))))
names(gene) <- values(gffgene)[,"locus_tag"]
pos <- values(gffgene[strand(gffgene) == "+"])[,"locus_tag"]
p1 <- translate(gene[names(gene) %in% pos])
names(p1) <- names(gene[names(gene) %in% pos])
neg <- values(gffgene[strand(gffgene) == "-"])[,"locus_tag"]
p2 <- translate(reverseComplement(gene[names(gene) %in% neg]))
names(p2) <- names(gene[names(gene) %in% neg])
writeXStringSet(c(p1, p2), "mypep.fasta")

## Task 2
reduced_ranges <- reduce(gffgene) # Collapses overlapping ranges to single ranges.
DNAStringSet(Views(chr[[1]], IRanges(start(reduced_ranges), end(reduced_ranges)))) 

## Task 3
intergenic <- gaps(reduced_ranges) # Returns uncovered regions
DNAStringSet(Views(chr[[1]], IRanges(start(intergenic), end(intergenic))))
