#!/bin/bash -l

######################################################
## HW02: Linux Basics Using Bioinformatics Examples ##
######################################################
## Author: First Last Name
## Last update: 06-Apr-2021

## Download Halobacterium proteome and inspect it
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Escherichia_coli/reference/GCA_000008865.2_ASM886v2/GCA_000008865.2_ASM886v2_protein.faa.gz 
gunzip GCA_000008865.2_ASM886v2_protein.faa.gz
mv GCA_000008865.2_ASM886v2_protein.faa ecoli.faa

# less ecoli.faa # press q to quit

## How many protein sequences are stored in the downloaded file?
grep '>' ecoli.faa | wc
grep '^>' ecoli.faa --count
# Answer: 5,153 protein sequences

## How many proteins contain the pattern "WxHxxH" or "WxHxxHH"?
egrep 'W.H..H{1,2}' ecoli.faa --count
# Answer: 15 matches

## Use less to find IDs for pattern matches or use awk
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' ecoli.faa | less
awk --posix -v RS='>' '/W.H..(H){1,2}/ { print ">" $0;}' ecoli.faa | grep '^>' | cut -c 2- | cut -f 1 -d\ > myIDs

## Create a BLASTable database with formatdb
module load ncbi-blast
makeblastdb -in ecoli.faa -out ecoli.faa -dbtype prot -hash_index -parse_seqids

## Query BLASTable database by IDs stored in a file (e.g. myIDs)
blastdbcmd -db ecoli.faa -dbtype prot -entry_batch myIDs -get_dups -out myseq.fasta

## Run BLAST search for sequences stored in myseq.fasta
blastp -query myseq.fasta -db ecoli.faa -outfmt 0 -evalue 1e-6 -out blastp.out
blastp -query myseq.fasta -db ecoli.faa -outfmt 6 -evalue 1e-6 -out blastp.tab
# Answer created file with -outfmt 6 format and uploaded to hw1 repos on GitHub classroom

## Return system time and host name
date 
hostname 
# Answer: Apr 06 16:25:58 PDT 2021
# Answer: parrot or some other host name of a node on HPCC

