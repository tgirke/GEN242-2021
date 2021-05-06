## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()
options(width=100, max.print=1000)
knitr::opts_chunk$set(
    eval=as.logical(Sys.getenv("KNITR_EVAL", "TRUE")),
    cache=as.logical(Sys.getenv("KNITR_CACHE", "TRUE")))


## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE------------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(fgsea)
})


## ----godb, eval=FALSE, warning=FALSE, message=FALSE---------------------------
## ## Load GOstats library
## library(GOstats); library(GO.db)
## ## Print complete GO term information for "GO:0003700"
## GOTERM$"GO:0003700"
## ## Print parent and children terms for a GO ID
## GOMFPARENTS$"GO:0003700"; GOMFCHILDREN$"GO:0003700"
## ## Print complete lineages of parents and children for a GO ID
## GOMFANCESTOR$"GO:0003700"; GOMFOFFSPRING$"GO:0003700"
## ## Print number of GO terms in each of the 3 ontologies
## zz <- eapply(GOTERM, function(x) x@Ontology); table(unlist(zz))
## ## Gene to GO mappings for an organism (here Arabidopsis)
## library(org.At.tair.db) # For human use org.Hs.eg.db
## xx <- as.list(org.At.tairGO2ALLTAIRS)


## ----keggdb, eval=FALSE, warning=FALSE, message=FALSE-------------------------
## ## KEGG Pathway ID to Hs Entrez ID list
## load_keggList <- function(org="ath") {
##     suppressMessages(suppressWarnings(library(KEGG.db)))
##     kegg_gene_list <- as.list(KEGGPATHID2EXTID) # All organisms in kegg
##     kegg_gene_list <- kegg_gene_list[grepl(org, names(kegg_gene_list))] # Only human
##     kegg_name_list <- unlist(as.list(KEGGPATHID2NAME)) # All organisms in kegg
##     kegg_name_list <- kegg_name_list[gsub(paste0("^", org), "", names(kegg_gene_list))]
##     names(kegg_gene_list) <- paste0(names(kegg_gene_list), " (", names(kegg_name_list), ") - ", kegg_name_list)
##     return(kegg_gene_list)
## }
## ## Usage:
## keggdb <- load_keggList(org="ath") # org can be: hsa, ath, dme, mmu, ...


## ----reactomedb, eval=FALSE, warning=FALSE, message=FALSE---------------------
## ## Reactome Pathway ID to Hs Entrez ID list
## load_reacList <- function(org="R-HSA") {
##     library(reactome.db)
##     reac_gene_list <- as.list(reactomePATHID2EXTID) # All organisms in reactome
##     reac_gene_list <- reac_gene_list[grepl(org, names(reac_gene_list))] # Only human
##     reac_name_list <- unlist(as.list(reactomePATHID2NAME)) # All organisms in reactome
##     reac_name_list <- reac_name_list[names(reac_gene_list)]
##     names(reac_gene_list) <- paste0(names(reac_gene_list), " (", names(reac_name_list), ") - ", gsub("^.*: ", "", reac_name_list))
##     return(reac_gene_list)
## }
## ## Usage:
## reacdb <- load_reacList(org="R-HSA")


## ----gostats, eval=FALSE, warning=FALSE, message=FALSE------------------------
## ## Load required packages
## library(GOstats); library(GO.db); library(org.At.tair.db)
## ## Define universe and test sample set
## geneUniverse <- keys(org.At.tairGENENAME)
## geneSample <- c("AT2G46210", "AT2G19880", "AT2G38910", "AT5G25140", "AT2G44525")
## ## Generate params object
## params <- new("GOHyperGParams", geneIds = geneSample,
##                 universeGeneIds = geneUniverse,
##                 annotation="org.At.tair", ontology = "MF", pvalueCutoff = 0.5,
##                 conditional = FALSE, testDirection = "over")
## ## Run enrichment test
## hgOver <- hyperGTest(params)
## ## Viewing of results
## summary(hgOver)[1:4,]
## htmlReport(hgOver, file = "MyhyperGresult.html") # html file will be written to current working directory


## ----gohypergall_catdb, eval=FALSE, warning=FALSE, message=FALSE--------------
## ## Create a custom genome-to-GO lookup table for enrichment testing
## library(systemPipeR); library(biomaRt)
## listMarts()  # To choose BioMart database
## listMarts(host = "plants.ensembl.org")
## ## Obtain annotations from BioMart
## listMarts() # To choose BioMart database
## m <- useMart("plants_mart", host = "plants.ensembl.org")
## listDatasets(m)
## m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
## listAttributes(m)  # Choose data types you want to download
## go <- getBM(attributes = c("go_id", "tair_locus", "namespace_1003"), mart = m)
## go <- go[go[, 3] != "", ]; go[, 3] <- as.character(go[, 3])
## go[go[, 3] == "molecular_function", 3] <- "F"; go[go[, 3] == "biological_process", 3] <- "P"; go[go[, 3] == "cellular_component", 3] <- "C"
## go[1:4, ]
## dir.create("./GO")
## write.table(go, "GO/GOannotationsBiomart_mod.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
## catdb <- makeCATdb(myfile = "GO/GOannotationsBiomart_mod.txt", lib = NULL, org = "", colno = c(1, 2, 3), idconv = NULL)
## save(catdb, file="GO/catdb.RData")


## ----gohypergall, eval=FALSE, warning=FALSE, message=FALSE--------------------
## ## Next time catDB can be loaded from file
## load("GO/catdb.RData")
## 
## ## Perform enrichment test on single gene set
## geneids <- unique(as.character(catmap(catdb)$D_MF[,"GeneID"]))
## gene_set_list <- sapply(c("Set1", "Set2", "Set3"), function(x) sample(geneids, 100), simplify=FALSE)
## GOHyperGAll(catdb=catdb, gocat="MF", sample=gene_set_list[[1]], Nannot=2)[1:20,]
## 
## ## Batch analysis of many gene sets for all and slim terms
## goall <- GOCluster_Report(catdb=catdb, setlist=gene_set_list, method="all", id_type="gene", CLSZ=2, cutoff=0.01, gocats=c("MF", "BP", "CC"), recordSpecGO = NULL)
## 
## ## GO Slim analysis by subsetting enrichment results accordingly
## m <- useMart("plants_mart", dataset = "athaliana_eg_gene", host = "plants.ensembl.org")
## goslimvec <- as.character(getBM(attributes = c("goslim_goa_accession"), mart = m)[, 1])
## goslim <- GOCluster_Report(catdb=catdb, setlist=gene_set_list, method="slim",id_type="gene", myslimv=goslimvec, CLSZ=2, cutoff=0.01, gocats = c("MF", "BP", "CC"), recordSpecGO = NULL)
## 
## ## Plot 'GOBatchResult' as bar plot
## goBarplot(goslim, gocat="MF")


## ----fgsea, eval=FALSE, warning=FALSE, message=FALSE--------------------------
## ## Load packages and create sample ranked gene list
## library(fgsea); library(data.table); library(ggplot2); library(org.At.tair.db)
## set.seed(42)
## 
## ## fgsea with KEGG (Arabidopsis)
## geneids <- mappedkeys(org.At.tairCHR)
## exampleRanks <- sort(setNames(sample(seq(-100,100, by=0.001), length(geneids)), geneids))
## fgseaResKegg <- fgsea(pathways=keggdb, stats=exampleRanks, minSize=15, maxSize=500)
## head(fgseaResKegg[order(pval), ])
## plotEnrichment(keggdb[["ath00052 (00052) - Galactose metabolism"]], exampleRanks) + labs(title="Galactose metabolism")
## 
## ## fgsea with Reactome (Human)
## geneids <- unique(as.character(unlist(reacdb)))
## exampleRanks <- sort(setNames(sample(seq(-100,100, by=0.001), length(geneids)), geneids))
## fgseaResReac <- fgsea(pathways=reacdb, stats=exampleRanks, minSize=15, maxSize=500)
## head(fgseaResReac[order(pval), ])
## plotEnrichment(reacdb[["R-HSA-3247509 (R-HSA-3247509) - Chromatin modifying enzymes"]], exampleRanks) + labs(title="Chromatin modifying enzymes")


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

