library(data.table)
library(tidyverse)
library(ape)

hmmer_res <- fread("output/08_hmmer/prodigalPFAM_domtblout.csv")
length(unique(hmmer_res$query))
# domains in 108/158 genes

# filter by individual e-value - 5e-04 (same as contig blast)
hmmscan_sig <- filter(hmmer_res, `domain_i-evalue`<5e-04) # only 75 sig hits
# how many genes have predicted domains - 43
length(unique(hmmscan_sig$query))
# how many unique domain predictions - 43
length(unique(hmmscan_sig$target))

##filter important columns out
hmmscan_sig_table <- hmmscan_sig[,c(4,2,23,18,19,13,12)]
fwrite(hmmscan_sig_table, "output/08_hmmer/hmmer_sig_domains.csv")

prodigal_gff <- data.table(data.frame(read.gff("output/06_prodigal/renamed_gene_predictions.gff")))
prodigal_gff$id <- tstrsplit(prodigal_gff$attributes, "ID=", keep=2)
prodigal_gff$id <- tstrsplit(prodigal_gff$id, ";", keep=1)
gff_table <- prodigal_gff[,c(10,4,5,7)]
fwrite(gff_table, "output/06_prodigal/gff_table.csv")


## do bro genes line up with repeats
sig_bro <- subset(hmmscan_sig_table, grepl("BRO", hmmscan_sig_table$description))
gff_bro <- subset(prodigal_gff, prodigal_gff$id %in% sig_bro$query)

