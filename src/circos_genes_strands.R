#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)
library(rtracklayer)

###########
# GLOBALS #
###########

gff <- snakemake@input[["gff"]]
annots <- snakemake@input[["annots"]]

########
# MAIN #
########

annots <- fread(annots)
prodigal_gff <- data.table(data.frame(readGFF(gff)))
prodigal_gff$ORF_ID <- tstrsplit(prodigal_gff$ID, "_", keep=2)
prodigal_gff$ORF_ID <- paste("ORF", prodigal_gff$ORF_ID, sep="")

# sort out bro genes
BRO_annot <- subset(annots, grepl("bro|baculovirus repeat ORF ", annotation, ignore.case=T))
BRO_desc <- subset(annots, grepl("bro|baculovirus repeat ORF ", description, ignore.case=T))
all_BRO <- full_join(BRO_annot, BRO_desc)
# positive
positive <- subset(prodigal_gff, strand=="+")
positive$color <- ifelse(positive$ORF_ID %in% all_BRO$`#gene_id`,
                         paste("color=viridis_dyellow"),  paste("color=viridis_yellow"))
positive_circos <- positive[,c(4,5,23)]
positive_circos$label <- paste("main")
positive_circos  <- positive_circos[,c(4,1,2,3)]
# negative
negative <- subset(prodigal_gff, strand=="-")
negative$color <- ifelse(negative$ORF_ID %in% all_BRO$`#gene_id`,
                         paste("color=viridis_dgreen"),  paste("color=viridis_green"))
negative_circos <- negative[,c(4,5,23)]
negative_circos$label <- paste("main")
negative_circos  <- negative_circos[,c(4,1,2,3)]
# all
circos_genes <- full_join(positive_circos, negative_circos)
write_tsv(circos_genes, snakemake@output[["circos_genes"]], col_names=F)

# write log
sessionInfo()