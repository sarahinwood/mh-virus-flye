library(data.table)
library(tidyverse)

blastn_res <- fread("output/10_other_assembled/blastp.outfmt6")

setnames(blastn_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blastn_res, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blastn_res[,.SD[which.min(evalue)], by=query]
min_evalues$contig <- tstrsplit(min_evalues$query, "_", keep=2)
min_evalues$contig <- paste("contig", min_evalues$contig, sep="_")

# bacterial
bacterial_taxids <- fread("output/taxids/bacterial_taxids.txt")
flye_bacterial <- subset(min_evalues, taxid %in% bacterial_taxids$V1)
flye_bacterial_notMhV1 <- subset(flye_bacterial, !(grepl("contig_859", query)))

viral_taxids <- fread("output/taxids/species_virus_taxids.txt")
flye_viral <- subset(min_evalues, taxid %in% viral_taxids$V1)
flye_viral_notMhV1 <- subset(flye_viral, !(grepl("contig_859", query)))

# circular
flye_assembled <- read_tsv("output/01_flye/assembly_info.txt")
assembled_circular <- subset(flye_assembled, `circ.`=="Y")
flye_circular <- subset(min_evalues, contig %in% assembled_circular$`#seq_name`)
# circular viral or bacterial
flye_circular_virus <- subset(flye_circular, taxid %in% viral_taxids$V1) # only 1 hit that isn't on MhV1
flye_circular_bacterial <- subset(flye_circular, taxid %in% bacterial_taxids$V1) # only 1 hit that isn't on MhV1

# adintovirus is the ony other species with multiple hits really