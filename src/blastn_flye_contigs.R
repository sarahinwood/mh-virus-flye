library(data.table)
library(tidyverse)

blastn_res <- fread("output/10_other_assembled/blastn.outfmt6")

setnames(blastn_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blastn_res, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blastn_res[,.SD[which.min(evalue)], by=query]

#need to make taxid list for bacteria that can be used to subset above
# already have viral but doubt many more viral hits - couple bracovirus hits?

flye_assembled <- read_tsv("output/01_flye/assembly_info.txt")
flye_circular <- subset(flye_assembled, `circ.`=="Y")

hits_for_circular <- subset(min_evalues, query %in% flye_circular$`#seq_name`)

bacterial_taxids <- fread("output/taxids/bacterial_taxids.txt")
flye_bacterial <- subset(min_evalues, taxid %in% bacterial_taxids$V1)

viral_taxids <- fread("output/taxids/species_virus_taxids.txt")
flye_viral <- subset(min_evalues, taxid %in% viral_taxids$V1)

flye_not_micro <- subset(min_evalues, !(taxid %in% viral_taxids$V1 & taxid %in% bacterial_taxids$V1))

## no bacterial hits, no other viral hits

#######################################################
## blastp predicted genes on viralflye viral contigs ##
#######################################################

## viralflye linears
blastp_res_linears <- fread("output/10_other_assembled/other_viruses/linears/prodigal_blastp.outfmt6")
setnames(blastp_res_linears, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blastp_res_linears, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues_linears <- blastp_res_linears[,.SD[which.min(evalue)], by=query]
viral_linear_minev <- subset(min_evalues_linears, taxid %in% viral_taxids$V1)

## viralflye components
blastp_res_comps <- fread("output/10_other_assembled/other_viruses/components//prodigal_blastp.outfmt6")
setnames(blastp_res_comps, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blastp_res_comps, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues_comps <- blastp_res_comps[,.SD[which.min(evalue)], by=query]
viral_comps_minev <- subset(min_evalues_comps, taxid %in% viral_taxids$V1)
## only viral hits are to Megastigmus wasp adintovirus

adintovirus_best_hits <- full_join(viral_linear_minev, viral_comps_minev)
adintovirus_best_hits$contig_id <- tstrsplit(adintovirus_best_hits$query, "2_", keep=1)
flye_contigs_adintovirus <- subset(flye_assembled, `#seq_name` %in% adintovirus_best_hits$contig_id)


# adintovirus genome is linear, 13,305 bp while the contig with hits to it is much longer