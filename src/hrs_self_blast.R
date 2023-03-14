library(data.table)
library(dplyr)

self_blast <- fread("output/09_circos/self_blast/self_blast.outfmt6")
setnames(self_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

self_blast <- subset(self_blast, alignment_length!=163432)
# seqid > 90 and length > 50 bp
hrs_regions <- subset(self_blast, `%_identical_matches`>90&alignment_length>50)

# remove first hit which is just whole genome against itself & last column
self_blast <- self_blast[-c(1),-c(13)]

## what genes are near these regions? too many to look at, no clear relationship with bro genes from first 5ish lines