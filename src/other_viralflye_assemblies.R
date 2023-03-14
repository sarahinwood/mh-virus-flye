library(data.table)
library(dplyr)

## read in
linear_blast <- fread('output/09_other_viruses/linears/prodigal_blastp.outfmt6')
component_blast <- fread('output/09_other_viruses/components/prodigal_blastp.outfmt6')

## set names
setnames(linear_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
setnames(component_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))

## sort hits
setorder(linear_blast, query, evalue, -bit_score, -`%_identical_matches`)
setorder(component_blast, query, evalue, -bit_score, -`%_identical_matches`)
## filter best hits
min_evalues_linear <- linear_blast[,.SD[which.min(evalue)], by=query]
min_evalues_component <- component_blast[,.SD[which.min(evalue)], by=query]
