library(data.table)
library(tidyverse)

virus_taxids <- fread("data/species_virus_taxids.txt", header=F)

## linear ##
linear_blast <- fread('output/09_other_viruses/linears/prodigal_blastp.outfmt6')
setnames(linear_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
          new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
                "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(linear_blast, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
linear_min_evalues <- linear_blast[,.SD[which.min(evalue)], by=query]


## components ##
component_blast <- fread('output/09_other_viruses/components/prodigal_blastp.outfmt6')
setnames(component_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(component_blast, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
component_min_evalues <- component_blast[,.SD[which.min(evalue)], by=query]

## other circular? ##
circular_blast <- fread("output/unpolished_prodigal_blastp/prodigal_blastp.outfmt6")
setnames(circular_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(circular_blast, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
circular_min_evalues <- circular_blast[,.SD[which.min(evalue)], by=query]
other_circular_min_evalues <- subset(circular_min_evalues, !grepl("contig_859", query)) # no others

## viral hits ##
linear_viral <- subset(linear_min_evalues, taxid %in% virus_taxids$V1)
component_viral <- subset(component_min_evalues, taxid %in% virus_taxids$V1)
## only Megastigmus wasp adintovirus hits
