library(data.table)
library(dplyr)
library(ape)

original_gff <- read.gff("output/06_prodigal/renamed_gene_predictions.gff")
viral_taxids <- fread("data/species_virus_taxids.txt", header=F)
blast_res <- fread("output/07_prodigal_blastp/prodigal_blastp.outfmt6")

setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))

lbfv <- subset(blast_res, grepl("Leptopilina boulardi filamentous virus", blast_res$annotation))
dafv <- subset(blast_res, grepl("Drosophila-associated filamentous virus", blast_res$annotation))

# make table of number of hits per taxid
blast_res$taxid <- tstrsplit(blast_res$taxid, ";", keep=1)
blast_res$species <- tstrsplit(blast_res$annotation, "[", fixed=T, keep=2)
blast_res$species <- tstrsplit(blast_res$species, "]", fixed=T, keep=1)
# unique gene hits by taxid
unique_genes_taxid <- unique(blast_res[,c('query','taxid', 'species')])
taxid_species <- unique(unique_genes_taxid[,c('taxid', 'species')])
# number of different genes per taxid with hits
taxid_occurances <- unique_genes_taxid %>% count(taxid)
taxid_occurances_species <- merge(taxid_occurances, taxid_species)

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res, query, evalue, -bit_score, -`%_identical_matches`)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blast_res[,.SD[which.min(evalue)], by=query]
sum(min_evalues$taxid=="552509") # 23/34 total lbfv hits are best hit
sum(min_evalues$taxid=="2743186") # 8/28 total dafv hits are best hit

# viral vs non-viral hits
min_evalues$viral_status <- ifelse(min_evalues$taxid %in% viral_taxids$V1, "Viral", "Non-viral")
sum(min_evalues$viral_status=="Viral") # 55
sum(min_evalues$viral_status=="Non-viral") # 19

fwrite(min_evalues, "output/07_prodigal_blastp/best_hits.csv")

## monooxygenase hits
bacterial_taxids <- fread("output/taxids/bacterial_taxids.txt")
mono_orfs <- c("ORF116", "ORF133")
orf116_orf133_res <- subset(blast_res, blast_res$query %in% mono_orfs)
orf116_orf133_res$bacterial <- ifelse(orf116_orf133_res$taxid %in% bacterial_taxids$V1, "bacterial", "not bacterial")

# hits that were viral but had a high number of nonviral

blast_res_viral_not <- blast_res
blast_res_viral_not$viral_status <- ifelse(blast_res_viral_not$taxid %in% viral_taxids$V1, "Viral", "Non-viral")
viral_not_by_gene <- blast_res_viral_not %>% count(query, viral_status)
viral_not_by_gene_wide <- dcast(viral_not_by_gene, query ~ viral_status)

viral_not_annot <- merge(viral_not_by_gene_wide, min_evalues, by="query")
viral_not_annot <- viral_not_annot[,c(1,2,3,16,17)]
