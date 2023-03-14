library("data.table")
library("dplyr")
library("rtracklayer")

prodigal_gff <- readGFF("output/06_prodigal/gene_predictions.gff")

# complete genes
sum(prodigal_gff$partial=="00")

# gene length
prodigal_gff$nt_length <- (prodigal_gff$end)-(prodigal_gff$start-1)
prodigal_gff$aa_length <- prodigal_gff$nt_length/3
# aa length
min(prodigal_gff$aa_length)
max(prodigal_gff$aa_length)
mean(prodigal_gff$aa_length)

# start type
sum(prodigal_gff$start_type=="ATG")

# stranded
stranded <- data.table()
stranded$positive <- sum(prodigal_gff$strand=="+")
stranded$negative <- sum(prodigal_gff$strand=="-")

prodigal_gff_score <- prodigal_gff[,-c(6)]

prodigal_gff_score %>% 
  group_by(strand) %>% 
  summarise(CDS_Length = sum(nt_length))


