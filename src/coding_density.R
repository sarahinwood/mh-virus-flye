library(data.table)
library(ape)

# genome length
genome_length <- 163432

prodigal_gff <- data.frame(read.gff("output/06_prodigal/renamed_gene_predictions.gff"))

##nt length ranges 102-4956 nts, mean nt length =  858.2bp
prodigal_gff$nt_length <- (prodigal_gff$end)-(prodigal_gff$start)+1
mean(prodigal_gff$nt_length)

## AA length range 34-1652, mean AA length = 286.1 aas
prodigal_gff$aa_length <- prodigal_gff$nt_length/3
mean(prodigal_gff$aa_length)

# 82.97% coding density in 163kb total (LbFV 80.0% in 111kb)
(sum(prodigal_gff$nt_length)/genome_length)*100

# no incomplete genes
sum(!(prodigal_gff$partial=="00"))

#37.76% GC content, 64.24% AT content (LbFV 78.7% AT)

# 92 genes on - strand
neg_strand <- subset(prodigal_gff, strand=="-")
#total length = 77.5kb
sum(neg_strand$nt_length)

# 66 genes on + strand
pos_strand <- subset(prodigal_gff, strand=="+")
#total length = 58.1kb
sum(pos_strand$nt_length)
