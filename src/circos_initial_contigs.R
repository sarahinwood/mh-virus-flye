library(data.table)
library(tidyverse)

## minimap_res

bed <- fread("output/09_circos/map_initial_MhV_MhV1/sorted_Mh_DNA_virus_contigs.bed")
bed$label <- paste("main")
bed$color <- ifelse(bed$V6 =="+", paste("color=red"), paste("color=black"))
circos_minimap <- bed[,c(7,2,3,8)]
write_tsv(circos_minimap, "output/09_circos/circos_initial_contigs/circos_initial_contigs_MhV1_mapped.txt", col_names=F)

  
bed_ex <- fread("output/09_circos/map_initial_MhV_MhV1/sorted_extended_Mh_DNA_virus_contigs.bed")
bed_ex$label <- paste("main")
bed_ex$color <- ifelse(bed_ex$V6 =="+", paste("color=orange"), paste("color=blue"))
circos_minimap_ex <- bed_ex[,c(7,2,3,8)]
write_tsv(circos_minimap_ex, "output/09_circos/circos_initial_contigs/circos_initial_contigs_extended_MhV1_mapped.txt", col_names=F)
