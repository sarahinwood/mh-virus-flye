library(data.table)
library(tidyverse)

self_blast <- fread("output/09_circos/self_blast/self_blast.outfmt6")

setnames(self_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "annotation"))
self_blast <- subset(self_blast, !(alignment_length == '163432'))
# each hit is present twice - once as query once as subject - so remove the 2nd hit
# should be for any hit where there is another of same length & %id remove second
self_blast$duplicated <- duplicated(self_blast$`%_identical_matches`)
dedup_hits <- subset(self_blast, duplicated==FALSE)

hrs_regions_putative <- subset(dedup_hits, `%_identical_matches`>90&alignment_length>50)
fwrite(hrs_regions_putative, "output/09_circos/self_blast/putative_hrs_regions.csv")

hrs_regions_putative_100_600b <- subset(hrs_regions_putative, alignment_length>100 & alignment_length<600)
mean(hrs_regions_putative_100_600b$alignment_length)

# format for circos
circos_repeats <- hrs_regions_putative_100_600b[,c(7,8,9,10)]
circos_repeats$label1 <- paste("main")
circos_repeats$label2 <- paste("main")
circos_repeats <- circos_repeats[,c(5,1,2,6,3,4)]

write_tsv(circos_repeats, "output/09_circos/circos_hrs_ribbons_100_600.txt", col_names=F)


# how to extract the regions that have all of these hits