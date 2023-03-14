#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)
library(zoo)

###########
# GLOBALS #
###########

nanopore <- snakemake@input[["nanopore"]]
illumina <- snakemake@input[["illumina"]]

########
# MAIN #
########

nanoporecov <- fread(nanopore)
setnames(nanoporecov, old=c("V1", "V2", "V3"), new=c("Chr", "locus", "depth"))

nanopore_sliding_window_cov <- nanoporecov %>% do(data.frame(window.start = rollapply(.$locus, width=500, by=100, FUN=min, align="left"),
                           window.end = rollapply(.$locus, width=500, by=100, FUN=max, align="left"),
                           coverage = rollapply(.$depth, width=500, by=100, FUN=mean, align="left")))
# calc mean for part between end of sliding window and end of chr
remainingnanopore <- subset(nanoporecov, locus > 163400)
nanopore_remaining_window <- remainingnanopore %>% do(data.frame(window.start = rollapply(.$locus, width=32, by=100, FUN=min, align="left"),
                                                             window.end = rollapply(.$locus, width=32, by=100, FUN=max, align="left"),
                                                             coverage = rollapply(.$depth, width=32, by=100, FUN=mean, align="left")))
nanopore_circos <- full_join(nanopore_sliding_window_cov, nanopore_remaining_window)
nanopore_circos$label <- paste("main")
nanopore_circos$fill_colour <- paste("fill_color=viridis_teal")
nanopore_circos <- nanopore_circos[,c(4,1,2,3,5)]
write_tsv(nanopore_circos, snakemake@output[["nanopore"]], col_names=F)

illuminacov <- fread(illumina)
setnames(illuminacov, old=c("V1", "V2", "V3"), new=c("Chr", "locus", "depth"))

illumina_sliding_window_cov <- illuminacov %>% do(data.frame(window.start = rollapply(.$locus, width=500, by=100, FUN=min, align="left"),
                                            window.end = rollapply(.$locus, width=500, by=100, FUN=max, align="left"),
                                            coverage = rollapply(.$depth, width=500, by=100, FUN=mean, align="left")))
# calc mean for part between end of sliding window and end of chr
remainingillumina <- subset(illuminacov, locus > 163400)
illumina_remaining_window <- remainingillumina %>% do(data.frame(window.start = rollapply(.$locus, width=32, by=100, FUN=min, align="left"),
                                                                 window.end = rollapply(.$locus, width=32, by=100, FUN=max, align="left"),
                                                                 coverage = rollapply(.$depth, width=32, by=100, FUN=mean, align="left")))
illumina_circos <- full_join(illumina_sliding_window_cov, illumina_remaining_window)
illumina_circos$label <- paste("main")
illumina_circos$fill_colour <- paste("fill_color=viridis_violet")
illumina_circos <- illumina_circos[,c(4,1,2,3,5)]
write_tsv(illumina_circos, snakemake@output[["illumina"]], col_names=F)

# write regions with no illumina coverage
illumina_nocov <- subset(illuminacov, depth==0)
fwrite(illumina_nocov, snakemake@output[["illumina_no_cov"]])

# write log
sessionInfo()
