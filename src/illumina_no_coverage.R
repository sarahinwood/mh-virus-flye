library(data.table)

no_coverage <- fread("output/09_circos/illumina_no_cov.csv")

no_cov <- no_coverage$locus

findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

nocov_runs <- data.table(paste0(findIntRuns(no_cov)))
nocov_runs$start <- tstrsplit(nocov_runs$V1, "-", keep=1)
nocov_runs$start <- as.numeric(nocov_runs$start)
nocov_runs$end <- tstrsplit(nocov_runs$V1, "-", keep=2)
nocov_runs$end <- as.numeric(nocov_runs$end)
# calculate length
nocov_runs$length <- (nocov_runs$end+1)-nocov_runs$start
mean(nocov_runs$length)
sum(nocov_runs$length)

circos_nocov <- nocov_runs[,c(2,3)]
circos_nocov$label <- paste("main")
circos_nocov$color <- paste("color=red")
circos_nocov <- circos_nocov[,c(3,1,2,4)]
write_tsv(circos_nocov, "output/09_circos/circos_nocov_regions.txt", col_names=F)
