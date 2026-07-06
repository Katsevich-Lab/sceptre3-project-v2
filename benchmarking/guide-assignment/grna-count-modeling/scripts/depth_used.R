#!/usr/bin/env Rscript
# Per-cell ambient depth USED by each method (the cell exposure it tests a count against):
#   fishash  = raw library size  N_{:,c} = colSums(counts)
#   fishash+ = denoised rank-1 depth  d_c  (cached fishash+ Cc = colSums of its rank-1 completion)
#   CLEANSER = L_i = sum of the cell's <=2 counts
# Overwrites the writeup depth CSVs so method_comparison.qmd renders depth-USED (not the rank-1 field
# fishash merely computes). Cheap: reuses cached fishash+ Cc, everything else is a colSums.
suppressPackageStartupMessages({ library(Matrix) })
set.seed(1)
W <- "results/ambient_ceiling/writeup"; CACHE <- "results/ambient_ceiling/fit_cache"
source("scripts/datasets.R")
paths <- dataset_paths()

depths <- list(); dsum <- list()
for (nm in names(paths)) {
  counts <- load_grna_matrix(paths[[nm]])
  rawlib <- Matrix::colSums(counts)                                    # fishash depth used
  l2 <- kappa_leq2(counts)                                             # CLEANSER depth used
  fit <- readRDS(file.path(CACHE, paste0(nm,"_fit.rds"))); Cc <- as.numeric(fit$Cc)  # fishash+ depth used
  Cc <- Cc[match(colnames(counts), names(fit$Cc))]
  keep <- rawlib>0 | Cc>0 | l2>0; idx <- which(keep); if (length(idx) > 6000) idx <- sample(idx, 6000)
  depths[[nm]] <- data.frame(dataset=nm, fishash=rawlib[idx], fishashplus=Cc[idx], cleanser=l2[idx])
  dsum[[nm]] <- data.frame(dataset=nm, guides=nrow(counts), cells=ncol(counts), med_lib=median(rawlib[rawlib>0]),
    med_depth_fishash=median(rawlib[rawlib>0]), med_depth_fishashplus=median(Cc[Cc>0]),
    med_depth_cleanser=median(l2[l2>0]),
    spearman_fp_clns=suppressWarnings(cor(Cc[Cc>0 & l2>0], l2[Cc>0 & l2>0], method="spearman")))
  cat(sprintf("%-20s rawlib=%.0f  rank1=%.0f  <=2=%.0f\n", nm, dsum[[nm]]$med_depth_fishash,
      dsum[[nm]]$med_depth_fishashplus, dsum[[nm]]$med_depth_cleanser)); flush.console()
  rm(counts); gc()
}
write.csv(do.call(rbind,depths), file.path(W,"depths_sampled.csv"), row.names=FALSE)
write.csv(do.call(rbind,dsum),   file.path(W,"depth_summary.csv"), row.names=FALSE)
cat("\nwrote depth-USED CSVs to", W, "\n")
