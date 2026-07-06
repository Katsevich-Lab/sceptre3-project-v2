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
