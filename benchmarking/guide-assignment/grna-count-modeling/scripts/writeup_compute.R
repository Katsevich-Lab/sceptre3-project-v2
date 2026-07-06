#!/usr/bin/env Rscript
# ============================================================================
# Compute backbone for the fishash / fishash+ / CLEANSER method-comparison writeup.
# Runs fishash (cell_margin="observed") and fishash+ (cell_margin="ambient") on every real dataset,
# plus CLEANSER's cheap <=2 ambient DEPTH (no MCMC). Real CLEANSER ASSIGNMENTS (MCMC) are produced
# separately by run_cleanser_batch.sh and read at render time from results/ambient_ceiling/cleanser_calls/.
#
# Saves to results/ambient_ceiling/writeup/ :
#   assign/<ds>_{fishash,fishashplus}.rds   sparse assignment matrices (for three-way Jaccard later)
#   depths_sampled.csv     per-cell depth (fishash rank-1, fishash+ rank-1, CLEANSER <=2), subsampled
#   depth_summary.csv      per-dataset median depths + spearman(f+ , cleanser)
#   assign_counts.csv      per-dataset per-method assignment counts + MOI
#   ambient_frac_by_count.csv   per (dataset,k,method): #ambient entries at count k / #entries at k
#   jaccard_ff.csv         per-dataset fishash-vs-fishash+ Jaccard (overall) + by count-bin
#   ceilings.csv           per-guide ambient ceiling (fishash, fishash+)
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats) })
source("scripts/contingency_method.R")
set.seed(1)
OUT <- "results/ambient_ceiling/writeup"; AOUT <- file.path(OUT,"assign")
dir.create(AOUT, recursive=TRUE, showWarnings=FALSE)

source("scripts/datasets.R")
paths <- dataset_paths()
