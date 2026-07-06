#!/usr/bin/env Rscript
# ============================================================================
# Per-cell ambient depth: fishash+ (rank-1 completion) vs CLEANSER (sum of <=2 UMIs).
# Scatter per dataset.
#   fishash+ depth  d_c^F = Cc[c]                 (cached rank-1 cell factor; colSums of noise)
#   CLEANSER depth  d_c^C = colSums(counts * (counts <= 2))   (sub-threshold, signal-free UMIs)
# fishash+ depth read from cache (no refit); CLEANSER depth is a colSums over the raw matrix.
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
set.seed(1)
OUT   <- "results/ambient_ceiling"; CACHE <- file.path(OUT, "fit_cache")

source("scripts/datasets.R")
paths <- dataset_paths()
