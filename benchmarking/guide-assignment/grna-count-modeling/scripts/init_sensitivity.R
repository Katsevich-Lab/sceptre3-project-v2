#!/usr/bin/env Rscript
# ============================================================================
# Sensitivity of fishash+ to the INITIALIZATION of the assigned set A.
#   default : A_0 = {}          -> background_1 = raw counts (initial depth = raw library size)
#   clns    : A_0 = {N_gc > 2}  -> background_1 masks >2 (initial depth = CLEANSER-style <=2 depth)
# Both then run the identical fishash+ iteration to a fixed point. We compare the FINAL assigned
# sets: if they coincide (Jaccard ~ 1) the fixed point is initialization-insensitive.
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2) })
source("scripts/contingency_method.R")
OUT <- "results/ambient_ceiling"

source("scripts/datasets.R")
paths <- dataset_paths()
