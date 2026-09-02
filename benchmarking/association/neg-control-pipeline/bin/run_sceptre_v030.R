#!/usr/bin/env Rscript
# Negative-control runner for the historical SCEPTRE v0.3.0
# (Katsevich-Lab/sceptre @ tag v0.3.0). Tests the pipeline's negative-control pairs,
# whose p-values should be approximately uniform if the method is calibrated.
#
# LOW MOI ONLY. High-MOI data is handled by sceptre_manuscript; this errors out on a
# high-MOI dataset rather than silently running the wrong analysis.
#
# v0.3.0 predates the modern sceptre object interface; it exposes a single call,
# run_sceptre_lowmoi(), which internally does gRNA assignment, the precomputation and
# the pairwise testing in one shot.
#
# WHY calibration_check = FALSE. v0.3.0 has a native calibration mode
# (calibration_check = TRUE) in which it constructs its OWN negative-control pairs
# from the gRNA groups. We deliberately do not use it: this pipeline supplies a fixed
# set of negative-control pairs that every method must test, and the modern sceptre
# negative-control run likewise uses its discovery analysis on those same pairs
# rather than its calibration check. Using v0.3.0's native mode would have it testing
# a different set of pairs from every other method.
#
# control_group = "complement" matches what the modern sceptre negative-control run
# forces for both MOI regimes (v0.3.0 would otherwise default to "nt_cells").
#
# NO QC: both pairwise thresholds are set to 0 so that every method sees exactly the
# same cells and genes. v0.3.0 applies no cell-level QC of its own.
#
# SERIAL. v0.3.0 has no internal parallelism (no parallel/n_processors argument, and
# no OpenMP in its C++), and its permutation resample indices are drawn once globally
# and shared across all pairs -- so splitting pairs across workers would change the
# draws. It therefore runs on a single core; allocate wall time accordingly
# (measured: roughly 0.13 s per pair).
#
# REPRODUCIBILITY: v0.3.0 takes no seed argument, so we set the RNG seed here.
#
# Inputs (from dataset_dir): response_matrix.rds, grna_matrix.rds,
# grna_target_data_frame.csv, cell_covariates.csv, formula_object.rds,
# discovery_pairs.rds (the negative-control pairs).
#
# Output: association_neg_control_sceptre_v030.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre) })

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]
seed <- as.integer(Sys.getenv("SCEPTRE_V030_SEED", "1234"))

cat("Running SCEPTRE v0.3.0 (negative control)\n")
cat("  dataset_dir:", dataset_dir, "\n  dataset_id:", dataset_id, "\n")
cat("  sceptre version:", as.character(packageVersion("sceptre")), " seed:", seed, "\n")

# --- MOI guard: this method is low-MOI only ---
DATASET_NAMES <- c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(n) grepl(n, dataset_id, ignore.case = TRUE))]
if (length(dataset_name) != 1) stop("Could not determine dataset from dataset_id: ", dataset_id)
moi <- list(gasperini = "high", replogle = "low")[[dataset_name]]
if (moi != "low") {
  stop("sceptre_v030 is a LOW-MOI method, but dataset '", dataset_id, "' is ", moi,
       " MOI. Use sceptre_manuscript for high-MOI data.")
}
cat("  detected dataset:", dataset_name, " MOI:", moi, "\n")

# --- load inputs ---
response_matrix <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna_matrix     <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
grna_target_df  <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(dataset_dir, "cell_covariates.csv"), stringsAsFactors = FALSE)
formula_object  <- as.formula(readRDS(file.path(dataset_dir, "formula_object.rds")))
discovery_pairs <- readRDS(file.path(dataset_dir, "discovery_pairs.rds"))

# v0.3.0 says "grna_group" where modern sceptre says "grna_target".
grna_group_data_frame <- data.frame(grna_id = as.character(grna_target_df$grna_id),
                                    grna_group = as.character(grna_target_df$grna_target),
                                    stringsAsFactors = FALSE)
response_grna_group_pairs <- data.frame(response_id = as.character(discovery_pairs$response_id),
                                        grna_group = as.character(discovery_pairs$grna_target),
                                        stringsAsFactors = FALSE)

cat("Data loaded:\n")
cat("  response_matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  grna_matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  covariates:", paste(colnames(cell_covariates), collapse = ", "), "\n")
cat("  formula:", deparse(formula_object), "\n")
cat("  pairs:", nrow(response_grna_group_pairs), "\n")
stopifnot(ncol(response_matrix) == ncol(grna_matrix),
          ncol(response_matrix) == nrow(cell_covariates))

# --- run (timed) ---
set.seed(seed)
start_time <- Sys.time()
results <- run_sceptre_lowmoi(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  covariate_data_frame = cell_covariates,
  grna_group_data_frame = grna_group_data_frame,
  response_grna_group_pairs = response_grna_group_pairs,
  formula_object = formula_object,
  calibration_check = FALSE,      # test the supplied negative-control pairs (see header)
  control_group = "complement",   # match the modern negative-control run
  side = "both",
  n_nonzero_trt_thresh = 0L,      # no QC
  n_nonzero_cntrl_thresh = 0L)    # no QC
end_time <- Sys.time()

write.csv(results, "association_neg_control_sceptre_v030.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("SCEPTRE v0.3.0 negative control complete!\n")
cat("Results saved to: association_neg_control_sceptre_v030.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
