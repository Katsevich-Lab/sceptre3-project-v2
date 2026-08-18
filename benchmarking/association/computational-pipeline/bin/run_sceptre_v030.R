#!/usr/bin/env Rscript
# Computational-pipeline runner for the historical SCEPTRE v0.3.0
# (Katsevich-Lab/sceptre @ tag v0.3.0). Runs the v0.3.0-era method over the
# discovery pairs and reports its runtime.
#
# LOW MOI ONLY. This is the historical SCEPTRE we benchmark on low-MOI data
# (Replogle); high-MOI data is benchmarked with the 2021 manuscript version
# instead (run_sceptre_manuscript.R). The script errors out on a high-MOI dataset
# rather than silently running the wrong analysis.
#
# v0.3.0 predates the modern sceptre object interface (import_data ->
# set_analysis_parameters -> assign_grnas -> run_qc -> run_discovery_analysis).
# Instead it exposes a single call, run_sceptre_lowmoi(), which internally does
# gRNA assignment, the precomputation, and the pairwise testing in one shot.
#
# NO QC: to keep every method in this benchmark operating on exactly the same
# cells and genes, both pairwise QC thresholds are set to 0. v0.3.0 applies no
# cell-level QC of its own, so this disables all of its filtering.
#
# Covariates: v0.3.0 builds its design matrix with
# stats::model.matrix(formula_object, covariate_data_frame) and does NOT compute
# any covariates of its own -- so it uses exactly the columns in
# cell_covariates.csv, and the log() transforms written into formula_object.rds
# are applied the same way as for the other methods.
#
# Naming: the method key is `sceptre_v030` (no dots) because Nextflow branch
# labels and process names must be valid identifiers.
#
# Inputs (from dataset_dir, same layout as the modern sceptre method):
# - response_matrix.rds          genes x cells sparse matrix
# - grna_matrix.rds              gRNAs x cells binary assignment matrix
# - grna_target_data_frame.csv   grna_id -> grna_target mapping
# - cell_covariates.csv          per-cell covariates
# - formula_object.rds           covariate formula string
# - discovery_pairs.rds          grna_target x response_id pairs to test
#
# Output:
# - association_computational_sceptre_v030.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre) })

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]

cat("Running SCEPTRE v0.3.0 (computational)\n")
cat("  dataset_dir:", dataset_dir, "\n  dataset_id:", dataset_id, "\n")
cat("  sceptre version:", as.character(packageVersion("sceptre")), "\n")

# --- determine dataset / MOI; this method is low-MOI only ---
DATASET_NAMES <- c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(name) grepl(name, dataset_id, ignore.case = TRUE))]
if (length(dataset_name) != 1) {
  stop("Could not determine dataset from dataset_id: ", dataset_id)
}
moi_lookup <- list(gasperini = "high", replogle = "low")
moi <- moi_lookup[[dataset_name]]
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
# calibration_check = FALSE -> discovery analysis (not the negative-control check).
# resampling_mechanism = "permutations" matches the modern sceptre computational run.
# control_group is left at the wrapper's own default (nt_cells for low MOI).
start_time <- Sys.time()
results <- run_sceptre_lowmoi(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  covariate_data_frame = cell_covariates,
  grna_group_data_frame = grna_group_data_frame,
  response_grna_group_pairs = response_grna_group_pairs,
  formula_object = formula_object,
  calibration_check = FALSE,
  resampling_mechanism = "permutations",
  side = "both",
  n_nonzero_trt_thresh = 0L,      # no QC
  n_nonzero_cntrl_thresh = 0L)    # no QC
end_time <- Sys.time()

write.csv(results, "association_computational_sceptre_v030.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("SCEPTRE v0.3.0 complete!\n")
cat("Results saved to: association_computational_sceptre_v030.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
