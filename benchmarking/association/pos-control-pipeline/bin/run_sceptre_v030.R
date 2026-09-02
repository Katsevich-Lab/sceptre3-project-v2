#!/usr/bin/env Rscript
# Positive-control runner for the historical SCEPTRE v0.3.0
# (Katsevich-Lab/sceptre @ tag v0.3.0). Tests each target against its OWN gene, so the
# p-values measure the method's power to detect a true effect.
#
# v0.3.0 has no dedicated power-check routine, but a power check is simply the
# ordinary association test applied to on-target pairs: we construct the
# positive-control pairs (every target that is itself a measured gene, paired with
# that gene) and run the standard discovery test on them. This mirrors what the
# modern sceptre positive-control run does with construct_positive_control_pairs().
#
# LOW MOI ONLY. High-MOI data is handled by sceptre_manuscript.
#
# calibration_check = FALSE because this is an association test on supplied pairs,
# not v0.3.0's negative-control calibration mode.
#
# control_group is left at the wrapper's default (nt_cells for low MOI), matching the
# modern positive-control run, which does not override it either. (Note this differs
# from the negative-control pipeline, where the modern run forces "complement".)
#
# NO QC: both pairwise thresholds are set to 0, and the modern positive-control run
# skips QC entirely, so all methods see the same cells and genes.
#
# SERIAL: v0.3.0 has no internal parallelism (no parallel/n_processors argument, no
# OpenMP), so this runs on one core regardless of the cpus allocation.
#
# REPRODUCIBILITY: v0.3.0 takes no seed argument, so we set the RNG seed here.
#
# Inputs (from dataset_dir): response_matrix.rds, grna_matrix.rds,
# grna_target_data_frame.csv, cell_covariates.csv. (No formula_object.rds or
# discovery_pairs.rds in this pipeline.)
#
# Output: association_on_target_sceptre_v030.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre) })

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]
seed <- as.integer(Sys.getenv("SCEPTRE_V030_SEED", "1234"))

cat("Running SCEPTRE v0.3.0 (positive control)\n")
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

# Covariate formula, matching the modern positive-control run for this dataset.
formula_object <- ~ log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) +
                    log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1)
cat("  formula:", deparse(formula_object), "\n")

# v0.3.0 says "grna_group" where modern sceptre says "grna_target".
grna_group_data_frame <- data.frame(grna_id = as.character(grna_target_df$grna_id),
                                    grna_group = as.character(grna_target_df$grna_target),
                                    stringsAsFactors = FALSE)

# Positive-control pairs: every target that is itself a measured gene, paired with it.
pc_targets <- intersect(unique(as.character(grna_target_df$grna_target)), rownames(response_matrix))
response_grna_group_pairs <- data.frame(response_id = pc_targets, grna_group = pc_targets,
                                        stringsAsFactors = FALSE)

cat("Data loaded:\n")
cat("  response_matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  grna_matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  covariates:", paste(colnames(cell_covariates), collapse = ", "), "\n")
cat("  positive control pairs:", nrow(response_grna_group_pairs), "\n")
if (nrow(response_grna_group_pairs) == 0L) {
  stop("No positive control pairs: no grna_target is a row of response_matrix.")
}
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
  calibration_check = FALSE,
  side = "both",
  n_nonzero_trt_thresh = 0L,      # no QC
  n_nonzero_cntrl_thresh = 0L)    # no QC
end_time <- Sys.time()

write.csv(results, "association_on_target_sceptre_v030.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("SCEPTRE v0.3.0 positive control complete!\n")
cat("Results saved to: association_on_target_sceptre_v030.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
