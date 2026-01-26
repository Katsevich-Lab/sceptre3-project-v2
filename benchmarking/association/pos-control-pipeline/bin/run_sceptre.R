#!/usr/bin/env Rscript
# Wrapper script for sceptre positive control power analysis
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - grna_matrix.rds: Binary gRNA assignment indicator matrix (gRNAs x cells)
# - grna_target_data_frame.csv: gRNA to target mapping
# - cell_covariates.csv: Cell-level metadata
#
# Outputs:
# - results_sceptre.rds: Power check results (data frame)

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id <- args[2]

cat("Running sceptre positive control power analysis\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Dataset ID:", dataset_id, "\n")

library(sceptre)

# Load inputs
cat("Loading data...\n")
response_matrix <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna_matrix <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
grna_target_df <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(dataset_dir, "cell_covariates.csv"), stringsAsFactors = FALSE)

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:", nrow(grna_target_df), "mappings\n")
cat("  Cell covariates:", nrow(cell_covariates), "cells x", ncol(cell_covariates), "covariates\n")

# Define association formula
# TODO this will throw an error if any of these covariates have exact 0s
assoc_fmla <- ~ log(response_n_nonzero) + log(response_n_umis) +
          log(grna_n_umis_subset) + log(grna_n_nonzero_subset)

cat("Association formula:", deparse(assoc_fmla), "\n")
cat("NOTE: batch not currently used for", dataset_id, "\n")

# Import data into sceptre object
cat("Importing data into sceptre object...\n")
scep <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = "low",
  extra_covariates = cell_covariates[,c("grna_n_umis_subset", "grna_n_nonzero_subset")]
)

# Construct positive control pairs (all gRNA -> target gene pairs)
cat("Constructing positive control pairs...\n")
pc_pairs <- construct_positive_control_pairs(scep)
cat("  Positive control pairs:", nrow(pc_pairs), "\n")

scep <- set_analysis_parameters(
  sceptre_object = scep,
  positive_control_pairs = pc_pairs,
  formula_object = assoc_fmla
)

scep <- assign_grnas(
  sceptre_object = scep,
  method = "thresholding",
  threshold = 1
)

# don't actually do any QC
scep <- run_qc(
  sceptre_object = scep,
  n_nonzero_trt_thresh = 0,
  n_nonzero_cntrl_thresh = 0,
  response_n_umis_range = c(0, 1),
  response_n_nonzero_range = c(0, 1),
  p_mito_threshold = 1
)

# Run power check on positive controls
cat("Running power check...\n")
scep <- run_power_check(sceptre_object = scep)

results <- get_result(scep, "run_power_check")
write.csv(results, "association_on_target_sceptre.csv", row.names = FALSE)

cat("Sceptre on-target association testing complete!\n")
cat("Results saved to: association_on_target_sceptre.csv\n")
