#!/usr/bin/env Rscript
# Wrapper script for sceptre computational benchmarking
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - grna_matrix.rds: Binary gRNA assignment indicator matrix (gRNAs x cells)
# - grna_target_data_frame.csv: gRNA to pseudo-target mapping
# - cell_covariates.csv: Cell-level metadata
# - formula_object.rds: Saved formula for association testing
# - discovery_pairs.rds: Pre-computed Cartesian product discovery pairs
#
# Outputs:
# - association_computational_sceptre.csv: Discovery analysis results

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id <- args[2]

cat("Running sceptre computational benchmarking\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Dataset ID:", dataset_id, "\n")

library(sceptre)

# Load inputs
cat("Loading data...\n")
response_matrix <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna_matrix <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
grna_target_df <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(dataset_dir, "cell_covariates.csv"), stringsAsFactors = FALSE)
formula_string <- readRDS(file.path(dataset_dir, "formula_object.rds"))
formula_object <- as.formula(formula_string)
discovery_pairs <- readRDS(file.path(dataset_dir, "discovery_pairs.rds"))

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:", nrow(grna_target_df), "mappings\n")
cat("  Cell covariates:", nrow(cell_covariates), "cells x", ncol(cell_covariates), "covariates\n")
cat("  Discovery pairs:", nrow(discovery_pairs), "\n")
cat("  Formula loaded from file:", formula_string, "\n")

# Determine which dataset this is (for MOI)
DATASET_NAMES <- c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(name) grepl(name, dataset_id, ignore.case = TRUE))]
if(length(dataset_name) != 1) {
  stop("Could not determine dataset from dataset_id: ", dataset_id)
}
cat("Detected dataset:", dataset_name, "\n")

# Set MOI based on dataset
moi_lookup <- list(
  gasperini = "high",
  replogle = "low"
)
moi <- moi_lookup[[dataset_name]]
cat("MOI:", moi, "\n")

# Import data into sceptre object
cat("Importing data into sceptre object...\n")
scep <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = moi,
  extra_covariates = cell_covariates
)

cat("Setting analysis parameters...\n")
scep <- set_analysis_parameters(
  sceptre_object = scep,
  discovery_pairs = discovery_pairs,
  formula_object = formula_object
)

scep <- assign_grnas(
  sceptre_object = scep,
  method = "thresholding",
  threshold = 1
)

# Very permissive QC for negative control
cat("Running QC (very permissive for computational benchmarking)...\n")
scep <- run_qc(
  sceptre_object = scep,
  n_nonzero_trt_thresh = 0,
  n_nonzero_cntrl_thresh = 0,
  response_n_umis_range = c(0, 1),
  response_n_nonzero_range = c(0, 1),
  p_mito_threshold = 1
)
# cat("Running QC (default settings)...\n")
# scep <- run_qc(
#   sceptre_object = scep
# )

cat("Running discovery analysis...\n")

# cat("SETTING B3=100!!\n")
# scep@B3=100L

start_time <- Sys.time()
scep <- run_discovery_analysis(
  sceptre_object = scep,
  parallel = FALSE,
  output_amount = 2
)
end_time <- Sys.time()

results <- get_result(scep, "run_discovery_analysis")
write.csv(results, "association_computational_sceptre.csv", row.names = FALSE)

cat("Sceptre computational benchmarking complete!\n")
cat("Results saved to: association_computational_sceptre.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total time:", round(elapsed, 2), "seconds\n")
