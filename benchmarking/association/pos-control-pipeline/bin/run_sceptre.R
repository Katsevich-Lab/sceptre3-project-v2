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

# Determine which dataset this is
DATASET_NAMES <- c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(name) grepl(name, dataset_id, ignore.case = TRUE))]
if(length(dataset_name) != 1) {
  stop("Could not determine dataset from dataset_id: ", dataset_id)
}
cat("Detected dataset:", dataset_name, "\n")

# Define association formula based on dataset
# Using *_full covariates with log(x + 1) transformation
assoc_fmla_lookup <- list(
  replogle = ~ log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) +
               log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1),
  gasperini = ~ log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) +
                log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + prep_batch
)
assoc_fmla <- assoc_fmla_lookup[[dataset_name]]
cat("Association formula:", deparse(assoc_fmla), "\n")

# Set MOI based on dataset
moi_lookup <- list(
  gasperini = "high",
  replogle = "low"
)
moi <- moi_lookup[[dataset_name]]
cat("MOI:", moi, "\n")

# Set covariates to use based on dataset
# Using *_full covariates (not *_subset)
covariates_lookup <- list(
  replogle = c("response_n_nonzero_full", "response_n_umis_full", "grna_n_nonzero_full", "grna_n_umis_full"),
  gasperini = c("response_n_nonzero_full", "response_n_umis_full", "grna_n_nonzero_full", "grna_n_umis_full", "prep_batch")
)
covariates_to_use <- covariates_lookup[[dataset_name]]
cat("Extra covariates:", paste(covariates_to_use, collapse = ", "), "\n")

# Import data into sceptre object
cat("Importing data into sceptre object...\n")
scep <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = moi,
  extra_covariates = cell_covariates[, covariates_to_use]
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
scep <- run_power_check(sceptre_object = scep, output_amount = 2)

results <- get_result(scep, "run_power_check")
write.csv(results, "association_on_target_sceptre.csv", row.names = FALSE)

cat("Sceptre on-target association testing complete!\n")
cat("Results saved to: association_on_target_sceptre.csv\n")
