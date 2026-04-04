#!/usr/bin/env Rscript
# validate_neg_control_output.R
# Quick validation script to check output consistency

rm(list=ls())
library(Matrix)

# Path to dataset
base_path <- "/home/josep/data/projects/sceptre3/benchmarking/association/neg-control/input_data/replogle-rd7_neg_control_100genes"

cat("Validating negative control dataset outputs...\n\n")

# Load SCEPTRE format files
cat("=== Loading SCEPTRE format ===\n")
response_mat <- readRDS(file.path(base_path, "sceptre/response_matrix.rds"))
grna_mat <- readRDS(file.path(base_path, "sceptre/grna_matrix.rds"))
grna_target_df <- read.csv(file.path(base_path, "sceptre/grna_target_data_frame.csv"))
discovery_pairs <- readRDS(file.path(base_path, "sceptre/discovery_pairs.rds"))
cell_covariates <- read.csv(file.path(base_path, "sceptre/cell_covariates.csv"))

cat("Response matrix:", nrow(response_mat), "genes ×", ncol(response_mat), "cells\n")
cat("gRNA matrix:", nrow(grna_mat), "guides ×", ncol(grna_mat), "cells\n")
cat("Target dataframe:", nrow(grna_target_df), "guides\n")
cat("Discovery pairs:", nrow(discovery_pairs), "pairs\n")
cat("Cell covariates:", nrow(cell_covariates), "cells\n\n")

# Validation 1: Same number of cells across matrices
cat("=== Validation 1: Cell count consistency ===\n")
stopifnot(ncol(response_mat) == ncol(grna_mat))
stopifnot(ncol(response_mat) == nrow(cell_covariates))
cat("✓ All matrices have same number of cells:", ncol(response_mat), "\n\n")

# Validation 2: Guides match between matrix and target dataframe
cat("=== Validation 2: Guide ID consistency ===\n")
stopifnot(nrow(grna_mat) == nrow(grna_target_df))
stopifnot(all(rownames(grna_mat) == grna_target_df$grna_id))
cat("✓ Guide IDs match between matrix and target dataframe\n")
cat("✓ Number of guides:", nrow(grna_mat), "\n\n")

# Validation 3: No guides with 0 cells (BUG FIX!)
cat("=== Validation 3: No empty guides (bug fix check) ===\n")
guide_cell_counts <- Matrix::rowSums(grna_mat)
if (any(guide_cell_counts == 0)) {
  cat("✗ FAILED: Found", sum(guide_cell_counts == 0), "guides with 0 cells!\n")
  cat("   This indicates the old bug is present.\n")
  stop("Validation failed: guides with 0 cells detected")
} else {
  cat("✓ All guides have at least 1 cell\n")
  cat("  Min cells per guide:", min(guide_cell_counts), "\n")
  cat("  Max cells per guide:", max(guide_cell_counts), "\n\n")
}

# Validation 4: Low MOI enforced (every cell has exactly 1 guide)
cat("=== Validation 4: Low MOI enforcement ===\n")
cells_per_guide <- Matrix::colSums(grna_mat)
if (!all(cells_per_guide == 1)) {
  cat("✗ FAILED: Some cells don't have exactly 1 guide\n")
  cat("  Unique values:", unique(cells_per_guide), "\n")
  stop("Validation failed: low MOI not enforced")
} else {
  cat("✓ Every cell has exactly 1 guide (low MOI enforced)\n\n")
}

# Validation 5: Discovery pairs reference valid IDs
cat("=== Validation 5: Discovery pairs validity ===\n")
valid_targets <- all(discovery_pairs$grna_target %in% grna_target_df$grna_target)
valid_genes <- all(discovery_pairs$response_id %in% rownames(response_mat))
stopifnot(valid_targets)
stopifnot(valid_genes)
cat("✓ All discovery pairs reference valid targets and genes\n")
cat("  Expected pairs (Cartesian product):", nrow(grna_target_df) * nrow(response_mat), "\n")
cat("  Actual pairs:", nrow(discovery_pairs), "\n\n")

# Validation 6: Check target names
cat("=== Validation 6: Pseudo-target naming ===\n")
cat("Sample target names:\n")
print(head(grna_target_df, 5))
cat("\nAll targets start with 'non-targeting':",
    all(grepl("^non-targeting", grna_target_df$grna_target)), "\n\n")

# Summary
cat("=== VALIDATION SUMMARY ===\n")
cat("✓ All validations passed!\n")
cat("✓ The refactored version is working correctly\n")
cat("✓ No dimension mismatches detected\n")
cat("✓ Bug fix confirmed: no guides with 0 cells\n")
