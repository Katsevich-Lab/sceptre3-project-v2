#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]  # Directory containing sceptre-compatible data

# input_dir <- "/home/josep/data/projects/sceptre3/benchmarking/guide_assignment/input_data/sceptre_example/sceptre_input/"

# Initialize renv and restore from lockfile if it exists
# Use a simple working directory approach for lockfile
lockfile_path <- "renv.lock"

if (file.exists(lockfile_path)) {
  # Restore from existing lockfile (reproducible mode)
  cat("Restoring R environment from lockfile...\n")
  renv::restore(lockfile = lockfile_path, prompt = FALSE)
} else {
  # First time setup: install sceptre and create lockfile
  cat("First time setup: installing sceptre and creating lockfile...\n")
  renv::init(bare = TRUE)
  
  # Install remotes first (needed for GitHub installations)
  install.packages("remotes", repos = "https://cran.r-project.org")
  
  # Install data.table first from CRAN (required dependency for sceptre)
  install.packages("data.table", repos = "https://cran.r-project.org")
  
  # Install sceptre with pinned version
  remotes::install_github("Katsevich-Lab/sceptre", ref = "v0.10.1")
  
  # Create lockfile for future reproducibility
  renv::snapshot(lockfile = lockfile_path, prompt = FALSE)
  cat("Lockfile created at:", lockfile_path, "\n")
}

# Load sceptre
library(sceptre)

# Load the sceptre-compatible data files
cat("Loading data from:", input_dir, "\n")

# Load matrices and data frames
response_matrix <- read.csv(file.path(input_dir, "response_matrix.csv"), row.names = 1)
grna_matrix <- read.csv(file.path(input_dir, "grna_matrix.csv"), row.names = 1)
grna_target_df <- read.csv(file.path(input_dir, "grna_target_data_frame.csv"))

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:", nrow(grna_target_df), "gRNAs mapped to targets\n")

# Create sceptre object
cat("Creating sceptre object...\n")
sceptre_object <- import_data(
  response_matrix = as.matrix(response_matrix),
  grna_matrix = as.matrix(grna_matrix),
  grna_target_data_frame = grna_target_df,
  moi = "low"  # Assume low MOI for now
)

# Set analysis parameters - use defaults
cat("Setting analysis parameters (using defaults)...\n")
sceptre_object <- set_analysis_parameters(sceptre_object)

# Assign gRNAs using thresholding method (more robust for dummy data)
# TODO: Change back to mixture method once realistic data is available
cat("Assigning gRNAs using thresholding method...\n")
sceptre_object <- assign_grnas(sceptre_object, method = "thresholding", threshold = 9)

# Extract guide assignments as sparse logical matrix
assignment_matrix <- get_grna_assignments(sceptre_object)
cat("Assignment matrix dimensions:", nrow(assignment_matrix), "gRNAs x", ncol(assignment_matrix), "cells\n")

# Debug the assignment matrix structure
debug_info <- paste(
  "Assignment matrix class:", paste(class(assignment_matrix), collapse = ", "), "\n",
  "Row names length:", length(rownames(assignment_matrix)), "\n",
  "Col names length:", length(colnames(assignment_matrix)), "\n", 
  "Sum of assignments:", sum(assignment_matrix), "\n",
  "First 5 row names:", paste(head(rownames(assignment_matrix), 5), collapse = ", "), "\n",
  "First 5 col names:", paste(head(colnames(assignment_matrix), 5), collapse = ", "), "\n"
)
writeLines(debug_info, "debug_assignment_matrix.txt")
cat(debug_info)

# Create standardized output using cell names from sceptre object
output_df <- data.frame(
  cell_id = colnames(grna_matrix),
  grna_id = rownames(assignment_matrix)[assignment_matrix |> apply(2, which)],
  stringsAsFactors = FALSE
)

# Remove any unassigned cells (marked as NA)
assigned_cells <- output_df[!is.na(output_df$grna_id), ]

# Handle case where no assignments were made
if(nrow(assigned_cells) == 0) {
  cat("Warning: No cells were assigned to any gRNA with current threshold.\n")
  assigned_cells <- data.frame(cell_id = character(0), grna_id = character(0))
}

output_df <- assigned_cells

cat("Guide assignment complete:", nrow(output_df), "cells assigned\n")

# Write standardized output
write.csv(output_df, "assignments_sceptre.csv", row.names = FALSE)