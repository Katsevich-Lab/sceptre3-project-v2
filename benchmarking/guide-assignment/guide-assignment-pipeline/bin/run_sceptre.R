#!/usr/bin/env Rscript

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]  # Directory containing sceptre-compatible data
dataset_id <- args[2] # Dataset identifier

# Load sceptre library (pre-installed in container)
library(sceptre)

# Load the sceptre-compatible data files
cat("Loading data from:", input_dir, "\n")

# Load sparse matrices and data frames
library(Matrix)
response_matrix <- readRDS(file.path(input_dir, "response_matrix.rds"))
grna_matrix <- readRDS(file.path(input_dir, "grna_matrix.rds"))
grna_target_df <- read.csv(file.path(input_dir, "grna_target_data_frame.csv"))
if(!any("non_targeting" %in% grna_target_df$grna_target)) {
	grna_target_df[1,"grna_target"] <- "non-targeting"
	warning("No NTs present. `grna_target_df` modified to have an NT")
}

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:", nrow(grna_target_df), "gRNAs mapped to targets\n")

# MOI lookup table based on dataset
moi_lookup <- list(
  gasperini = "high",
  replogle = "low"
)

# Determine MOI from dataset_id
dataset_key <- if (grepl("gasperini", dataset_id, ignore.case = TRUE)) {
  "gasperini"
} else if (grepl("replogle", dataset_id, ignore.case = TRUE)) {
  "replogle"
} else {
  stop("Unknown dataset '", dataset_id, "'. Expected 'replogle' or 'gasperini' in dataset name.")
}

moi <- moi_lookup[[dataset_key]]
cat("Dataset:", dataset_id, "-> MOI:", moi, "\n")


# Create sceptre object
cat("Creating sceptre object...\n")
sceptre_object <- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = moi  # Assume low MOI for now
)

# Set analysis parameters - use defaults
cat("Setting analysis parameters (using defaults)...\n")
# dummy formula since we're only doing assignment
# TODO: update formula when we actually do assoc testing
sceptre_object <- set_analysis_parameters(sceptre_object, formula = ~1) 

# Assign gRNAs using thresholding method (more robust for dummy data)
cat("Assigning gRNAs using maximum method...\n")
sceptre_object <- assign_grnas(sceptre_object, method = "mixture", formula_object = ~ grna_n_nonzero + grna_n_umis)

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
cells_with_assignments <- assignment_matrix |> apply(2, which)
output_df <- data.frame(
  cell_id = colnames(grna_matrix),
  grna_id = rownames(assignment_matrix)[cells_with_assignments],
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
