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
extra_covariates <- data.frame()


is_simulated_data <- grepl("simulated", input_dir)
if(is_simulated_data) {
	# there should be a covariate data frame then
	extra_covariates <- read.csv(file.path(input_dir, "covariate_data_frame.csv"))
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
  moi = moi,
  extra_covariates = extra_covariates 
)

# Set analysis parameters - use defaults
cat("Setting analysis parameters (using defaults)...\n")
# dummy formula since we're only doing assignment
# TODO: update formula when we actually do assoc testing
sceptre_object <- set_analysis_parameters(sceptre_object, formula = ~1) 

# Build formula dynamically with log transformations and conditional pseudocounts
cat("Building mixture model formula...\n")
has_zero_nonzero <- any(sceptre_object@covariate_data_frame$grna_n_nonzero == 0)
has_zero_umis <- any(sceptre_object@covariate_data_frame$grna_n_umis == 0)

nonzero_term <- if(has_zero_nonzero) "log(grna_n_nonzero + 1)" else "log(grna_n_nonzero)"
umis_term <- if(has_zero_umis) "log(grna_n_umis + 1)" else "log(grna_n_umis)"

# TODO batch is not here right now!
formula_str <- paste0("~ ", nonzero_term, " + ", umis_term)
assign_grnas_formula_object <- as.formula(formula_str)

cat("Using formula:", formula_str, "\n")
cat("  has_zero_nonzero:", has_zero_nonzero, "\n")
cat("  has_zero_umis:", has_zero_umis, "\n")


if(is_simulated_data) {
   cat("Ignore all those other messages. This is simulated data.")
   formula_str <- "~ batch + log(true_grna_n_nonzero + 1) + log(true_grna_n_umis + 1)"
   cat("Actually we will use this formula:", formula_str)
   assign_grnas_formula_object <- as.formula(formula_str)
}

# Assign gRNAs using mixture method with log-transformed covariates
cat("Assigning gRNAs using mixture method...\n")
sceptre_object <- assign_grnas(sceptre_object, method = "mixture", formula_object = assign_grnas_formula_object)

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

# Save assignment matrix as RDS (preserves sparse matrix structure with row/col names)
saveRDS(assignment_matrix, "assignment_matrix_sceptre.rds")
cat("Guide assignment complete. Assignment matrix saved:", nrow(assignment_matrix), "gRNAs x", ncol(assignment_matrix), "cells\n")
cat("Total assignments:", sum(assignment_matrix), "\n")
