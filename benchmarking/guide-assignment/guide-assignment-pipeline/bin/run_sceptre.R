#!/usr/bin/env Rscript

library(sceptre)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]  # Directory containing sceptre-compatible data
dataset_id <- args[2] # Dataset identifier


# determine which dataset this is
DATASET_NAMES = c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(name) grepl(name, dataset_id, ignore.case = TRUE))]
if(length(dataset_name) != 1) {
	stop("Could not determine dataset from the input path.")
}

# Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Loading data from:", input_dir, "\n")

response_matrix <- readRDS(file.path(input_dir, "response_matrix.rds"))
grna_matrix <- readRDS(file.path(input_dir, "grna_matrix.rds"))
grna_target_df <- read.csv(file.path(input_dir, "grna_target_data_frame.csv"))
if(!any(grna_target_df$grna_target == "non-targeting")) {
	grna_target_df[1,"grna_target"] <- "non-targeting"
	warning("No NTs present. `grna_target_df` modified to have an NT")
}

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:", sum(grna_target_df$grna_target != "non-targeting"),
   "targeting and", sum(grna_target_df$grna_target == "non-targeting"),
  "non-targeting.\n")


# load extra covariates if they exist; keep at sceptre default otherwise
if(file.exists(covariates_fp <- file.path(input_dir, "covariate_data_frame.csv"))) {
	extra_covariates <- read.csv(covariates_fp, row.names = 1)
	cat("  extra_covariates loaded from file.\n")
} else {
	extra_covariates <- data.frame()
	cat("  no existing extra_covariates detected.\n")
}


# Load other sceptre params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load assign_grnas() formula if it exists else set it with pre-determined formula
if(file.exists(assign_grnas_fmla_fp <- file.path(input_dir, "assign_grnas_formula.rds"))) {
	assign_grnas_fmla <- readRDS(assign_grnas_fmla_fp) |> as.formula()
} else {
	# construct it
	fmla_lookup <- list(
		replogle = ~ 1 + log(grna_n_nonzero+1) + log(grna_n_umis+1),
		gasperini = ~ prep_batch + log(grna_n_nonzero+1) + log(grna_n_umis+1)
			    )
	assign_grnas_fmla <- fmla_lookup[[dataset_name]]

}
cat("  Using assign_grnas(..., formula_object =", as.character(assign_grnas_fmla), "\b).\n")

use_simulation_moi <- grepl("simulat", dataset_id, ignore.case = TRUE)
if(use_simulation_moi) {
	moi <- "high"
} else {
moi_lookup <- list(
  gasperini = "high",
  replogle = "low"
)
if(dataset_name %in% names(moi_lookup)) {
        moi <- moi_lookup[[dataset_name]]
} else {
        stop("MOI lookup missing ", dataset_name, ".")
}
}
cat("  MOI = ", moi, "\b.\n")

# Create sceptre object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# Assign gRNAs using mixture method with log-transformed covariates
cat("Assigning gRNAs using mixture method...\n")
sceptre_object <- assign_grnas(sceptre_object, method = "mixture", formula_object = assign_grnas_fmla)

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
