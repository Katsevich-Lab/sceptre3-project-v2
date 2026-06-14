#!/usr/bin/env Rscript

# Dispatcher for "script_*" guide-assignment methods.
#
# Strips the "script_" prefix from the method name to find the matching variant
# file under bin/script/{variant}.R. Each variant must define:
#
#   assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
#                                    extra_covariates, formula, moi, cpus)
#
# returning a list:
#   list(
#     assignment_matrix = <sparse logical Matrix; rows = grnas, cols = cells,
#                          rownames + colnames set>,
#     extras            = <anything to save for inspection, or NULL>
#   )
#
# `assignment_matrix` is saved to assignment_matrix_script.rds (required).
# `extras` is saved to assignment_extras_script.rds when non-NULL (optional).

# Resolve own script path so we can locate the bin/script/ sibling directory
# regardless of cwd. commandArgs(trailingOnly = FALSE) is robust under Rscript.
self_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(self_arg) != 1) {
  stop("Could not determine script path via --file= argument.")
}
self_path <- sub("^--file=", "", self_arg)
bin_dir <- dirname(normalizePath(self_path))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: run_script.R <input_dir> <dataset_id> <method_name> <cpus>")
}
input_dir   <- args[1]
dataset_id  <- args[2]
method_name <- args[3]
cpus        <- as.integer(args[4])

# Determine variant from method name
if (!grepl("^script_", method_name)) {
  stop("Method name must start with 'script_'. Got: ", method_name)
}
variant <- sub("^script_", "", method_name)
variant_file <- file.path(bin_dir, "script", paste0(variant, ".R"))
if (!file.exists(variant_file)) {
  stop("Variant file not found: ", variant_file)
}

# Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Loading data from:", input_dir, "\n")
response_matrix <- readRDS(file.path(input_dir, "response_matrix.rds"))
grna_matrix     <- readRDS(file.path(input_dir, "grna_matrix.rds"))
grna_target_df  <- read.csv(file.path(input_dir, "grna_target_data_frame.csv"))
if (!any(grna_target_df$grna_target == "non-targeting")) {
  grna_target_df[1, "grna_target"] <- "non-targeting"
  warning("No NTs present. `grna_target_df` modified to have an NT.")
}

cat("Data loaded:\n")
cat("  Response matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  gRNA matrix:    ", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  gRNA targets:   ", sum(grna_target_df$grna_target != "non-targeting"),
    "targeting and", sum(grna_target_df$grna_target == "non-targeting"),
    "non-targeting.\n")

covariates_fp <- file.path(input_dir, "cell_covariates.csv")
if (file.exists(covariates_fp)) {
  extra_covariates <- read.csv(covariates_fp)
  cat("  extra_covariates loaded from file.\n")
} else {
  stop("covariates should be found at: ", covariates_fp)
}

formula_fp <- file.path(input_dir, "formula_object.rds")
if (file.exists(formula_fp)) {
  assign_grnas_fmla <- as.formula(readRDS(formula_fp))
} else {
  stop("formula_object.rds should be found at: ", formula_fp)
}
cat("  formula:", as.character(assign_grnas_fmla), "\n")

# MOI lookup mirrors run_sceptre.R
use_simulation_moi <- grepl("simulat", dataset_id, ignore.case = TRUE) |
  grepl("^sim", dataset_id, ignore.case = TRUE)
if (use_simulation_moi) {
  moi <- "high"
} else {
  DATASET_NAMES <- c("gasperini", "replogle")
  dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(name) {
    grepl(name, dataset_id, ignore.case = TRUE)
  })]
  if (length(dataset_name) != 1) {
    stop("Could not determine dataset from dataset_id: ", dataset_id)
  }
  moi_lookup <- list(gasperini = "high", replogle = "low")
  if (!(dataset_name %in% names(moi_lookup))) {
    stop("MOI lookup missing ", dataset_name, ".")
  }
  moi <- moi_lookup[[dataset_name]]
}
cat("  MOI:", moi, "\n")

# Dispatch to variant ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Sourcing variant:", variant_file, "\n")
source(variant_file)
if (!exists("assign_grnas_script", inherits = FALSE) ||
    !is.function(get("assign_grnas_script"))) {
  stop("Variant file did not define a function `assign_grnas_script`: ", variant_file)
}

cat("Running assign_grnas_script() with cpus =", cpus, "\n")
ret <- assign_grnas_script(
  response_matrix   = response_matrix,
  grna_matrix       = grna_matrix,
  grna_target_df    = grna_target_df,
  extra_covariates  = extra_covariates,
  formula           = assign_grnas_fmla,
  moi               = moi,
  cpus              = cpus
)

if (!is.list(ret) || is.null(ret$assignment_matrix)) {
  stop("assign_grnas_script() must return a list with `$assignment_matrix` ",
       "(and optionally `$extras`).")
}
assignment_matrix <- ret$assignment_matrix
extras            <- ret$extras

if (is.null(rownames(assignment_matrix)) || is.null(colnames(assignment_matrix))) {
  stop("assignment_matrix must have rownames (gRNAs) and colnames (cells).")
}

cat("Assignment matrix:", nrow(assignment_matrix), "gRNAs x",
    ncol(assignment_matrix), "cells; total assignments =",
    sum(assignment_matrix), "\n")

saveRDS(assignment_matrix, "assignment_matrix_script.rds")
cat("Wrote assignment_matrix_script.rds\n")

if (!is.null(extras)) {
  saveRDS(extras, "assignment_extras_script.rds")
  cat("Wrote assignment_extras_script.rds\n")
}
