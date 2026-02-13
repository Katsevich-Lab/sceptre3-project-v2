#!/usr/bin/env Rscript
# make_neg_control_replogle-rd7.R
# Script to create negative control datasets for Replogle RD7 data
# Uses modular functions from neg_control_functions.R and utils_data_prep.R

library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

# Source utility functions
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "utils_data_prep.R"))
source(file.path(script_dir, "neg_control_functions.R"))

# Set up paths
path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/replogle-rd7/sceptre-pipeline"
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/replogle-rd7/sceptre-pipeline"
)

# Load data
cat("Loading Replogle RD7 data...\n")
scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

# Run negative control data preparation
num_genes <- 5000
dataset_name <- paste0("replogle-rd7_neg_control_", num_genes, "genes")

cat("\nCreating negative control dataset:", dataset_name, "\n")
cat("=========================================\n\n")

make_neg_control_replogle(
  dataset_name = dataset_name,
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame,
  num_genes = num_genes,
  include_batch = FALSE,
  random_seed = 243535
)

cat("\n=========================================\n")
cat("Done!\n")


# Uncomment to test with smaller dataset:
# num_genes_test <- 100
# make_neg_control_replogle(
#   dataset_name = paste0("replogle-rd7_neg_control_", num_genes_test, "genes_TEST"),
#   response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
#   grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
#   cell_covariates = scep@covariate_data_frame,
#   scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
#   grna_target_df = scep@grna_target_data_frame,
#   num_genes = num_genes_test,
#   include_batch = FALSE,
#   random_seed = 243535
# )
