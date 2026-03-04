#!/usr/bin/env Rscript
# make_computational_replogle-rd7.R
# Script to create computational benchmarking datasets for Replogle RD7 data
# Creates multiple dataset sizes to test computational scalability
# Uses modular functions from neg_control_functions.R and utils_data_prep.R

library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

# Source utility functions
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "utils_data_prep.R"))
source(file.path(script_dir, "neg_control_functions.R"))
source(file.path(script_dir, "computational_benchmarking_functions.R"))

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

response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
cell_covariates = scep@covariate_data_frame
scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
grna_target_df = scep@grna_target_data_frame

# Create datasets at multiple scales for computational benchmarking
dataset_params = data.frame(
  num_genes =     c(100,   500,    600),
  num_targets =   c(100,   200,    800), 
  max_num_cells = c(50000, 100000, 250000)
) |>
  mutate(
    dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
                          "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k")
  )

for(i in 1:nrow(dataset_params)) {
  dataset_name <- dataset_params$dataset_name[i]

  cat("\n=========================================\n")
  cat("Creating computational dataset:", dataset_name, "\n")
  cat("=========================================\n\n")

  make_computational_replogle(
    dataset_name = dataset_name,
    response_odm = response_odm,
    grna_odm = grna_odm,
    cell_covariates = cell_covariates,
    scep_assn_mat = scep_assn_mat,
    grna_target_df = grna_target_df,
    num_genes = dataset_params$num_genes[i],
    num_targets  = dataset_params$num_targets[i],
    max_num_cells  = dataset_params$max_num_cells[i],
    random_seed = i
  )
  
  
  
}

cat("\n=========================================\n")
cat("All computational datasets created!\n")
cat("=========================================\n")
