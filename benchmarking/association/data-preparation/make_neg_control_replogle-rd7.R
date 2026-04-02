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

response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
grna_target_df = scep@grna_target_data_frame
cell_covariates = scep@covariate_data_frame


num_genes <- 1000
gene_qc_thresh = 7
random_seed = 54654

dataset_name = paste0("replogle-rd7_neg_ctrl_ngenes=", num_genes, "_gene_thresh=", gene_qc_thresh)


moi_observed <- Matrix::colSums(scep_assn_mat) |> mean()

total_num_guides <- grna_target_df |>
  filter(! grna_target %in% c("nt_off_target", "unknown")) |>
  pull(grna_id) |>
  unique() |>
  length()

avg_num_guides_per_target <- grna_target_df |>
  filter(! grna_target %in% c("nt_off_target", "unknown")) |>
  group_by(grna_target) |>
  summarize(n=n()) |>
  pull(n) |>
  mean()

gene_summary_stats = read_csv(file.path(path_to_data, "gene_summary_stats.csv"), show_col_types = FALSE)

genes_passing_qc <- gene_summary_stats |>
  filter(gene_n_nonzero * moi_observed / total_num_guides * avg_num_guides_per_target >= gene_qc_thresh) |>
  pull(gene)

cat("There are", length(genes_passing_qc), "genes passing QC.\n")

cat("\nCreating negative control dataset:", dataset_name, "\n")
cat("=========================================\n\n")

make_neg_control_replogle(
  dataset_name = dataset_name,
  response_odm = response_odm,
  grna_odm = grna_odm,
  cell_covariates = cell_covariates,
  scep_assn_mat = scep_assn_mat,
  grna_target_df = grna_target_df,
  num_genes = num_genes,
  genes_passing_qc = genes_passing_qc,
  include_batch = FALSE,
  random_seed = random_seed
)

cat("\n=========================================\n")
cat("Done!\n")
