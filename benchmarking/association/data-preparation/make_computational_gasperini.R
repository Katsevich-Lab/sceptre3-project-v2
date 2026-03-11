#!/usr/bin/env Rscript
# make_computational_gasperini.R
# Script to create computational benchmarking datasets for gasperini data
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
  "guide_assignment/input_data/gasperini/sceptre-pipeline"
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/gasperini/sceptre-pipeline"
)

# Load data
cat("Loading Gasperini data...\n")
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
# dataset_params = data.frame(
#   num_genes =     c(30,   100,   500,    600),
#   num_targets =   c(50 ,  100,   200,    800),
#   max_num_cells = c(5000, 50000, 100000, 250000)
# ) |>
#   mutate(
#     dataset_name = paste0("gasperini_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k")
#   )


dataset_params = data.frame(
  num_genes =     c(100),
  num_targets =   c(100),
  max_num_cells = c(50000),
  # we sample randomly from genes with the
  gene_n_nonzero_p = c(0, .25, .5, .6, .75, .9, .97)
) |>
  mutate(
    dataset_name = paste0("gasperini_comp_ngenes=", num_genes,
                          "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_n_nonzero_p=", gene_n_nonzero_p)
  )

# cell idx is for the 50 one
gene_summary_stats = read_csv(file.path(path_to_data, "gene_summary_stats.csv"))
# top_genes = gene_summary_stats |> arrange(-gene_at_cell_idx_n_nonzero) |> head(100) |> pull(gene)


for(i in 1:nrow(dataset_params)) {
  dataset_name <- dataset_params$dataset_name[i]

  cat("\n=========================================\n")
  cat("Creating computational dataset:", dataset_name, "\n")
  cat("=========================================\n\n")

  num_genes = dataset_params$num_genes[i]
  all_valid_genes = gene_summary_stats |>
    filter(gene_n_nonzero >= quantile(gene_n_nonzero, dataset_params$gene_n_nonzero_p[i])) |>
    pull(gene)
  if(length(all_valid_genes) < num_genes) {
    stop('Not enough genes for run ', i)
  } else {
    set.seed(i + 123)
    genes = sample(all_valid_genes, num_genes)
  }
  
  make_computational_gasperini(
      dataset_name = dataset_name,
      response_odm=response_odm,
      grna_odm=grna_odm,
      cell_covariates=cell_covariates,
      scep_assn_mat=scep_assn_mat,
      grna_target_df=grna_target_df,
      genes=genes,
      num_targets  = dataset_params$num_targets[i],
      max_num_cells  = dataset_params$max_num_cells[i],
      nt_name = "non-targeting",
      force_nt_inclusion = FALSE,
      fr_perturb_concat_string = "@",
      random_seed = i
  )
}

cat("\n=========================================\n")
cat("All computational datasets created!\n")
cat("=========================================\n")

# scp -r gasperini_comp_ngenes=100_ntargets=100_ncells=50k_n_nonzero_p={0,0.25,0.5,0.6,0.75,0.9,0.97} hpc3:...data/projects/sceptre3/benchmarking/association/computational/input_data