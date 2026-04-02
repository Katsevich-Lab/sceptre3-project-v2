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

gene_qc_thresh <- 7

gene_summary_stats = read_csv(file.path(path_to_data, "gene_summary_stats.csv"), show_col_types = FALSE)

genes_passing_qc <- gene_summary_stats |>
  filter(gene_n_nonzero * moi_observed / total_num_guides * avg_num_guides_per_target >= gene_qc_thresh) |>
  pull(gene)

cat(length(genes_passing_qc), " genes (out of ", nrow(gene_summary_stats), ") pass gene QC.\n", sep="")


dataset_params = data.frame(
  num_genes =     rep(c(250, 375, 500), each = 3),
  num_targets =   200,
  max_num_cells = rep(c(50000, 75000, 100000), times = 3)
) |>
  mutate(
    dataset_name = paste0("gasperini_comp_ngenes=", num_genes,
                          "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
  )



for(i in 1:nrow(dataset_params)) {
  dataset_name <- dataset_params$dataset_name[i]

  cat("\n=========================================\n")
  cat("Creating computational dataset:", dataset_name, "\n")
  cat("=========================================\n\n")

  num_genes = dataset_params$num_genes[i]
  curr_genes = sample(genes_passing_qc, num_genes)
  
  make_computational_gasperini(
      dataset_name = dataset_name,
      response_odm=response_odm,
      grna_odm=grna_odm,
      cell_covariates=cell_covariates,
      scep_assn_mat=scep_assn_mat,
      grna_target_df=grna_target_df,
      genes=curr_genes,
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


## making sceptre-pipeline/cells_to_remove.csv -----------------------------------------------
for(i in 1:nrow(dataset_params)) {
  dataset_name <- dataset_params$dataset_name[i]
  dataset_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/computational/input_data", dataset_name
  )
  
  cell_info <- read.csv(file.path(dataset_fp, "cell_info.csv"))
  stopifnot(!any(duplicated(cell_info$cell_idx)))
  stopifnot(max(cell_info$cell_idx) <= ncol(scep_assn_mat))
  
  sceptre_pipeline_dir <- file.path(dataset_fp, "sceptre-pipeline")
  dir.create(sceptre_pipeline_dir, showWarnings = FALSE, recursive = TRUE)
  
  cells_to_remove <- setdiff(seq_len(ncol(scep_assn_mat)), cell_info$cell_idx)
  writeLines(as.character(cells_to_remove), file.path(sceptre_pipeline_dir, "cells_to_remove.csv"))
  cat("scp -r ", sceptre_pipeline_dir, " hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/association/computational/input_data/", dataset_name, "\n", sep="")
}





# scp'ing all of them to the server ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# for(i in 1:nrow(dataset_params)) {
#   dataset_name <- dataset_params$dataset_name[i]
#   dataset_fp <- file.path(
#     .get_config_path("LOCAL_BENCHMARKING_DIR"),
#     "association/computational/input_data", dataset_name
#   )
#   cat("scp -r ", dataset_fp, " hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/association/computational/input_data/", dataset_name, "\n", sep="")
# }

# make sceptre-pipeline datasets_config.csv
data.frame(
  dataset_name = dataset_params$dataset_name,
  npairs = rep(c(50000, 75000, 100000), each = 3) |> as.character()
) |>
  apply(1, paste0, collapse=",") |>
  cat(sep="\n")



config_df_to_file <- function(config_df) {
  header = "method,dataset,cpus,memory"
  cat(header, "\n")
  config_df|>
    mutate(
      memory = paste0(memory, ".GB")
    ) |>
    apply(1, paste0, collapse=",") |>
    
    cat(sep="\n")
}

config_df = rbind(
  data.frame(
    method = "sceptre",
    dataset = dataset_params$dataset_name,
    cpus=1,
    memory = 15
  ),
  data.frame(
    method = "frperturb",
    dataset = dataset_params$dataset_name,
    cpus=1,
    memory = 15
  )
)
sum(config_df$memory)

config_df_to_file(config_df)





