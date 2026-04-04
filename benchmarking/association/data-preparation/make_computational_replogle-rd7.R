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

# getting gene QC parameters
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

## datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this is the first part
# dataset_params = data.frame(
#   num_genes =     rep(c(250, 375, 500), each = 3),
#   num_targets =   200,
#   max_num_cells = rep(c(50000, 75000, 98900), times = 3)
# ) |>
#   mutate(
#     # dataset_name = "replogle-test"
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
#   )
# 
# dataset_params = data.frame(
#   num_genes =     rep(c(222, 333, 444), each = 3),
#   num_targets =   225,
#   max_num_cells = rep(c(50000, 75000, 98900), times = 3)
# ) |>
#   mutate(
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
#   )

# dataset_params = data.frame(
#   num_genes =     rep(c(222, 333, 444), each = 3),
#   num_targets =   225,
#   max_num_cells = rep(c(50000, 75000, 98900), times = 3)
# ) |>
#   mutate(
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
#   )



# only varying ntargets. leaving cells alone
# dataset_params = data.frame(
#   num_genes =     1000,
#   num_targets =   c(50, 75, 100, 200),
#   max_num_cells = Inf
# ) |>
#   mutate(
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_gene_thresh=", gene_qc_thresh)
#   )

# ## this one is for scep-pipe only
# dataset_params = data.frame(
#   num_genes =     c(1000),
#   num_targets =   1600,#c(400,800,1600),
#   max_num_cells = Inf
# ) |>
#   mutate(
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_gene_thresh=", gene_qc_thresh)
#   )
# METHODS_TO_SKIP = c("mixscale", "frperturb")

gene_summary_stats = read_csv(file.path(path_to_data, "gene_summary_stats.csv"))

genes_passing_qc <- gene_summary_stats |>
  filter(gene_n_nonzero * moi_observed / total_num_guides * avg_num_guides_per_target >= gene_qc_thresh) |>
  pull(gene)

cat(length(genes_passing_qc), " genes (out of ", nrow(gene_summary_stats), ") pass gene QC.\n", sep="")
# paste0("rm -r ", strsplit(a, "\n")[[1]]) |> cat(sep="\n")
# top_genes = replogle_gene_summary_stats |> arrange(-gene_at_cell_idx_n_nonzero) |> head(100) |> pull(gene)

for(i in 1:nrow(dataset_params)) {
  dataset_name <- dataset_params$dataset_name[i]

  cat("\n=========================================\n")
  cat("Creating computational dataset:", dataset_name, "\n")
  cat("=========================================\n\n")

  num_genes = dataset_params$num_genes[i]
  # all_valid_genes = replogle_gene_summary_stats |>
  #   filter(gene_n_nonzero >= quantile(gene_n_nonzero, dataset_params$gene_n_nonzero_p[i])) |>
  #   pull(gene)
  all_valid_genes = genes_passing_qc
  if(length(all_valid_genes) < num_genes) {
    stop('Not enough genes for run ', i)
  } else {
    set.seed(i + 123)
    genes = sample(all_valid_genes, num_genes)
  }

  func_args = list(
    dataset_name = dataset_name,
    response_odm = response_odm,
    grna_odm = grna_odm,
    cell_covariates = cell_covariates,
    scep_assn_mat = scep_assn_mat,
    grna_target_df = grna_target_df,
    genes = genes,
    num_targets  = dataset_params$num_targets[i],
    max_num_cells  = dataset_params$max_num_cells[i],
    random_seed = i
  )
  if(exists("METHODS_TO_SKIP")) {
    func_args[["methods_to_skip"]] = METHODS_TO_SKIP
    cat("Skipping", METHODS_TO_SKIP, "\n")
  }

  do.call(make_computational_replogle, func_args)



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
  # cat("scp -r ", sceptre_pipeline_dir, " hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/association/computational/input_data/", dataset_name, "\n", sep="")
}




# making just what scep-pipe needs for some bigger ones ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TO RUN sceptre-pipeline WE ACTUALLY ONLY NEED THE FOLLOWING:
# 1. cells_to_remove.csv [located in <dataset>/sceptre-pipline/]
# 2. discovery_pairs.rds [located in <dataset>/sceptre/]
#
# So we address this as follows:
# i.  we will use the same ntargets as a previous run, so we copy the appropriate cells_to_remove.csv
# ii. we increase ngenes, so we need to make a new disc pairs
#     - load the old one, get all unique targets
#     - sample new genes, make cartesian product, write to new dataset file

# completed_run = "replogle-rd7_comp_ngenes=1000_ntargets=1600_gene_thresh=7"
# scep_pipe_params <- data.frame(
#   num_genes = c(5000, 10000)
# ) |>
#   mutate(
#   dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                         "_ntargets=1600_gene_thresh=", gene_qc_thresh)
# )
# base_comp_path = file.path(
#   .get_config_path("LOCAL_BENCHMARKING_DIR"),
#   "association/computational/input_data"
# )
# old_disc_pairs = read_rds(file.path(base_comp_path, completed_run, "sceptre/discovery_pairs.rds"))
# targets_to_use = unique(old_disc_pairs$grna_target) # these do not change
# 
# 
# # idea: we exactly copy cells_to_remove.csv, because ntargets is fixed so our cells is fixed
# # But we increase ngenes, so we need to make new discovery_pairs.rds
# for(i in 1:nrow(scep_pipe_params)) {
#   set.seed(i+10)
#   dataset_name = scep_pipe_params$dataset_name[i]
# 
#   path_to_new_dataset = file.path(
#     base_comp_path, dataset_name
#   )
# 
#   path_to_new_scep = file.path(path_to_new_dataset, "sceptre")
#   path_to_new_pipe = file.path(path_to_new_dataset, "sceptre-pipeline")
#   
#   dir.create(path_to_new_scep, showWarnings = FALSE, recursive = TRUE) # also makes dataset dir
#   dir.create(path_to_new_pipe, showWarnings = FALSE, recursive = TRUE)
#   
#   file.copy(
#     from = file.path(base_comp_path, completed_run, "sceptre-pipeline/cells_to_remove.csv"),
#     to = file.path(path_to_new_pipe, "cells_to_remove.csv")
#   ) 
#   
#   new_genes = sample(genes_passing_qc, scep_pipe_params$num_genes[i])
#   new_disc_pairs = expand.grid(
#     grna_target = targets_to_use,
#     response_id = new_genes
#   )
#   stopifnot(nrow(new_disc_pairs) == 1600 * scep_pipe_params$num_genes[i])
#   write_rds(new_disc_pairs, file.path(path_to_new_scep, "discovery_pairs.rds"))
#   
# }
# 







# here's what to scp
paste0(
  "scp -r /home/josep/data/projects/sceptre3/benchmarking/association/computational/input_data/",
  dataset_params$dataset_name,
  " hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/association/computational/input_data"
) |>
  cat(sep="\n")

# method,dataset,cpus,memory
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



# this became run_repl_3_cells_2_pairs_extra_mixscale_config.csv
# dataset_params = data.frame(
#   num_genes =     rep(c(250, 375, 500), each = 3),
#   num_targets =   200,
#   max_num_cells = rep(c(50000, 75000, 98900), times = 3)
# ) |>
#   mutate(
#     # dataset_name = "replogle-test"
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
#   )
# 
# config_df = rbind(
#   data.frame(
#     method = "mixscale",
#     dataset = dataset_params$dataset_name[c(1:6, 9)],
#     cpus=1,
#     memory = c(100,100,110,100,100,120,148)
#   ),
#   data.frame(
#     method = "sceptre",
#     dataset = dataset_params$dataset_name[1:6],
#     cpus=1,
#     memory = c(15,20,25,15,20,30)
#   ),
#   data.frame(
#     method = "frperturb",
#     dataset = dataset_params$dataset_name[1:6],
#     cpus=1,
#     memory = c(20)
#   )
# ) 
# sum(config_df$memory)

# this becomes part of  run_repl_remaining_runs_200_or_225_targets_config.csv
# config_df_pt1 = rbind(
#   data.frame(
#     method = "mixscale",
#     dataset = dataset_params$dataset_name[7:8],
#     cpus=1,
#     memory = c(130)
#   ),
#   data.frame(
#     method = "sceptre",
#     dataset = dataset_params$dataset_name[7:9],
#     cpus=1,
#     memory = c(15)
#   ),
#   data.frame(
#     method = "frperturb",
#     dataset = dataset_params$dataset_name[7:9],
#     cpus=1,
#     memory = c(10)
#   )
# )
# sum(config_df_pt1$memory)


# dataset_params = data.frame(
#   num_genes =     rep(c(222, 333, 444), each = 3),
#   num_targets =   225,
#   max_num_cells = rep(c(50000, 75000, 98900), times = 3)
# ) |>
#   mutate(
#     dataset_name = paste0("replogle-rd7_comp_ngenes=", num_genes,
#                           "_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)
#   )
# # this became run_repl_3_cells_2_pairs_more_targets_config.csv
# config_df = rbind(
#   data.frame(
#     method = "mixscale",
#     dataset = dataset_params$dataset_name[1:6],
#     cpus=1,
#     memory = c(130,135,145,150,150,160)
#   ),
#   data.frame(
#     method = "sceptre",
#     dataset = dataset_params$dataset_name[1:6],
#     cpus=1,
#     memory = 15
#   ),
#   data.frame(
#     method = "frperturb",
#     dataset = dataset_params$dataset_name[1:6],
#     cpus=1,
#     memory = 10
#   )
# )
# sum(config_df$memory)

# this becomes the rest of run_repl_remaining_runs_200_or_225_targets_config.csv witj config_df_pt1
# config_df_pt2 = rbind(
#   data.frame(
#     method = "mixscale",
#     dataset = dataset_params$dataset_name[7:9],
#     cpus=1,
#     memory = c(140)
#   ),
#   data.frame(
#     method = "sceptre",
#     dataset = dataset_params$dataset_name[7:9],
#     cpus=1,
#     memory = c(15)
#   ),
#   data.frame(
#     method = "frperturb",
#     dataset = dataset_params$dataset_name[7:9],
#     cpus=1,
#     memory = c(10)
#   )
# )
# sum(config_df_pt2$memory)
# config_df = rbind(config_df_pt1, config_df_pt2)
# sum(config_df$memory)

# config_df_to_file(config_df)
# config_df_to_file(rbind(config_df, config_df2))
