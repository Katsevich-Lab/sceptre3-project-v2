#!/usr/bin/env Rscript
# build_pos_control_replogle.R
# Thin driver: Replogle (low-MOI) positive-control (on-target) dataset.
# Replaces legacy make_pos-control_replogle-rd7.R.

library(tidyverse)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) source(file.path(lib, f))

# ---- config ---------------------------------------------------------------
dataset        <- dataset_replogle
gene_qc_thresh <- 7
num_targets    <- 750
max_num_cells  <- Inf # 104235  # this is the exact number that come up for this seed, so actually no downsampling 
only_consider_targets_with_this_many_cells <- 25
random_seed    <- 876568
dataset_name   <- paste0(dataset$name, "_max_pos_ctrl_ntargets=", num_targets,
                         "_gene_thresh=", gene_qc_thresh)   # no ncells: full data, induced by ntargets

# ---- load -----------------------------------------------------------------
src <- read_source_data(dataset$name)
stopifnot(src$low_moi == dataset$moi$low_moi)
al  <- load_assignment(dataset$name, dataset$assignment_method, src$cell_covariates,
                       excluded_required = dataset$excluded_required)
genes_qc <- compute_genes_passing_qc(src$gene_summary_stats, al$assn_raw,
                                     src$grna_target_df, gene_qc_thresh)

# ---- select + build -------------------------------------------------------
spec <- select_pos(al$assn, src$grna_target_df, genes_qc, num_targets, max_num_cells,
                   only_consider_targets_with_this_many_cells, dataset,
                   nt_name = "non-targeting", random_seed = random_seed)

output_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                         "association/pos-control/input_data")

build_dataset(spec, dataset, src, assn = al$assn, excluded_idx = al$excluded_idx,
              dataset_name = dataset_name, output_root = output_root,
              random_seed = random_seed, concat_string = "@")
