#!/usr/bin/env Rscript
# build_pos_control_gasperini.R
# Thin driver: Gasperini (high-MOI) positive-control (on-target) dataset.
# Replaces legacy make_pos_control_gasperini.R.
#
# NOTE: high-MOI => NT is NOT force-included. Legacy pos-Gasperini once
# force-included NT to boost Gasperini's small positive-control N; to re-enable
# that, pass force_nt = TRUE to select_pos() below.

library(tidyverse)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) source(file.path(lib, f))

# ---- config ---------------------------------------------------------------
dataset        <- dataset_gasperini
gene_qc_thresh <- 7
num_targets    <- 377            # the max available
max_num_cells  <- 100000
only_consider_targets_with_this_many_cells <- 100
random_seed    <- 8765623
dataset_name   <- paste0(dataset$name, "_t5_pos_ctrl_ntargets=", num_targets,
                         "_ncells=", max_num_cells / 1000, "k_gene_thresh=", gene_qc_thresh)

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
