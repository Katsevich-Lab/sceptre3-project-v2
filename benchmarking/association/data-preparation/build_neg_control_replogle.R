#!/usr/bin/env Rscript
# build_neg_control_replogle.R
# Thin driver: Replogle (low-MOI) negative-control dataset on the new lib.
# Replaces legacy make_neg_control_replogle-rd7.R.

library(tidyverse)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) {
  source(file.path(lib, f))
}

# ---- parameters -----------------------------------------------------------
dataset        <- dataset_replogle
num_genes      <- 1000
gene_qc_thresh <- 7
random_seed    <- 54654
dataset_name   <- paste0(dataset$name, "_max_neg_ctrl_ngenes=", num_genes,
                         "_gene_thresh=", gene_qc_thresh)

# ---- load -----------------------------------------------------------------
cat("Loading", dataset$name, "source data...\n")
src <- read_source_data(dataset$name)
stopifnot(src$low_moi == dataset$moi$low_moi)   # regime sanity check

al <- load_assignment(dataset$name, dataset$assignment_method, src$cell_covariates,
                      excluded_required = dataset$excluded_required)

genes_qc <- compute_genes_passing_qc(src$gene_summary_stats, al$assn_raw,
                                     src$grna_target_df, gene_qc_thresh)

# ---- select + build -------------------------------------------------------
spec <- select_neg(al$assn, src$grna_target_df, genes_qc, num_genes, dataset,
                   nt_name = "non-targeting", random_seed = random_seed)

output_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                         "association/neg-control/input_data")

build_dataset(spec, dataset, src, assn = al$assn, excluded_idx = al$excluded_idx,
              dataset_name = dataset_name, output_root = output_root,
              random_seed = random_seed, concat_string = "@")
