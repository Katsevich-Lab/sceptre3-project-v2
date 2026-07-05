#!/usr/bin/env Rscript
# build_computational_gasperini.R
# ONE driver for ALL Gasperini (high-MOI) computational-benchmark datasets. The
# param table's `pipeline_only` flag selects, per dataset:
#   FALSE -> full subset export (for in-memory methods) + cells_to_remove (pipeline)
#   TRUE  -> lightweight pipeline-only inputs (discovery_pairs + cells_to_remove,
#            NO matrices / no ODM read)
# Replaces legacy make_computational_gasperini.R.
#
# High-MOI: writes sceptre + frperturb (no mixscale). complement control => every
# sceptre pair runs over ALL kept cells, so sceptre-direct is the time limiter
# (not memory); verify its runtime at the biggest in-memory point before trusting
# the full grid. max_num_cells = Inf (no downsampling) -- same call as Replogle.

library(tidyverse)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) source(file.path(lib, f))

# ---- config ---------------------------------------------------------------
dataset         <- dataset_gasperini
gene_qc_thresh  <- 7
METHODS_TO_SKIP <- character(0)          # high-MOI writers = sceptre + frperturb

# One param table for EVERY computational dataset. `pipeline_only`:
#   FALSE = 6 in-memory datasets (full export): ngenes {250,500} x ntargets {100,200,400}
#   TRUE  = 4 pipeline-only extension datasets: ngenes {250,500} x ntargets {800,1600}
# max_num_cells = Inf everywhere -> no downsampling (full cell union per target set).
# Dimension-based seeding (below) -> cells at a given ntargets are identical across ngenes.
dataset_params <- rbind(
  expand.grid(num_genes = c(250, 500), num_targets = c(100, 200, 400), pipeline_only = FALSE),
  expand.grid(num_genes = c(250, 500), num_targets = c(800, 1600),     pipeline_only = TRUE)
) |>
  mutate(max_num_cells = Inf,
         dataset_name = paste0(dataset$name, "_t5_comp_ngenes=", num_genes,
                               "_ntargets=", num_targets, "_gene_thresh=", gene_qc_thresh))

# ---- load -----------------------------------------------------------------
src <- read_source_data(dataset$name)
stopifnot(src$low_moi == dataset$moi$low_moi)
al  <- load_assignment(dataset$name, dataset$assignment_method, src$cell_covariates,
                       excluded_required = dataset$excluded_required)
genes_qc <- compute_genes_passing_qc(src$gene_summary_stats, al$assn_raw,
                                     src$grna_target_df, gene_qc_thresh)

output_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                         "association/computational/input_data")

# ---- build each dataset ---------------------------------------------------
for (i in seq_len(nrow(dataset_params))) {
  ng <- dataset_params$num_genes[i]
  nt <- dataset_params$num_targets[i]
  dn <- dataset_params$dataset_name[i]
  po <- dataset_params$pipeline_only[i]
  cat("\n================ ", dn, if (po) " [pipeline-only]" else "", " ================\n")

  # Seed by DIMENSION, not row index, so the cell set is IDENTICAL across ngenes
  # at a fixed ntargets (clean ngenes axis; full & extension share cells per ntargets).
  set.seed(700L + ng)                                 # genes: same set per ngenes level
  curr_genes <- sample(genes_qc, ng)

  spec <- select_comp(al$assn, src$response_odm, src$grna_target_df, curr_genes,
                      num_targets = nt, max_num_cells = dataset_params$max_num_cells[i],
                      dataset = dataset, random_seed = 30000L + nt)   # targets/cells: same per ntargets

  if (po) {
    write_pipeline_inputs(spec, ncol(al$assn), dn, output_root)              # lightweight: no matrices
  } else {
    res <- build_dataset(spec, dataset, src, assn = al$assn, excluded_idx = al$excluded_idx,
                         dataset_name = dn, output_root = output_root,
                         random_seed = 30000L + nt, concat_string = "@", methods_to_skip = METHODS_TO_SKIP)
    write_cells_to_remove(res$path, res$cell_idx, res$n_cells_total)         # full export
  }
}

cat("\nAll computational datasets created.\n")
