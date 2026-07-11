#!/usr/bin/env Rscript
# build_computational_replogle.R
# ONE driver for ALL Replogle (low-MOI) computational-benchmark datasets. The
# param table's `pipeline_only` flag selects, per dataset:
#   FALSE -> full subset export (for in-memory methods) + cells_to_remove (pipeline)
#   TRUE  -> lightweight pipeline-only inputs (discovery_pairs + cells_to_remove,
#            NO matrices / no ODM read)
# Replaces legacy make_computational_replogle-rd7.R.

library(tidyverse)
rm(list = ls())
source("~/.Rprofile")

script_dir <- dirname(sys.frame(1)$ofile)
lib <- file.path(script_dir, "lib")
for (f in c("io.R", "regime.R", "select.R", "validate.R", "assemble.R")) source(file.path(lib, f))

# ---- config ---------------------------------------------------------------
dataset         <- dataset_replogle
gene_qc_thresh  <- 7
METHODS_TO_SKIP <- character(0)          # e.g. c("frperturb", "mixscale")

# One param table for EVERY computational dataset (ngenes in {500, 1000}). `pipeline_only`:
#   FALSE = 6 in-memory datasets (full export): ngenes {500,1000} x ntargets {200,400,800}
#           (all 3 methods run; sceptre/mixscale/frperturb).
#   TRUE  = 4 pipeline-only datasets: ngenes {500,1000} x ntargets {1600,2383}
#           (nt1600 is pipeline-only FOR NOW -- in-memory nt1600 deferred to a later round;
#            flip it to FALSE + rebuild to add the matrices. 2383 = all targets, the full screen.)
# max_num_cells = Inf everywhere -> no downsampling (full cell union per target set).
# Dimension seeding (below): cells at a given ntargets are IDENTICAL across ngenes (clean ngenes axis).
# Cell counts (aligned sweep): 200->37k, 400->62k, 800->111k, 1600->215k, 2383(all)->~308k.
dataset_params <- rbind(
  expand.grid(num_genes = c(500, 1000), num_targets = c(200, 400, 800),  pipeline_only = FALSE),
  expand.grid(num_genes = c(500, 1000), num_targets = c(1600, 2383),     pipeline_only = TRUE)
) |>
  mutate(max_num_cells = Inf,
         dataset_name = paste0(dataset$name, "_max_comp_ngenes=", num_genes,
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

# ---- scp command to sync exactly these datasets to the cluster ------------
# Printed so it can be pasted straight into a terminal after every build.
# $LOCAL_BENCHMARKING_DIR / $REMOTE_BENCHMARKING_DIR are expanded by the shell.
.ds_names <- dataset_params$dataset_name
.ds_arg   <- if (length(.ds_names) > 1) paste0("{", paste(.ds_names, collapse = ","), "}") else .ds_names
cat("\n# ---- copy these datasets to the cluster (paste into a terminal) ----\n")
cat(sprintf(
  "scp -r ${LOCAL_BENCHMARKING_DIR}/association/computational/input_data/%s ${REMOTE_BENCHMARKING_DIR}/association/computational/input_data\n\n",
  .ds_arg))
