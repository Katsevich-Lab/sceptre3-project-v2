#!/usr/bin/env Rscript
# golden_master_pos.R  -- characterization test for the POSITIVE-CONTROL builders.
# Legacy vs new lib on identical fake data + assignment. Shared plumbing in
# golden_master_helpers.R.
#
# NOT byte-clean: pos gains the always-on filter (legacy didn't run it). To keep
# this a REFACTOR check, the new selector runs with force_nt=TRUE here so the
# (separate, intended) pos-Gasperini NT-removal policy change doesn't swamp the
# signal. The production driver keeps force_nt=FALSE for high-MOI.
#
# EXPECTED (intended) diffs:
#   both regimes : sceptre/grna_matrix.rds + grna_target_data_frame.csv show
#                  "new SUBSET (filter)" (zero-cell guides dropped; bug #1);
#                  response_matrix may shrink (all-zero genes; bug #2)
#   Replogle     : cell_covariates.csv (top + sceptre) + frperturb gain `batch`
#   (no discovery_pairs -- sceptre builds positive_control_pairs itself)

rm(list = ls())
suppressMessages({ library(tidyverse); library(Matrix) })
source("~/.Rprofile")
pdf(NULL)   # swallow the ggplot the legacy pos-replogle function prints
.gm_dir <- local({
  a <- commandArgs(FALSE); fa <- grep("^--file=", a, value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  for (i in seq_len(sys.nframe())) { of <- sys.frame(i)$ofile; if (!is.null(of)) return(dirname(normalizePath(of))) }
  normalizePath(getwd())
})
source(file.path(.gm_dir, "golden_master_helpers.R"))

OUT <- make_out("golden_master_pos")
P <- list(num_targets = 3, max_cells = 40, min_cells = 1, seed = 7)

run_legacy <- function(regime, dn) {
  fake <- make_fake(regime, for_pos = TRUE); source_legacy(OUT)
  if (regime == "low") {
    make_pos_control_replogle_rd7(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes_passing_qc = fake$genes, num_targets = P$num_targets, max_num_cells = P$max_cells,
      only_consider_targets_with_this_many_cells = P$min_cells, random_seed = P$seed)
  } else {
    make_pos_control_gasperini(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes_passing_qc = fake$genes, num_targets = P$num_targets, max_num_cells = P$max_cells,
      only_consider_targets_with_this_many_cells = P$min_cells,
      fr_perturb_concat_string = "@", random_seed = P$seed)
  }
}
run_new <- function(regime, dn) {
  fake <- make_fake(regime, for_pos = TRUE); source_new()
  ds <- if (regime == "low") dataset_replogle else dataset_gasperini
  spec <- select_pos(fake$assn, fake$gtdf, fake$genes, P$num_targets, P$max_cells,
                     P$min_cells, ds, random_seed = P$seed, force_nt = TRUE)   # isolate refactor
  src <- list(response_odm = fake$resp, cell_covariates = fake$cov)
  build_dataset(spec, ds, src, assn = fake$assn, excluded_idx = integer(0), dataset_name = dn,
                output_root = file.path(OUT, "new/association/pos-control/input_data"),
                random_seed = P$seed, concat_string = "@")
}

for (regime in c("low", "high")) {
  cat(sprintf("\n=============== POSITIVE-CONTROL (%s-MOI, force_nt=TRUE) ===============\n", regime))
  dn <- paste0("gm_pos_", regime)
  run_legacy(regime, dn); run_new(regime, dn)
  compare_dirs(file.path(OUT, "legacy/association/pos-control/input_data", dn),
               file.path(OUT, "new/association/pos-control/input_data", dn))
}
cat("\nExpected: grna_matrix/grna_target_data_frame 'new SUBSET (filter)'; response may shrink;",
    "Replogle cell_covariates+frperturb +batch. Different cell count / targets = REAL regression.\n")
