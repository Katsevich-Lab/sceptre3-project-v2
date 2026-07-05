#!/usr/bin/env Rscript
# golden_master_comp.R  -- characterization test for the COMPUTATIONAL builders.
# Legacy vs new lib on identical fake data + assignment. Shared plumbing in
# golden_master_helpers.R.
#
# NOT byte-clean: comp gains the always-on filter (legacy didn't run it). The
# test confirms SELECTION matches legacy (same cells/targets/genes) so the only
# diffs are the intended ones. A different CELL COUNT or target set = regression.
#
# EXPECTED (intended) diffs:
#   both regimes : sceptre/grna_matrix.rds + grna_target_data_frame.csv show
#                  "new SUBSET (filter)" (zero-cell guides dropped; bug #1);
#                  response_matrix may shrink (all-zero genes; bug #2);
#                  discovery_pairs may drop now-empty targets (bug #1)
#   Replogle     : cell_covariates.csv (top + sceptre) + frperturb gain `batch`

rm(list = ls())
suppressMessages({ library(tidyverse); library(Matrix) })
source("~/.Rprofile")
.gm_dir <- local({
  a <- commandArgs(FALSE); fa <- grep("^--file=", a, value = TRUE)
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  for (i in seq_len(sys.nframe())) { of <- sys.frame(i)$ofile; if (!is.null(of)) return(dirname(normalizePath(of))) }
  normalizePath(getwd())
})
source(file.path(.gm_dir, "golden_master_helpers.R"))

OUT <- make_out("golden_master_comp")
P <- list(num_targets = 3, num_genes = 15, max_cells = 40, seed = 7)

run_legacy <- function(regime, dn) {
  fake <- make_fake(regime); source_legacy(OUT)
  set.seed(P$seed); genes <- sample(fake$genes, P$num_genes)   # legacy gene-sampling seed
  if (regime == "low") {
    make_computational_replogle(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes = genes, num_targets = P$num_targets, max_num_cells = P$max_cells, random_seed = P$seed)
  } else {
    make_computational_gasperini(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes = genes, num_targets = P$num_targets, max_num_cells = P$max_cells,
      force_nt_inclusion = FALSE, random_seed = P$seed, methods_to_skip = "frperturb")
  }
}
run_new <- function(regime, dn) {
  fake <- make_fake(regime); source_new()
  ds <- if (regime == "low") dataset_replogle else dataset_gasperini
  set.seed(P$seed); genes <- sample(fake$genes, P$num_genes)
  spec <- select_comp(fake$assn, fake$resp, fake$gtdf, genes, P$num_targets, P$max_cells,
                      ds, random_seed = P$seed)
  src <- list(response_odm = fake$resp, cell_covariates = fake$cov)
  skip <- if (regime == "high") "frperturb" else character(0)
  build_dataset(spec, ds, src, assn = fake$assn, excluded_idx = integer(0), dataset_name = dn,
                output_root = file.path(OUT, "new/association/computational/input_data"),
                random_seed = P$seed, concat_string = "@", methods_to_skip = skip)
}

for (regime in c("low", "high")) {
  cat(sprintf("\n=============== COMPUTATIONAL (%s-MOI) ===============\n", regime))
  dn <- paste0("gm_comp_", regime)
  run_legacy(regime, dn); run_new(regime, dn)
  compare_dirs(file.path(OUT, "legacy/association/computational/input_data", dn),
               file.path(OUT, "new/association/computational/input_data", dn))
}
cat("\nExpected: grna_matrix/grna_target_data_frame 'new SUBSET (filter)'; response/discovery may shrink;",
    "Replogle cell_covariates+frperturb +batch. Different cell count / targets = REAL regression.\n")
