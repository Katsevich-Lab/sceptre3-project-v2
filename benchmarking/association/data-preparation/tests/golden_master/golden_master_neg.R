#!/usr/bin/env Rscript
# golden_master_neg.R  -- characterization test for the NEGATIVE-CONTROL builders.
# Legacy vs new lib on identical fake data + assignment; diffs every output.
# Shared plumbing is in golden_master_helpers.R.
#
# EXPECTED (intended) diffs -- everything else must read [same]:
#   Replogle : sceptre/cell_covariates.csv (+batch), frperturb_covariates (+batch),
#              mixscale/mixscale_nt_targets.rds (ONLY LEGACY -- file dropped)
#   Gasperini: sceptre/discovery_pairs.rds (new character vs legacy factor)

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

OUT <- make_out("golden_master_neg")
NUM_GENES <- 20; SEED <- 7

run_legacy <- function(regime, dn) {
  fake <- make_fake(regime); source_legacy(OUT)
  if (regime == "low") {
    make_neg_control_replogle(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes_passing_qc = fake$genes, num_genes = NUM_GENES, include_batch = FALSE, random_seed = SEED)
  } else {
    make_neg_control_gasperini(dataset_name = dn, response_odm = fake$resp, grna_odm = fake$grna,
      cell_covariates = fake$cov, scep_assn_mat = fake$assn, grna_target_df = fake$gtdf,
      genes_passing_qc = fake$genes, num_genes = NUM_GENES, nt_name = "non-targeting",
      fr_perturb_concat_string = "@", random_seed = SEED)
  }
}
run_new <- function(regime, dn) {
  fake <- make_fake(regime); source_new()
  ds <- if (regime == "low") dataset_replogle else dataset_gasperini
  spec <- select_neg(fake$assn, fake$gtdf, fake$genes, NUM_GENES, ds,
                     nt_name = "non-targeting", random_seed = SEED)
  src <- list(response_odm = fake$resp, cell_covariates = fake$cov)
  build_dataset(spec, ds, src, assn = fake$assn, excluded_idx = integer(0), dataset_name = dn,
                output_root = file.path(OUT, "new/association/neg-control/input_data"),
                random_seed = SEED, concat_string = "@")
}

for (regime in c("low", "high")) {
  cat(sprintf("\n=============== NEG-CONTROL (%s-MOI) ===============\n", regime))
  dn <- paste0("gm_neg_", regime)
  run_legacy(regime, dn); run_new(regime, dn)
  compare_dirs(file.path(OUT, "legacy/association/neg-control/input_data", dn),
               file.path(OUT, "new/association/neg-control/input_data", dn))
}
cat("\nExpected diffs only: Replogle cell_covariates/frperturb +batch, mixscale_nt_targets ONLY-LEGACY;",
    "Gasperini discovery_pairs character-vs-factor. Everything else [same].\n")
