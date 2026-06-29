#!/usr/bin/env Rscript
# Generic sceptre runner for survey datasets (no dataset_id branching).
# Same input contract as the lab's benchmarking/association/{pos,neg}-control-
# pipeline/bin/run_sceptre.R, but uses a generic formula + MOI from a meta.json
# alongside the sceptre/ input directory.
#
# Inputs (in dataset_dir = <out>/sceptre/):
#   response_matrix.rds         genes x cells
#   grna_matrix.rds             guides x cells  (binary; our method's calls)
#   grna_target_data_frame.csv  grna_id, grna_target
#   cell_covariates.csv         per-cell covariates (the _full ones)
#   discovery_pairs.rds         (neg only)
#   formula_object.rds          (neg only; pos uses the generic formula here)
#   meta.json                   { kind: "pos"|"neg", moi: "high"|"low" }
#
# Outputs (in dataset_dir):
#   pos: association_on_target_sceptre.csv
#   neg: association_neg_control_sceptre.csv
#
# Usage: Rscript run_sceptre_survey.R <dataset_dir> <pos|neg>
suppressPackageStartupMessages({library(sceptre); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]; kind <- args[2]
stopifnot(kind %in% c("pos", "neg"))

# read meta (kind + moi)
meta_f <- file.path(dataset_dir, "meta.json")
moi <- if (file.exists(meta_f)) {
  m <- jsonlite::fromJSON(meta_f); if (!is.null(m$moi)) m$moi else "low"
} else "low"
cat(sprintf("[survey runner] kind=%s moi=%s\n", kind, moi))

resp <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
tdf  <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cov  <- read.csv(file.path(dataset_dir, "cell_covariates.csv"),       stringsAsFactors = FALSE)
cat(sprintf("  resp %dx%d  grna %dx%d  targets %d  cov %dx%d\n",
            nrow(resp), ncol(resp), nrow(grna), ncol(grna), nrow(tdf), nrow(cov), ncol(cov)))

# generic formula: log-transformed cell-level covariates that exist
fmla_str <- "~ log(response_n_umis_full + 1) + log(response_n_nonzero_full + 1) + log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1)"
if (kind == "neg" && file.exists(file.path(dataset_dir, "formula_object.rds"))) {
  fmla_str <- readRDS(file.path(dataset_dir, "formula_object.rds"))
}
fmla <- as.formula(fmla_str)
cat("  formula:", deparse(fmla), "\n")

scep <- import_data(response_matrix = resp, grna_matrix = grna,
                    grna_target_data_frame = tdf, moi = moi, extra_covariates = cov)

if (kind == "pos") {
  pcp <- construct_positive_control_pairs(scep)
  cat(sprintf("  positive control pairs: %d\n", nrow(pcp)))
  scep <- set_analysis_parameters(scep, positive_control_pairs = pcp, formula_object = fmla)
  scep <- assign_grnas(scep, method = "thresholding", threshold = 1)
  # run_qc is permissive but sceptre's NT-cell update bombs on tiny NT sets;
  # skip it and let run_power_check apply default filtering.
  scep <- tryCatch(run_qc(scep, n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0,
                          response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
                          p_mito_threshold = 1),
                   error = function(e) { cat("  (run_qc skipped:", conditionMessage(e), ")\n"); scep })
  scep <- run_power_check(scep, output_amount = 2)
  res <- get_result(scep, "run_power_check")
  out <- file.path(dataset_dir, "..", "association_on_target_sceptre.csv")
  write.csv(res, out, row.names = FALSE)
} else {
  disc <- readRDS(file.path(dataset_dir, "discovery_pairs.rds"))
  scep <- set_analysis_parameters(scep, discovery_pairs = disc, formula_object = fmla,
                                  control_group = "complement")
  scep <- assign_grnas(scep, method = "thresholding", threshold = 1)
  scep <- tryCatch(run_qc(scep, n_nonzero_trt_thresh = 0, n_nonzero_cntrl_thresh = 0,
                          response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
                          p_mito_threshold = 1),
                   error = function(e) { cat("  (run_qc skipped:", conditionMessage(e), ")\n"); scep })
  ncpus <- as.integer(Sys.getenv("NCPUS", "0"))
  if (ncpus >= 2) {
    scep <- run_discovery_analysis(scep, parallel = TRUE, n_processors = ncpus, output_amount = 2)
  } else {
    scep <- run_discovery_analysis(scep, parallel = FALSE, output_amount = 2)
  }
  res <- get_result(scep, "run_discovery_analysis")
  out <- file.path(dataset_dir, "..", "association_neg_control_sceptre.csv")
  write.csv(res, out, row.names = FALSE)
}
cat(sprintf("[survey runner] wrote %s (%d rows)\n", out, nrow(res)))
