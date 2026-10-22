#!/usr/bin/env Rscript
# sceptre positive-control run on a full-Gasperini dataset from build-dataset.R,
# for one assignment method and one cell-QC arm.
#
#   OUT=<LOCAL_BENCHMARKING_DIR>/association/assign-pos-control/outputs/<dataset>/<method>_<qc>
#   mkdir -p $OUT && SIF=<pipeline-images>/sceptre.sif \
#     sbatch -c 2 -t 1:00:00 -o $OUT/sceptre-%j.out \
#     ../../guide-assignment/grna-simulator/run-in-image.sh run-sceptre.R <dataset> <method> <qc>
#
#   <dataset> gasperini_allcells        377 TSS pairs (target vs its own gene), built
#                                       by construct_positive_control_pairs(), tested
#                                       with run_power_check()
#             gasperini_allcells_globin 130 globin control pairs from
#                                       shared/discovery_pairs.rds, tested with
#                                       run_discovery_analysis() (same test; the
#                                       output adds BH-adjusted columns)
#   <method>  fishash | fishashplus | thresh3 | thresh5 | thresh10
#   <qc>      none    : no cell QC (as pos-control-pipeline/bin/run_sceptre.R)
#             default : sceptre's default cell-QC rule, applied in build-dataset.R
#                       to the FULL-transcriptome covariates; passed here as
#                       additional_cells_to_remove. sceptre's own cell filters are
#                       set to remove nothing, because they would use totals over
#                       only the 377 genes in this response matrix.
#
# Settings otherwise follow pos-control-pipeline/bin/run_sceptre.R (Gasperini):
# high MOI, *_full covariates in the formula, two-sided, guides pooled per
# target, no pairwise QC. Every argument is written out, defaults included.
#
# Writes <OUT>/association_on_target.csv: the sceptre result (power check or
# discovery analysis) plus n_trt_cells (cells in use assigned any guide of the
# target).

suppressPackageStartupMessages({ library(Matrix); library(sceptre) })
source("~/.Rprofile")                                    # .get_config_path()

args <- commandArgs(trailingOnly = TRUE)
DATASETS <- c(gasperini_allcells = 377L, gasperini_allcells_globin = 130L)   # dataset = n pairs
METHODS  <- c("fishash", "fishashplus", "thresh3", "thresh5", "thresh10")
if (length(args) != 3 || !args[1] %in% names(DATASETS) || !args[2] %in% METHODS ||
    !args[3] %in% c("none", "default")) {
  stop("usage: run-sceptre.R <", paste(names(DATASETS), collapse = "|"), "> <",
       paste(METHODS, collapse = "|"), "> <none|default>")
}
dataset <- args[1]; method <- args[2]; qc <- args[3]
use_discovery <- dataset == "gasperini_allcells_globin"

bench   <- .get_config_path("LOCAL_BENCHMARKING_DIR")
in_root <- file.path(bench, "association", "assign-pos-control", "input_data", dataset)
out_dir <- file.path(bench, "association", "assign-pos-control", "outputs", dataset, paste0(method, "_", qc))
out_csv <- file.path(out_dir, "association_on_target.csv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(out_csv)) stop("result already exists, not overwriting: ", out_csv)
cat("dataset:", dataset, "| method:", method, "| qc:", qc, "\nout:", out_dir, "\n\n")

N_CELLS <- 207324L; N_PAIRS <- DATASETS[[dataset]]

# ---- inputs ------------------------------------------------------------------
response_matrix <- readRDS(file.path(in_root, "shared", "response_matrix.rds"))
grna_matrix     <- readRDS(file.path(in_root, "grna", paste0(method, ".rds")))
grna_target_df  <- read.csv(file.path(in_root, "shared", "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(in_root, "shared", "cell_covariates.csv"), stringsAsFactors = FALSE)
qc_excl         <- readRDS(file.path(in_root, "shared", "qc_excluded_cells.rds"))

stopifnot(ncol(response_matrix) == N_CELLS, ncol(grna_matrix) == N_CELLS,
          identical(colnames(response_matrix), cell_covariates$cell_id),
          identical(colnames(grna_matrix), cell_covariates$cell_id),
          identical(rownames(grna_matrix), grna_target_df$grna_id),
          all(grna_matrix@x == 1))
cat("response:", nrow(response_matrix), "x", ncol(response_matrix),
    "| grna:", nrow(grna_matrix), "x", ncol(grna_matrix), ";", sum(grna_matrix), "assignments\n")

covariates_to_use <- c("response_n_nonzero_full", "response_n_umis_full",
                       "grna_n_nonzero_full", "grna_n_umis_full", "prep_batch")
assoc_fmla <- ~ log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) +
                log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + prep_batch
cells_to_remove <- if (qc == "default") qc_excl$excluded_idx else integer(0)

# ---- sceptre -----------------------------------------------------------------
scep <- import_data(
  response_matrix        = response_matrix,
  grna_matrix            = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi                    = "high",
  extra_covariates       = cell_covariates[, covariates_to_use],
  response_names         = NA_character_
)

no_pairs <- data.frame(grna_target = character(0), response_id = character(0))
if (use_discovery) {
  disc_pairs <- readRDS(file.path(in_root, "shared", "discovery_pairs.rds"))
  pc_pairs   <- no_pairs
  cat("discovery pairs:", nrow(disc_pairs), "\n")
  stopifnot(nrow(disc_pairs) == N_PAIRS,
            all(disc_pairs$grna_target %in% grna_target_df$grna_target),
            all(disc_pairs$response_id %in% rownames(response_matrix)))
} else {
  disc_pairs <- no_pairs
  pc_pairs   <- construct_positive_control_pairs(scep)
  cat("positive control pairs:", nrow(pc_pairs), "\n")
  stopifnot(nrow(pc_pairs) == N_PAIRS, all(pc_pairs$grna_target == pc_pairs$response_id))
}

scep <- set_analysis_parameters(
  sceptre_object            = scep,
  discovery_pairs           = disc_pairs,
  positive_control_pairs    = pc_pairs,
  side                      = "both",
  grna_integration_strategy = "union",
  formula_object            = assoc_fmla,
  resampling_approximation  = "skew_normal",
  control_group             = "complement",
  resampling_mechanism      = "crt",
  multiple_testing_method   = "BH",
  multiple_testing_alpha    = 0.1
)

# binary input + threshold 1 = sceptre's assignment is exactly grna_matrix
scep <- assign_grnas(sceptre_object = scep, method = "thresholding", threshold = 1)

scep <- run_qc(
  sceptre_object            = scep,
  n_nonzero_trt_thresh      = 0L,           # no pairwise QC
  n_nonzero_cntrl_thresh    = 0L,
  response_n_umis_range     = c(0, 1),      # sceptre's own cell filters remove nothing
  response_n_nonzero_range  = c(0, 1),
  p_mito_threshold          = 1,
  additional_cells_to_remove = cells_to_remove
)
crm <- scep@cell_removal_metrics
cat("\ncell_removal_metrics:\n"); print(crm)
cat("cells in use:", length(scep@cells_in_use), "\n")
stopifnot(crm[["n_cells_rm_user_specified"]] == length(cells_to_remove),
          crm[["n_cells_rm_total"]] == length(cells_to_remove),
          length(scep@cells_in_use) == N_CELLS - length(cells_to_remove))
print(scep)

analysis_fn <- if (use_discovery) run_discovery_analysis else run_power_check
scep <- analysis_fn(
  sceptre_object = scep,
  output_amount  = 2,
  print_progress = TRUE,
  parallel       = FALSE,
  n_processors   = "auto",
  log_dir        = tempdir()
)
results <- get_result(scep, if (use_discovery) "run_discovery_analysis" else "run_power_check")
stopifnot(nrow(results) == N_PAIRS)

# ---- treatment-cell counts (cells in use, any guide of the target) -----------
in_use <- scep@cells_in_use
n_trt <- vapply(results$grna_target, function(g) {
  rows <- which(grna_target_df$grna_target == g)
  sum(Matrix::colSums(grna_matrix[rows, in_use, drop = FALSE]) > 0)
}, integer(1))
results$n_trt_cells <- unname(n_trt)
results$dataset <- dataset; results$method <- method; results$qc <- qc

write.csv(results, out_csv, row.names = FALSE)
cat("\nwrote", out_csv, "\n")
cat("p_value: NA", sum(is.na(results$p_value)), "; < 1e-10", sum(results$p_value < 1e-10, na.rm = TRUE),
    "; > 0.05", sum(results$p_value > 0.05, na.rm = TRUE), "\n")
