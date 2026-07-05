#!/usr/bin/env Rscript
# Try the simple nonparametric thresholding schemes on a (simulated or real)
# guide-assignment dataset and, when ground truth is available, score them.
#
# Usage:
#   Rscript run_thresholds.R <dataset_id> [n_example_grnas]
#
# Reads input_data/<dataset_id>/sceptre/grna_matrix.rds and (if present)
# input_data/<dataset_id>/true_perturbation_status.rds. Writes metrics +
# diagnostics + plots to results/<dataset_id>/.

suppressPackageStartupMessages({library(Matrix)})
HERE     <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE))))
DATA_DIR <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(HERE, "..", "guide-assignment-pipeline", "bin", "script", "lib",
                 "threshold_methods.R"))

args       <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(args) >= 1) args[1] else "replogle-rd7_simulated_medium"
n_example  <- if (length(args) >= 2) as.integer(args[2]) else 9L
outdir     <- file.path(HERE, "results", dataset_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Dataset:", dataset_id, "\n")
grna_matrix <- readRDS(file.path(DATA_DIR, dataset_id, "sceptre", "grna_matrix.rds"))
cat("  gRNA matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")

gt_fp <- file.path(DATA_DIR, dataset_id, "true_perturbation_status.rds")
has_truth <- file.exists(gt_fp)
truth <- if (has_truth) readRDS(gt_fp) else NULL
if (has_truth) cat("  ground truth:", nrow(truth), "x", ncol(truth),
                   "| total true perturbations =", sum(truth != 0), "\n")

methods <- list(
  otsu_log1p      = otsu_threshold_log1p,
  smoothed_valley = smoothed_valley_threshold
)

# ---- run both methods -------------------------------------------------------
results <- list()
for (m in names(methods)) {
  cat("Running", m, "...\n")
  t0 <- proc.time()[3]
  results[[m]] <- assign_by_threshold(grna_matrix, methods[[m]], verbose = TRUE)
  cat("  done in", round(proc.time()[3] - t0, 1), "s;",
      "total assignments =", sum(results[[m]]$assignment_matrix), "\n")
  saveRDS(results[[m]]$assignment_matrix,
          file.path(outdir, paste0("assignment_matrix_", m, ".rds")))
  write.csv(results[[m]]$diagnostics,
            file.path(outdir, paste0("diagnostics_", m, ".csv")), row.names = FALSE)
}

# ---- scoring against ground truth (per gRNA) --------------------------------
score_grna <- function(pred_row, truth_row) {
  pred  <- pred_row != 0
  true  <- truth_row != 0
  tp <- sum(pred & true); fp <- sum(pred & !true)
  fn <- sum(!pred & true); tn <- sum(!pred & !true)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
            2 * prec * sens / (prec + sens) else NA_real_
  jac  <- if (tp + fp + fn > 0) tp / (tp + fp + fn) else NA_real_
  c(n_true = sum(true), tp = tp, fp = fp, fn = fn,
    sensitivity = sens, specificity = spec, precision = prec, f1 = f1, jaccard = jac)
}

if (has_truth) {
  # align gRNA order
  common <- intersect(rownames(grna_matrix), rownames(truth))
  if (length(common) == 0) {            # fall back to positional match
    common <- seq_len(min(nrow(grna_matrix), nrow(truth)))
    rn <- function(x, idx) x[idx, , drop = FALSE]
  } else {
    rn <- function(x, idx) x[idx, , drop = FALSE]
  }
  truth_m <- as.matrix(truth[common, , drop = FALSE])

  summary_rows <- list()
  per_grna_all <- list()
  for (m in names(methods)) {
    pred_m <- as.matrix(results[[m]]$assignment_matrix[common, , drop = FALSE])
    sc <- t(vapply(seq_along(common),
                   function(i) score_grna(pred_m[i, ], truth_m[i, ]),
                   numeric(9)))
    sc <- as.data.frame(sc); sc$grna_id <- common; sc$method <- m
    per_grna_all[[m]] <- sc
    summary_rows[[m]] <- data.frame(
      method          = m,
      n_grna          = length(common),
      n_abstained     = sum(!results[[m]]$diagnostics$ok),
      mean_sens       = mean(sc$sensitivity, na.rm = TRUE),
      mean_spec       = mean(sc$specificity, na.rm = TRUE),
      mean_f1         = mean(sc$f1, na.rm = TRUE),
      mean_jaccard    = mean(sc$jaccard, na.rm = TRUE),
      median_jaccard  = median(sc$jaccard, na.rm = TRUE),
      total_tp        = sum(sc$tp), total_fp = sum(sc$fp), total_fn = sum(sc$fn)
    )
  }
  summary_df <- do.call(rbind, summary_rows)
  per_grna   <- do.call(rbind, per_grna_all)
  write.csv(summary_df, file.path(outdir, "score_summary.csv"), row.names = FALSE)
  write.csv(per_grna,   file.path(outdir, "score_per_grna.csv"), row.names = FALSE)
  cat("\n===== SCORE SUMMARY (", dataset_id, ") =====\n", sep = "")
  print(summary_df, row.names = FALSE)
}

# ---- diagnostic plots: per-gRNA log1p histograms with chosen thresholds -----
gm  <- as(grna_matrix, "RsparseMatrix")
set.seed(1)
ng  <- min(n_example, nrow(grna_matrix))
sel <- sort(sample(seq_len(nrow(grna_matrix)), ng))
png(file.path(outdir, "example_thresholds.png"),
    width = 1200, height = 1000, res = 110)
op <- par(mfrow = c(ceiling(sqrt(ng)), ceiling(ng / ceiling(sqrt(ng)))),
          mar = c(3.2, 3.2, 2, 0.6), mgp = c(1.9, 0.6, 0))
for (g in sel) {
  counts <- as.numeric(gm[g, ])
  z <- log1p(counts)
  hist(z, breaks = 40, col = "grey85", border = "grey60",
       main = paste0(substr(rownames(grna_matrix)[g], 1, 16),
                     "  (nz=", sum(counts > 0), ")"),
       xlab = "log(1 + count)", cex.main = 0.85)
  to <- otsu_threshold_log1p(counts)$t
  tv <- smoothed_valley_threshold(counts)$t
  if (is.finite(to)) abline(v = log1p(to), col = "firebrick", lwd = 2)
  if (is.finite(tv)) abline(v = log1p(tv), col = "royalblue", lwd = 2, lty = 2)
}
par(op); dev.off()
cat("\nWrote results to", outdir, "\n")
