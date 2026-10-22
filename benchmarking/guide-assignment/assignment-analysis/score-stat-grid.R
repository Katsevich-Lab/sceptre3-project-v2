#!/usr/bin/env Rscript
# Score the MOI x signal statistical-benchmarking grid (grna-simulator/simulate-stat-grid.R).
#
#   sbatch -c 2 ../grna-simulator/run-in-image.sh score-stat-grid.R
#
# Metric: the fishash authors' (fishash_analysis bin/process_assignments.R, subset
# "full"). Truth is the simulator's ground_truth, 1{integrations >= 1}, so a true
# integration with zero counts is a false negative for every method. Over the full
# guides x cells matrix, per dataset and method:
#   Precision = TP/(TP+FP), Recall = TP/(TP+FN), F1 = their harmonic mean.
# Also recall_obs2: recall over true entries with count >= 2 only. No fishash variant
# can call a count below 2 (min_count = 2), so this separates dropout (low counts
# no method sees) from separation.
#
# Reads, under <LOCAL_BENCHMARKING_DIR>/guide_assignment/:
#   input_data/stat_grid_manifest.csv               datasets and design
#   input_data/<dataset>/true_pert_matrix.rds       truth, guides x cells, with dimnames
#   input_data/<dataset>/fishashplus/grna_matrix.rds counts, with dimnames
#   outputs/<run>/assignments_<method>_<dataset>.csv, for run in RUNS:
#     fishash      cell_id, grna_id      1-based column and row INDICES
#     fishashplus  cell_id, grna_id      NAMES
#     crispat      cell_id, grna_id      NAMES
#     pertpy       cell_id, grna_id      NAMES; grna_id lists a cell's guides joined by "+",
#                                        or is "negative" when none is assigned
#     cleanser     grna_id, cell_id, posterior   1-based INDICES, one row per nonzero
#                                        entry; assigned if posterior >= CLEANSER_CUTOFF
#   outputs/<run>/iterations.tsv (optional)  method, dataset, n_iter, converged
# Missing assignment files are listed and skipped, so a run still in progress can be
# scored and scored again later.
# Writes to outputs/stat_grid/: stat_grid_scores.csv, stat_grid_summary.csv

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()

# Searched in this order; the first run holding a method's file for a dataset wins.
# stat_grid_refit100 reruns fishash and fishash+ with the refit cap raised from 10
# to 100 (it was hit before convergence at MOI 30), superseding stat_grid_fast's.
# The _ext runs hold the MOI 10 and Phi_noise = 1 datasets (refit cap 100 already).
RUNS    <- c("stat_grid_refit100", "stat_grid_fast", "stat_grid_slow",
             "stat_grid_ext_fast", "stat_grid_ext_slow")
METHODS <- c("fishash", "fishashplus", "crispat", "pertpy", "cleanser")
# CLEANSER runs with --dc on these datasets (run_cleanser.py); 0.5 is CLEANSER's
# stated cutoff for direct capture (0.8 for CROP-seq), the values that reproduce
# the fishash paper's barnyard Table 2.
CLEANSER_CUTOFF <- 0.5

base    <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
in_root <- file.path(base, "input_data")
out_dir <- file.path(base, "outputs", "stat_grid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- read.csv(file.path(in_root, "stat_grid_manifest.csv"))

find_assignments <- function(method, ds) {
  fp <- file.path(base, "outputs", RUNS, sprintf("assignments_%s_%s.csv", method, ds))
  fp[file.exists(fp)][1]
}

#' (row, column) indices of assigned entries on the truth's grid.
read_pairs <- function(fp, method, truth) {
  a <- read.csv(fp, colClasses = "character")
  if (method == "fishash") {
    i <- as.integer(a$grna_id); j <- as.integer(a$cell_id)
  } else if (method == "cleanser") {
    keep <- as.numeric(a$posterior) >= CLEANSER_CUTOFF
    i <- as.integer(a$grna_id[keep]); j <- as.integer(a$cell_id[keep])
  } else {
    if (method == "pertpy") {
      # pertpy labels a cell with no guide assigned "negative"; drop those rows.
      a <- a[a$grna_id != "negative", , drop = FALSE]
      g <- strsplit(a$grna_id, "+", fixed = TRUE)
      a <- data.frame(cell_id = rep(a$cell_id, lengths(g)), grna_id = unlist(g))
    }
    i <- match(a$grna_id, rownames(truth)); j <- match(a$cell_id, colnames(truth))
    if (anyNA(i) || anyNA(j))
      stop(method, ": names not in the dataset, e.g. ",
           paste(head(unique(c(a$grna_id[is.na(i)], a$cell_id[is.na(j)])), 3), collapse = ", "))
  }
  stopifnot(all(i >= 1 & i <= nrow(truth)), all(j >= 1 & j <= ncol(truth)),
            !anyDuplicated(cbind(i, j)))
  cbind(i, j)
}

rows <- list(); missing <- character(0)
for (ds in manifest$dataset) {
  truth  <- as(readRDS(file.path(in_root, ds, "true_pert_matrix.rds")), "CsparseMatrix")
  counts <- readRDS(file.path(in_root, ds, "fishashplus", "grna_matrix.rds"))
  stopifnot(identical(dimnames(truth), dimnames(counts)))
  obs2   <- truth & (counts >= 2)
  n_true <- sum(truth); n_obs2 <- sum(obs2)
  for (m in METHODS) {
    fp <- find_assignments(m, ds)
    if (is.na(fp)) { missing <- c(missing, sprintf("%s %s", m, ds)); next }
    ij <- read_pairs(fp, m, truth)
    tp      <- sum(truth[ij]); n_assigned <- nrow(ij)
    tp_obs2 <- sum(obs2[ij])
    P <- tp / n_assigned; R <- tp / n_true
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = ds, method = m, n_assigned = n_assigned, TP = tp, FP = n_assigned - tp,
      FN = n_true - tp, Precision = P, Recall = R, F1 = 2 * P * R / (P + R),
      recall_obs2 = tp_obs2 / n_obs2, frac_true_obs2 = n_obs2 / n_true)
  }
}
if (length(missing)) cat("missing (skipped):\n ", paste(missing, collapse = "\n  "), "\n")

scores <- merge(manifest[, c("dataset", "phi_noise", "moi", "signal_per_integration")],
                do.call(rbind, rows), by = "dataset")
it <- lapply(file.path(base, "outputs", RUNS, "iterations.tsv"), function(f)
  if (file.exists(f)) read.delim(f) else NULL)
it <- do.call(rbind, it)
if (!is.null(it)) it <- it[!duplicated(it[, c("method", "dataset")]), ]   # same run order as above
if (!is.null(it)) scores <- merge(scores, it, by = c("method", "dataset"), all.x = TRUE)
scores <- scores[order(scores$phi_noise, scores$moi, -scores$signal_per_integration,
                       scores$method, scores$dataset), ]
write.csv(scores, file.path(out_dir, "stat_grid_scores.csv"), row.names = FALSE)

key  <- scores[, c("phi_noise", "moi", "signal_per_integration", "method")]
vals <- c("Precision", "Recall", "F1", "recall_obs2")
summ <- aggregate(scores[, vals], by = key, FUN = mean)
summ$F1_sd <- aggregate(scores$F1, by = key, FUN = sd)$x
summ$n     <- aggregate(scores$F1, by = key, FUN = length)$x
summ <- summ[order(summ$phi_noise, summ$moi, -summ$signal_per_integration, summ$method), ]
write.csv(summ, file.path(out_dir, "stat_grid_summary.csv"), row.names = FALSE)

show <- summ
num  <- vapply(show, is.numeric, logical(1))
show[num] <- lapply(show[num], signif, 3)
print(show, row.names = FALSE)
cat("\nscores  ->", file.path(out_dir, "stat_grid_scores.csv"),
    "\nsummary ->", file.path(out_dir, "stat_grid_summary.csv"), "\n")
