#!/usr/bin/env Rscript
# Score the 2x2 statistical-benchmarking run (MOI x signal per integration).
#
#   sbatch ../grna-simulator/run-in-image.sh score-stat-2x2.R [run_id]   # default run_id stat2x2
#
# Metric: the fishash authors' (fishash_analysis bin/process_assignments.R,
# subset "full"). Truth is the simulator's ground_truth, 1{integrations >= 1},
# so a true integration with zero counts is a false negative for every method.
# Over the full guides x cells matrix, per dataset and method:
#   TP = assigned and true, FP = assigned and not true, FN = true and not assigned,
#   Precision = TP/(TP+FP), Recall = TP/(TP+FN), F1 = their harmonic mean.
#
# Reads, under <LOCAL_BENCHMARKING_DIR>/guide_assignment/:
#   input_data/stat2x2_manifest.csv                 datasets and design (simulate-stat-2x2.R)
#   input_data/<dataset>/true_pert_matrix.rds       truth, guides x cells, with dimnames
#   outputs/<run_id>/assignments_<method>_<dataset>.csv
#       fishash:     cell_id, grna_id as 1-based column and row INDICES
#       fishashplus: cell_id, grna_id as NAMES (CELL_<j>, grna_<i>)
#   outputs/<run_id>/iterations.tsv (optional)      method, dataset, n_iter, converged,
#       collected from the task logs after the run (the work directory is not
#       mounted in the container)
# Writes to outputs/<run_id>/:
#   stat2x2_scores.csv    one row per dataset and method
#   stat2x2_summary.csv   per design cell and method: means over replicates, and
#                         dF1 = F1(fishashplus) - F1(fishash) per replicate

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()

args   <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args) >= 1) args[1] else "stat2x2"
METHODS <- c("fishash", "fishashplus")

base     <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
in_root  <- file.path(base, "input_data")
out_dir  <- file.path(base, "outputs", run_id)
manifest <- read.csv(file.path(in_root, "stat2x2_manifest.csv"))

#' Assigned entries as a logical guides x cells matrix on the truth's grid.
read_assigned <- function(fp, method, truth) {
  a <- read.csv(fp, colClasses = "character")
  if (method == "fishash") {
    i <- as.integer(a$grna_id); j <- as.integer(a$cell_id)
  } else {
    i <- match(a$grna_id, rownames(truth)); j <- match(a$cell_id, colnames(truth))
  }
  stopifnot(!anyNA(i), !anyNA(j), all(i <= nrow(truth)), all(j <= ncol(truth)),
            !anyDuplicated(cbind(i, j)))
  sparseMatrix(i = i, j = j, x = TRUE, dims = dim(truth), dimnames = dimnames(truth))
}

rows <- list()
for (ds in manifest$dataset) {
  truth <- as(readRDS(file.path(in_root, ds, "true_pert_matrix.rds")), "CsparseMatrix")
  n_true <- sum(truth)
  for (m in METHODS) {
    fp <- file.path(out_dir, sprintf("assignments_%s_%s.csv", m, ds))
    if (!file.exists(fp)) { cat("missing:", fp, "\n"); next }
    A  <- read_assigned(fp, m, truth)
    tp <- sum(A & truth); n_assigned <- sum(A)
    P  <- tp / n_assigned; R <- tp / n_true
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = ds, method = m, TP = tp, FP = n_assigned - tp, FN = n_true - tp,
      Precision = P, Recall = R, F1 = 2 * P * R / (P + R))
  }
}
scores <- merge(manifest[, c("dataset", "moi", "signal_per_integration")],
                do.call(rbind, rows), by = "dataset")

it_fp <- file.path(out_dir, "iterations.tsv")
if (file.exists(it_fp)) {
  scores <- merge(scores, read.delim(it_fp), by = c("method", "dataset"), all.x = TRUE)
} else {
  cat("no iterations.tsv; iteration counts left out\n")
}
scores <- scores[order(scores$moi, scores$signal_per_integration, scores$dataset, scores$method), ]
write.csv(scores, file.path(out_dir, "stat2x2_scores.csv"), row.names = FALSE)

# ---- summary per design cell -------------------------------------------------
cells <- unique(scores[, c("moi", "signal_per_integration")])
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
  d <- merge(scores, cells[k, ])
  w <- merge(d[d$method == "fishash", c("dataset", "F1")],
             d[d$method == "fishashplus", c("dataset", "F1")],
             by = "dataset", suffixes = c("_fishash", "_fishashplus"))
  dF1 <- w$F1_fishashplus - w$F1_fishash
  per <- lapply(METHODS, function(m) {
    x <- d[d$method == m, ]
    out <- data.frame(n = nrow(x), Precision = mean(x$Precision), Recall = mean(x$Recall),
                      F1 = mean(x$F1))
    if ("n_iter" %in% names(x)) out$max_iter <- max(x$n_iter)
    names(out) <- paste(m, names(out), sep = "_")
    out
  })
  cbind(cells[k, ], do.call(cbind, per),
        dF1_mean = mean(dF1), dF1_sd = sd(dF1), dF1_min = min(dF1), dF1_max = max(dF1),
        dF1_n_positive = sum(dF1 > 0), dF1_n = length(dF1))
}))
write.csv(summ, file.path(out_dir, "stat2x2_summary.csv"), row.names = FALSE)

show <- summ
num  <- vapply(show, is.numeric, logical(1))
show[num] <- lapply(show[num], signif, 4)
print(t(show), quote = FALSE)
cat("\nscores  ->", file.path(out_dir, "stat2x2_scores.csv"),
    "\nsummary ->", file.path(out_dir, "stat2x2_summary.csv"), "\n")
