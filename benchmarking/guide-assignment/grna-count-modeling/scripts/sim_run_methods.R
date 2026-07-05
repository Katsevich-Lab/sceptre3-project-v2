#!/usr/bin/env Rscript
# =============================================================================
# Run the in-process (R) method panel on every Model A/B dataset, score vs the
# injected ground truth, and dump the count matrix for the out-of-process
# methods (geomux, crispat).  Methods:
#   fixed thresholds (thresh3/10/20), ambient test, Otsu, smoothed valley,
#   fishash, depth_fix (ours), sceptre (current/mixture).
# Truth = "integrated AND observed" (the recoverable target).
# Outputs: results/sim_framework/scores_R.csv,  manifest_all.csv,
#          and per-dataset ext_counts.mtx / ext_guides.txt / ext_cells.txt.
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(sceptre) })
source(file.path(getwd(), "scripts", "sim_lib.R"))
OUT <- SIMFW()
R_METHODS <- c("thresh3","thresh10","thresh20","ambient","otsu","valley","fishash","depth_fix","sceptre")

# unified manifest (common columns across A and B)
manA <- read.csv(file.path(OUT, "manifest_modelA.csv"))
manB <- read.csv(file.path(OUT, "manifest_modelB.csv"))
unify <- function(m, regime_col) data.frame(model = m$model, id = m$id, chemistry = m$chemistry,
  regime = m[[regime_col]], purity = m$purity, mu_level = m$mu_level, mu_pert = m$mu_pert,
  moi_level = m$moi_level, moi = m$moi)
man <- rbind(unify(manA, "sample"), unify(manB, "regime"))
write.csv(man, file.path(OUT, "manifest_all.csv"), row.names = FALSE)
cat("unified manifest:", nrow(man), "datasets\n")

rows <- list()
for (k in seq_len(nrow(man))) {
  id <- man$id[k]; d <- load_dataset(id)
  counts <- as(d$counts, "CsparseMatrix")
  if (is.null(rownames(counts))) rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("c", seq_len(ncol(counts)))
  tr <- truth_observed(d$Z, counts); dimnames(tr) <- dimnames(counts)
  cat(sprintf("[%2d/%2d] %s  (G=%d C=%d)\n", k, nrow(man), id, nrow(counts), ncol(counts)))
  # sceptre's "low"/"high" MOI knob: lab convention is moi==1 -> low, else high.
  sc_moi <- if (is.finite(man$moi[k]) && man$moi[k] >= 5) "high" else "low"
  A <- run_methods(counts, methods = R_METHODS, sceptre_moi = sc_moi, verbose = FALSE)
  sc <- score_panel(A, tr)
  sc <- cbind(man[k, ], sc[, c("method","jaccard","precision","recall","fdr_pooled","recall_pooled","n_pred","n_true")])
  rows[[id]] <- sc
  # dump inputs for out-of-process methods; skip if fresh wrt counts.rds
  ext_mtx <- file.path(d$dir, "ext_counts.mtx")
  if (!file.exists(ext_mtx) ||
      file.mtime(ext_mtx) < file.mtime(file.path(d$dir, "counts.rds"))) {
    writeMM(counts, ext_mtx)
    writeLines(rownames(counts), file.path(d$dir, "ext_guides.txt"))
    writeLines(colnames(counts), file.path(d$dir, "ext_cells.txt"))
  }
}
scores <- do.call(rbind, rows); rownames(scores) <- NULL
write.csv(scores, file.path(OUT, "scores_R.csv"), row.names = FALSE)
cat("\nwrote scores_R.csv (", nrow(scores), "rows ) + per-dataset ext_counts.mtx\n")
