#!/usr/bin/env Rscript
# =============================================================================
# Sensitivity sweep: rebuild Model B regimes at SIM_D_SIGMA=0.3 (vs the primary
# 0.05), re-run the in-process R panel + geomux + crispat on just those regimes,
# and report the change in the depth_fix - fishash gap as the structural
# fairness check the review flagged.
#
# The primary regimes (d_sigma=0.05) live under datasets/B__regime__<ds>/.
# This driver:
#   1. Snapshots the primary regimes.csv + manifest_modelB.csv + the Model B
#      rows of scores_all.csv to *_dsig005.* so we can compare later.
#   2. Renames datasets/B__regime__<ds>/  ->  datasets/B__regime_dsig005__<ds>/
#      (preserves primary results untouched).
#   3. Re-runs sim_regimes.R with SIM_D_SIGMA=0.3 to produce a fresh set of
#      datasets/B__regime__<ds>/ at the wider DGP, plus a new regimes.csv and
#      manifest_modelB.csv.  These are SHIPPED only for the duration of this run.
#   4. (Then the caller re-runs the panel/external/scoring on the new datasets
#      and produces the comparison table.)
#
# Run from the folder root after the primary pipeline has fully completed.
# =============================================================================
suppressPackageStartupMessages(library(Matrix))
source(file.path(getwd(), "scripts", "sim_lib.R"))
OUT <- SIMFW()

# ---- 1. snapshot primary outputs -------------------------------------------
files_to_snap <- c("regimes.csv", "manifest_modelB.csv")
for (f in files_to_snap) {
  src <- file.path(OUT, f); dst <- file.path(OUT, sub("\\.csv$", "_dsig005.csv", f))
  if (file.exists(src)) { file.copy(src, dst, overwrite = TRUE); cat("snapshot:", basename(dst), "\n") }
}
# snapshot the Model B rows of scores_all.csv
if (file.exists(file.path(OUT, "scores_all.csv"))) {
  s <- read.csv(file.path(OUT, "scores_all.csv"))
  sb <- s[s$model == "B", ]
  write.csv(sb, file.path(OUT, "scores_B_dsig005.csv"), row.names = FALSE)
  cat("snapshot: scores_B_dsig005.csv (", nrow(sb), "rows )\n")
}

# ---- 2. rename primary regime dataset dirs (preserve primary on disk) -------
prim_dirs <- list.dirs(DSDIR(), recursive = FALSE)
prim_dirs <- prim_dirs[grepl("/B__regime__", prim_dirs)]
for (d in prim_dirs) {
  new <- sub("/B__regime__", "/B__regime_dsig005__", d)
  if (!dir.exists(new)) file.rename(d, new)
}
cat(sprintf("renamed %d primary regime dirs -> *_dsig005__*\n", length(prim_dirs)))
cat("ready to run: SIM_D_SIGMA=0.3 Rscript scripts/sim_regimes.R\n")
