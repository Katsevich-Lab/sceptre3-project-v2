#!/usr/bin/env Rscript
# How much do nonzeros per cell vary between GuideBender draws with the same
# parameters? Diagnostic.
#
#   Rscript check-guidebender-seed-variation.R
#
# For Replogle rungs G = 299 and G = 558, simulates 10 seeds of 10,000 cells
# each, with the snr and lambda that tonight's calibration chose
# (sim_scaling_calibration_replogle.csv), and reports nonzeros per cell and
# thresholded MOI per seed. Cells are independent given the guide-level draw
# (guide abundances and size factors), so 10,000 cells pin a draw's nonzeros per
# cell to well under 1%; spread across seeds beyond that is guide-level.

suppressPackageStartupMessages({ library(Matrix); library(SummarizedExperiment) })
source("~/.Rprofile")                                    # .get_config_path()
.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "grna-sim-utils.R"))

reg <- SCALING_REGIMES$replogle
tg  <- measure_real(reg)
cal <- read.csv(file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                          "guide_assignment", "input_data", "sim_scaling_calibration_replogle.csv"))
SEEDS <- 1:10; N_CELLS <- 10000L

res <- do.call(rbind, lapply(c(299L, 558L), function(G) {
  row <- cal[cal$role == "rung" & cal$n_guides == G, ]
  stopifnot(nrow(row) == 1)
  do.call(rbind, lapply(SEEDS, function(s) {
    sim <- simulate_regime(tg, G, N_CELLS, row$lambda, row$snr, s)
    st  <- cal_stats(assay(sim, "counts"), reg$threshold)
    data.frame(G = G, seed = s, lambda = row$lambda, snr = row$snr,
               nnz_per_cell = st[["nnz_per_cell"]], moi = st[["moi"]])
  }))
}))

cat(sprintf("\nreal nonzeros per cell %.3f, real MOI %.3f\n\n", tg$nnz_per_cell, tg$moi))
print(res, row.names = FALSE, digits = 4)
cat("\nper rung across seeds:\n")
for (G in unique(res$G)) {
  v <- res$nnz_per_cell[res$G == G]
  cat(sprintf("  G = %d: nonzeros per cell mean %.2f, sd %.2f (CV %.1f%%), range %.2f to %.2f; mean vs real %+.1f%%\n",
              G, mean(v), sd(v), 100 * sd(v) / mean(v), min(v), max(v),
              100 * (mean(v) / tg$nnz_per_cell - 1)))
}
