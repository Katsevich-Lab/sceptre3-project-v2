#!/usr/bin/env Rscript
# Is thresholded MOI monotone in lambda? Diagnostic for the calibration failure
# in tests/test-scaling-utils.R [9] (G = 60, starting snr = 5).
#
#   srun -p genoa-std-mem -c 2 -t 0:30:00 ./run-in-image.sh check-moi-vs-lambda.R
#
# Same settings as that test. Prints thresholded MOI and nonzeros per cell over a
# grid of lambda, at several snr and two numbers of guides. Writes nothing.

suppressPackageStartupMessages({ library(Matrix); library(SummarizedExperiment) })
.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "grna-sim-utils.R"))

tg <- list(count_per_cell = 200, hurdle_prob = 0.1, alpha = 2)
THRESHOLD <- 5L; N_CELLS <- 4000L; SEED <- 11L
TARGET_MOI <- 3.303                                       # test [9]'s target

grid <- expand.grid(lambda = c(0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64),
                    snr = c(5, 10, 20, 40), G = c(60L, 100L))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g  <- grid[i, ]
  st <- cal_stats(assay(simulate_regime(tg, g$G, N_CELLS, g$lambda, g$snr, SEED), "counts"),
                  THRESHOLD)
  data.frame(g, moi = st[["moi"]], nnz_per_cell = st[["nnz_per_cell"]])
}))

for (G in unique(res$G)) {
  cat(sprintf("\n=== G = %d: thresholded MOI (target %.3f); rows lambda, columns snr ===\n",
              G, TARGET_MOI))
  sub <- res[res$G == G, ]
  print(round(xtabs(moi ~ lambda + snr, sub), 3))
  cat("\nnonzeros per cell:\n")
  print(round(xtabs(nnz_per_cell ~ lambda + snr, sub), 2))
}
