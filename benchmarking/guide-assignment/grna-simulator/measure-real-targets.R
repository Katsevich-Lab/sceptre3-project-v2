#!/usr/bin/env Rscript
# Summarize a real gRNA count matrix. Run once per dataset; the printed numbers
# get pasted into the REGIMES block of simulate-scaling-datasets.R.
#
#   Rscript measure-real-targets.R gasperini
#   Rscript measure-real-targets.R replogle-rd7
#
# Needs the R_461 renv library:
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# Only n_guides, moi and count_per_cell are consumed by the simulation right now.
# zero_frac and pert_rate are printed as a record: they are the targets a later
# calibration of snr / frac_noise_endo would aim at, and are worth having on file
# before those parameters get tuned.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                  # .get_config_path()

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else
  stop("usage: measure-real-targets.R <gasperini|replogle-rd7>")

# Cells above this count are treated as genuinely perturbed. Higher for replogle
# because its guides carry more UMIs per cell.
THRESHOLD <- if (grepl("replogle", DATASET)) 10L else 5L
HURDLE    <- 0.1

# simulate_guidebender2 takes the PRE-selection Poisson rate. After zero-truncation
# and the hurdle the realized mean is lambda/(1 - exp(-lambda)) * (1 - hurdle_prob),
# so invert that to get the lambda which reproduces an observed MOI.
lambda_from_moi <- function(m, hurdle = HURDLE) {
  uniroot(function(l) l / (1 - exp(-l)) * (1 - hurdle) - m,
          interval = c(1e-8, 500), tol = 1e-10)$root
}

mtx_fp <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", DATASET,
                    "cleanser", "grna_matrix.mtx")
cat("\n=== real-matrix summary:", DATASET, "===\n")
cat("reading", mtx_fp, "\n")
counts <- as(Matrix::readMM(mtx_fp), "CsparseMatrix")

n_guides <- nrow(counts)
n_cells  <- ncol(counts)
above    <- counts >= THRESHOLD

moi            <- mean(Matrix::colSums(above))
zero_frac      <- 1 - Matrix::nnzero(counts) / (as.numeric(n_guides) * n_cells)
pert_rate      <- mean(Matrix::rowSums(above)) / n_cells
count_per_cell <- median(Matrix::colSums(counts))
lambda         <- lambda_from_moi(moi)

cat(sprintf("\n  %d guides x %d cells, threshold %d\n", n_guides, n_cells, THRESHOLD))
cat(sprintf("  moi            %.3f  (mean guides/cell with count >= %d)\n", moi, THRESHOLD))
cat(sprintf("  count_per_cell %.0f  (median guide UMIs per cell)\n", count_per_cell))
cat(sprintf("  zero_frac      %.5f\n", zero_frac))
cat(sprintf("  pert_rate      %.3e  (%.0f cells/guide at this size)\n",
            pert_rate, pert_rate * n_cells))

cat("\n=== paste into simulate-scaling-datasets.R ===\n")
cat(sprintf("    n_guides       = %d,\n", n_guides))
cat(sprintf("    moi            = %.4f,   # lambda; reproduces observed moi %.2f\n",
            lambda, moi))
cat(sprintf("    count_per_cell = %.0f,\n", count_per_cell))
cat(sprintf("    threshold      = %d,\n", THRESHOLD))
cat(sprintf("\n# for the record, not yet targeted: zero_frac %.5f, pert_rate %.3e\n",
            zero_frac, pert_rate))
