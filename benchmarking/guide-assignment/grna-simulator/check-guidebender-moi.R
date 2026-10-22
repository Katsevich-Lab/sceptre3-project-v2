#!/usr/bin/env Rscript
# One GuideBender Gasperini dataset: its true MOI (from the ground truth) against
# the real Gasperini MOI estimated with a threshold of 5.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()
root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data")
ds <- "sim_gasperini_g694_c11000"
THR <- 5

p <- readRDS(file.path(root, ds, "sim_params.rds"))
cat("== parameters of", ds, "\n")
str(p[setdiff(names(p), "targets")], give.attr = FALSE)

truth  <- as(as(readRDS(file.path(root, ds, "true_pert_matrix.rds")), "CsparseMatrix"), "generalMatrix")
counts <- as(as(Matrix::readMM(file.path(root, ds, "cleanser", "grna_matrix.mtx")), "CsparseMatrix"), "generalMatrix")
stopifnot(identical(dim(truth), dim(counts)))
truth <- truth != 0                                      # logical, TRUE = truly perturbed
n_cells <- ncol(counts)

# Counts at the truly perturbed pairs. A true pair with count 0 has no entry in
# `counts`, so zeros are the true pairs missing from the product.
at_true <- counts * truth
n_true  <- sum(truth)
x_true  <- at_true@x[at_true@x > 0]
n_zero  <- n_true - length(x_true)

# The real matrix, thresholded exactly as the targets were.
real <- readRDS(file.path(root, "gasperini", "fishashplus", "grna_matrix.rds"))
real_thr_moi <- sum(real@x >= THR) / ncol(real)

cat(sprintf("\n== MOI (per cell)\n"))
cat(sprintf("  real Gasperini, entries >= %d      %7.2f   (%d guides x %d cells)\n",
            THR, real_thr_moi, nrow(real), ncol(real)))
cat(sprintf("  simulated, TRUE MOI (ground truth) %7.2f   (%d guides x %d cells)\n",
            n_true / n_cells, nrow(counts), n_cells))
cat(sprintf("  simulated, entries >= %d (reference) %7.2f\n", THR, sum(counts@x >= THR) / n_cells))

cat(sprintf("\n== counts at the simulation's %d true perturbations (%.2f per cell)\n",
            n_true, n_true / n_cells))
for (b in list(c(0, 0), c(1, 1), c(2, 4), c(THR, Inf))) {
  n <- if (b[1] == 0) n_zero else sum(x_true >= b[1] & x_true <= b[2])
  lab <- if (b[1] == 0) "0 (no entry)" else if (is.infinite(b[2])) sprintf(">= %d", b[1]) else
    if (b[1] == b[2]) sprintf("%d", b[1]) else sprintf("%d-%d", b[1], b[2])
  cat(sprintf("  count %-13s %9d  %5.1f%%   %6.2f per cell\n", lab, n, 100 * n / n_true, n / n_cells))
}
cat(sprintf("  mean count at true pairs %.2f; median %.0f (zeros included)\n",
            sum(x_true) / n_true, stats::median(c(rep(0, n_zero), x_true))))

# What the sub-threshold nonzeros are made of, and how they compare with real.
sub_sim <- counts@x > 0 & counts@x < THR
sim_sub_true <- sum(at_true@x > 0 & at_true@x < THR)
cat(sprintf("\n== nonzero entries below %d, per cell\n", THR))
cat(sprintf("  simulated %6.2f, of which truly perturbed %6.2f (%.1f%%)\n",
            sum(sub_sim) / n_cells, sim_sub_true / n_cells, 100 * sim_sub_true / sum(sub_sim)))
cat(sprintf("  real      %6.2f\n", sum(real@x > 0 & real@x < THR) / ncol(real)))
cat(sprintf("  share equal to 1: simulated %.3f, real %.3f\n",
            mean(counts@x[sub_sim] == 1), mean(real@x[real@x > 0 & real@x < THR] == 1)))
