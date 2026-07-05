#!/usr/bin/env Rscript
# Calibrate a Gasperini-like simulation: reuse the paper's sum-process simulator
# (grna_simulator_iid_sum_process) with REAL Gasperini per-cell library sizes,
# and grid-search perturbation params to match Gasperini's measured separation
# (median mode-gap ~2.2 on log1p, perturbed mode ~6 counts, shallow valley ~0.46).

suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
source(file.path(GA, "grna-simulator", "sims-for-paper.R"))   # grna_simulator_iid_sum_process

# real Gasperini per-cell library sizes (high-MOI structure)
gm <- readRDS(file.path(DATA, "gasperini", "sceptre", "grna_matrix.rds"))
set.seed(1)
cell_sel <- sort(sample(ncol(gm), 50000))
n_nonzero <- Matrix::colSums(gm[, cell_sel] > 0)
lib <- n_nonzero / mean(n_nonzero)
cat("Gasperini lib sizes: mean nnz/cell =", round(mean(n_nonzero), 1),
    "(high MOI); scaled lib range", round(range(lib), 2), "\n")

measure <- function(params, num_guides = 40, seed = 10) {
  sim <- grna_simulator_iid_sum_process(num_guides, lib, params, seed)
  gmr <- as(sim$grna_matrix_sim, "RsparseMatrix")
  tru <- sim$true_perts
  do.call(rbind, lapply(seq_len(num_guides), function(g) {
    counts <- as.numeric(gmr[g, ]); v <- smoothed_valley_threshold(counts)
    tr <- as.numeric(tru[g, ]) > 0 & counts > 0
    ov <- if (sum(tr) > 0 && is.finite(v$t)) mean(counts[tr] < v$t) else NA_real_
    data.frame(bimodal = isTRUE(v$ok),
               gap   = if (isTRUE(v$ok)) log1p(v$mode2) - log1p(v$mode1) else NA_real_,
               depth = if (isTRUE(v$ok)) 1 - v$valley_h / v$mode_h else NA_real_,
               m2    = if (isTRUE(v$ok)) v$mode2 else NA_real_, overlap = ov)
  }))
}

# low per-guide background (Gasperini, like Replogle, is mostly 0/1 per guide);
# the distinguishing feature is a WEAK perturbed signal sitting just above it.
bg <- list(prob_endo = 0.01, mu_endo = 2, theta_endo = 3, mu_exo = 0.05)
grid <- expand.grid(mu_pert = c(2, 3, 4, 6, 8), theta_pert = c(1, 2, 3))
PROB_PERT <- 0.03   # higher rate so the small perturbed bump stays detectable

cat("\nTARGET (real Gasperini): gap~2.2  depth~0.46  m2~6  bimodal~97%\n\n")
cat(sprintf("%6s %6s | %5s %6s %6s %6s %7s\n",
            "muP", "thP", "%bim", "gap", "depth", "m2", "overlap"))
res <- list()
for (i in seq_len(nrow(grid))) {
  params <- c(bg, list(prob_pert = PROB_PERT,
                       mu_pert = grid$mu_pert[i], theta_pert = grid$theta_pert[i]))
  m <- measure(params)
  s <- c(mu_pert = grid$mu_pert[i], theta_pert = grid$theta_pert[i],
         pct_bimodal = round(100 * mean(m$bimodal), 0),
         gap = round(median(m$gap, na.rm = TRUE), 2),
         depth = round(median(m$depth, na.rm = TRUE), 2),
         m2 = round(median(m$m2, na.rm = TRUE)),
         overlap = round(median(m$overlap, na.rm = TRUE), 2))
  res[[i]] <- s
  cat(sprintf("%6.0f %6.0f | %5.0f %6.2f %6.2f %6.0f %7.2f\n",
              s["mu_pert"], s["theta_pert"], s["pct_bimodal"], s["gap"],
              s["depth"], s["m2"], s["overlap"]))
}
