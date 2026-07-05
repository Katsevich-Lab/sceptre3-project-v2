#!/usr/bin/env Rscript
# Generate a Gasperini-calibrated simulation that fills the low-separation regime
# the Replogle-calibrated sims miss. Uses the paper's sum-process simulator with
# REAL Gasperini per-cell library sizes; 3 perturbation levels spanning a
# difficulty gradient anchored at real Gasperini's separation (gap ~2.2).
# Saves grna_matrix.rds + true_pert_matrix.rds in the framework's layout, then
# verifies separation (sep_row) and scores Otsu/valley against ground truth.

suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
source(file.path(GA, "grna-simulator", "sims-for-paper.R"))

DS_NAME <- "sims_gasperini_calibrated"
out_dir <- file.path(DATA, DS_NAME); dir.create(file.path(out_dir, "sceptre"), recursive = TRUE, showWarnings = FALSE)

# real Gasperini library sizes
gm_real <- readRDS(file.path(DATA, "gasperini", "sceptre", "grna_matrix.rds"))
set.seed(1)
cell_sel  <- sort(sample(ncol(gm_real), 50000))
n_nonzero <- Matrix::colSums(gm_real[, cell_sel] > 0)
lib <- n_nonzero / mean(n_nonzero)
cell_names <- paste0("cell_idx_", cell_sel)

bg <- list(prob_endo = 0.01, mu_endo = 2, theta_endo = 3, mu_exo = 0.05)
params_list <- list(
  Gasp_P1 = c(bg, list(prob_pert = 0.04, mu_pert = 3,  theta_pert = 2)),  # hardest
  Gasp_P2 = c(bg, list(prob_pert = 0.04, mu_pert = 8,  theta_pert = 2)),  # ~ real Gasperini
  Gasp_P3 = c(bg, list(prob_pert = 0.04, mu_pert = 15, theta_pert = 3))   # a bit easier
)
NG <- 100

pieces_g <- list(); pieces_t <- list()
for (pn in names(params_list)) {
  s <- grna_simulator_iid_sum_process(NG, lib, params_list[[pn]], seed = 100 + which(names(params_list) == pn))
  gn <- paste0("sim_Gasp_", sub("Gasp_", "", pn), "_", seq_len(NG))   # sim_Gasp_P1_1 ...
  pieces_g[[pn]] <- `dimnames<-`(s$grna_matrix_sim, list(gn, cell_names))
  pieces_t[[pn]] <- `dimnames<-`(s$true_perts,      list(gn, cell_names))
}
grna_matrix     <- as(do.call(rbind, pieces_g), "CsparseMatrix")
true_pert_matrix <- as(do.call(rbind, pieces_t), "CsparseMatrix")

saveRDS(grna_matrix,      file.path(out_dir, "sceptre", "grna_matrix.rds"))
saveRDS(true_pert_matrix, file.path(out_dir, "true_pert_matrix.rds"))
saveRDS(list(params_list = params_list, num_guides_per_param = NG, num_cells = 50000),
        file.path(out_dir, "sim_params.rds"))
write.csv(data.frame(grna_id = rownames(grna_matrix),
                     grna_target = c("non-targeting", paste0("t_", seq_len(nrow(grna_matrix) - 1)))),
          file.path(out_dir, "sceptre", "grna_target_data_frame.csv"), row.names = FALSE)
cat("Saved", DS_NAME, ":", nrow(grna_matrix), "guides x", ncol(grna_matrix), "cells\n\n")

# ---- verify separation + score Otsu/valley vs ground truth ------------------
gmr <- as(grna_matrix, "RsparseMatrix")
extract_p <- function(x) vapply(strsplit(x, "_", fixed = TRUE), `[[`, character(1), 3)
score <- function(pred, truth) { tp <- sum(pred & truth)
  prec <- if (sum(pred)) tp / sum(pred) else NA; rec <- if (sum(truth)) tp / sum(truth) else NA
  jac <- if (sum(pred | truth)) tp / sum(pred | truth) else NA; c(prec, rec, jac) }
rows <- list()
for (g in seq_len(nrow(gmr))) {
  counts <- as.numeric(gmr[g, ]); tr <- as.numeric(true_pert_matrix[g, ]) > 0
  v <- smoothed_valley_threshold(counts); o <- otsu_threshold_log1p(counts)
  sv <- if (is.finite(v$t)) counts >= v$t else rep(FALSE, length(counts))
  ot <- if (is.finite(o$t)) counts >= o$t else rep(FALSE, length(counts))
  rows[[g]] <- data.frame(P = extract_p(rownames(gmr)[g]),
    bimodal = isTRUE(v$ok),
    gap = if (isTRUE(v$ok)) log1p(v$mode2) - log1p(v$mode1) else NA_real_,
    otsu_jac = score(ot, tr)[3], valley_jac = score(sv, tr)[3],
    otsu_prec = score(ot, tr)[1], valley_prec = score(sv, tr)[1])
}
df <- do.call(rbind, rows)
agg <- aggregate(cbind(gap, otsu_jac, valley_jac, otsu_prec, valley_prec) ~ P, df,
                 function(x) round(mean(x, na.rm = TRUE), 3))
bim <- aggregate(bimodal ~ P, df, function(x) round(100 * mean(x)))
cat("Per-level separation + Otsu/valley accuracy (vs ground truth):\n")
print(merge(agg, bim, by = "P"), row.names = FALSE)
cat("\nReal Gasperini reference: gap~2.2, ~98% bimodal.\n")

# ---- visual: simulated guide vs a real Gasperini guide ----------------------
real_umis <- as.numeric(gm_real[order(Matrix::rowSums(gm_real > 0), decreasing = TRUE)[500], ])
gi <- which(extract_p(rownames(gmr)) == "P2")[1]
p <- plot_umi_histogram_real_vs_sim(
  umis_real = real_umis, umis_sim = as.numeric(gmr[gi, ]),
  is_pert = as.numeric(true_pert_matrix[gi, ]) > 0,
  title = "Real Gasperini guide  vs  Gasperini-calibrated sim guide (P2)")
ggsave(file.path(HERE, "results", "gasperini_sim_validation.png"), p, width = 9, height = 4, dpi = 120)
cat("\nWrote results/gasperini_sim_validation.png\n")
