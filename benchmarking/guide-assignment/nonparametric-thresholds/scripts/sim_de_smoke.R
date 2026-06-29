#!/usr/bin/env Rscript
# =============================================================================
# DE pushthrough smoke test: run all 8 in-process methods through the full
# pos-control (power) + neg-control (calibration) pipeline on a subsample of
# Gasperini (300 guides x 30k cells).  Confirms the cross-method comparison
# works end-to-end and produces realistic-looking results before scaling up.
# Outputs: results/sim_framework/de/smoke_gasperini/{<method>_{pos,neg}/}
#          results/sim_framework/de/smoke_gasperini.csv  (one row per method)
# =============================================================================
suppressPackageStartupMessages({library(Matrix); library(sceptre)})
source(file.path(getwd(), "scripts", "sim_de_lib.R"))

real <- load_real_dataset("gasperini", max_cells = 30000, max_guides = 300, seed = 1)
methods <- c("thresh3", "thresh10", "ambient", "otsu", "valley", "fishash", "depth_fix", "sceptre")
A <- run_methods_on_real(real, methods = methods, verbose = TRUE)
out_root <- file.path(DE_OUT, "smoke_gasperini")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

rows <- list()
for (m in names(A)) {
  d_pos <- file.path(out_root, paste0(m, "_pos"))
  d_neg <- file.path(out_root, paste0(m, "_neg"))
  unlink(c(d_pos, d_neg), recursive = TRUE)
  cat(sprintf("\n=== %s ===\n", m))
  cat("[build pos]\n"); pinfo <- tryCatch(build_pos_control_input(real, A[[m]], d_pos,
                            n_targets = 20, max_cells = 8000),
                            error = function(e) { cat("  ERR ", conditionMessage(e), "\n"); NULL })
  cat("[build neg]\n"); ninfo <- tryCatch(build_neg_control_input(real, A[[m]], d_neg,
                            n_genes = 30, max_cells = 6000),
                            error = function(e) { cat("  ERR ", conditionMessage(e), "\n"); NULL })
  ps <- if (!is.null(pinfo)) tryCatch({csv <- run_de_sceptre(d_pos, "pos", "gasperini"); score_de(csv, "pos") },
                                       error = function(e) { cat("  pos ERR ", conditionMessage(e), "\n"); NULL })
  ns <- if (!is.null(ninfo)) tryCatch({csv <- run_de_sceptre(d_neg, "neg", "gasperini"); score_de(csv, "neg") },
                                       error = function(e) { cat("  neg ERR ", conditionMessage(e), "\n"); NULL })
  rows[[m]] <- data.frame(method = m,
    n_calls = sum(A[[m]]),
    pos_pairs = ifelse(is.null(ps), NA, ps$n_pairs),
    power_q10 = ifelse(is.null(ps), NA, ps$power_q),
    pos_median_p = ifelse(is.null(ps), NA, ps$median_p),
    neg_pairs = ifelse(is.null(ns), NA, ns$n_pairs),
    realized_t1e_005 = ifelse(is.null(ns), NA, ns$realized_t1e_005),
    realized_t1e_001 = ifelse(is.null(ns), NA, ns$realized_t1e_001),
    ks_uniform = ifelse(is.null(ns), NA, as.numeric(ns$ks_uniform)))
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write.csv(res, file.path(out_root, "smoke_gasperini.csv"), row.names = FALSE)
cat("\n=== DE smoke summary (Gasperini 300 guides x 30k cells) ===\n")
print(format(res, digits = 3), row.names = FALSE)
cat("\nwrote", file.path(out_root, "smoke_gasperini.csv"), "\n")
