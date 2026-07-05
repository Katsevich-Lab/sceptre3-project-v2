#!/usr/bin/env Rscript
# ============================================================================
# Reproduce fishash's reported Simpson's-paradox phenomenon (paper Figure 7).
#
# Claim under test: a one-sided Fisher test on the full (signal+noise) margins
# produces FALSE POSITIVES via Simpson's paradox, and the failure is governed
# by the correlation between the SIGNAL guide composition and the NOISE (ambient)
# guide composition (Proposition 1: paradox vanishes when they are equal).
#   - LOW correlation  (noise composition unrelated to signal): paradox strongest,
#     uncorrected Fisher fails FDR control (paper: 20/20 replicates).
#   - HIGH correlation (noise composition == signal):           paradox absent.
#
# We use the authors' OWN simulator (simulate_guidebender2) and the REAL package:
#   - uncorrected = fishash(refit = 0)  -> plain Fisher on full margins
#                   (this is also exactly our ambient/geomux-core test in spirit)
#   - corrected   = fishash(refit = 10) -> iterative Simpson's-paradox correction
# The knob `endo_shape_flat` sets the regime:
#   1   -> uniform Dirichlet noise comp  (LOW corr w/ signal)
#   0.5 -> Dirichlet centred on signal   (MEDIUM corr)
#   0   -> Dirichlet proportional to signal (HIGH corr)
#
# Figure-7 settings: 20 guides, median 100 UMIs/cell, SNR 4, MOI 0.3,
# hurdle 0.1, exogenous noise = 0 (frac_noise_endo = 1).
# ============================================================================

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})

HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)                       # scripts/ -> folder root
RES  <- file.path(HERE, "results"); dir.create(RES, showWarnings = FALSE)

regimes <- c(low = 1, medium = 0.5, high = 0)   # endo_shape_flat
N_REP   <- 10
N_CELLS <- 4000
FDR_Q   <- 0.05

# realized error metrics of a logical assignment matrix vs logical ground truth
score <- function(assigned, truth) {
  assigned <- as(assigned, "lgCMatrix"); truth <- as(truth, "lgCMatrix")
  n_assigned <- sum(assigned); n_truth <- sum(truth)
  tp <- sum(assigned & truth); fp <- n_assigned - tp
  c(fdr       = if (n_assigned > 0) fp / n_assigned else 0,
    precision = if (n_assigned > 0) tp / n_assigned else NA,
    recall    = if (n_truth   > 0) tp / n_truth   else NA,
    n_assigned = n_assigned)
}

rows <- list()
for (rg in names(regimes)) {
  for (rep in seq_len(N_REP)) {
    set.seed(1000 * which(names(regimes) == rg) + rep)
    sim <- simulate_guidebender2(
      n_guides = 20, n_cells = N_CELLS,
      moi = 0.3, hurdle_prob = 0.1,
      snr = 4, count_per_cell = 100,
      frac_noise_endo = 1,                 # exogenous noise = 0
      endo_shape_flat = regimes[[rg]],
      return_sparse_only = TRUE
    )
    counts <- assay(sim, "counts")
    truth  <- assay(sim, "ground_truth")

    unc <- score(assay(fishash(counts, refit = 0,  padj_cutoff = FDR_Q), "assigned"), truth)
    cor <- score(assay(fishash(counts, refit = 10, padj_cutoff = FDR_Q), "assigned"), truth)

    rows[[length(rows) + 1]] <- data.frame(
      regime = rg, rep = rep,
      unc_fdr = unc["fdr"], unc_prec = unc["precision"], unc_rec = unc["recall"],
      cor_fdr = cor["fdr"], cor_prec = cor["precision"], cor_rec = cor["recall"])
    cat(sprintf("  %-6s rep %2d | uncorrected FDR %.3f rec %.3f | corrected FDR %.3f rec %.3f\n",
                rg, rep, unc["fdr"], unc["recall"], cor["fdr"], cor["recall"]))
  }
}
df <- do.call(rbind, rows); rownames(df) <- NULL
write.csv(df, file.path(RES, "simpson_paradox_repro.csv"), row.names = FALSE)

cat("\n================ SUMMARY (mean over", N_REP, "replicates) ================\n")
cat(sprintf("%-8s | %-28s | %-28s\n", "regime",
            "UNCORRECTED (refit=0)", "CORRECTED (refit=10)"))
cat(sprintf("%-8s | %-28s | %-28s\n", "",
            "FDR    prec   rec  #fail", "FDR    prec   rec  #fail"))
for (rg in names(regimes)) {
  s <- df[df$regime == rg, ]
  cat(sprintf("%-8s | %.3f  %.3f  %.3f  %2d/%d | %.3f  %.3f  %.3f  %2d/%d\n", rg,
      mean(s$unc_fdr), mean(s$unc_prec), mean(s$unc_rec), sum(s$unc_fdr > FDR_Q), N_REP,
      mean(s$cor_fdr), mean(s$cor_prec), mean(s$cor_rec), sum(s$cor_fdr > FDR_Q), N_REP))
}
cat("\n#fail = replicates whose realized FDR exceeds the nominal", FDR_Q, "\n")
cat("wrote", file.path(RES, "simpson_paradox_repro.csv"), "\n")
