#!/usr/bin/env Rscript
# Summarize a real gRNA count matrix, and test the assumptions the scaling
# simulations rest on. A diagnostic report: simulate-scaling-datasets.R measures
# the real matrix itself and reads nothing written here.
#
#   Rscript measure-real-targets.R gasperini
#   Rscript measure-real-targets.R replogle-rd7
#
# Needs the R_461 renv library:
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# The statistics come from guide_matrix_stats() in grna-sim-utils.R, which
# computes the same quantities on any counts matrix. Running it on a simulated
# matrix produces a row directly comparable to the one written here, so
# calibration is a comparison of measured numbers rather than of a derivation.
# The row is written to <dataset>_real_stats.csv for that purpose.
#
# Sections 3-5 exist to test the assumptions the parameter choices rest on. Each
# prints a measured value next to the value some assumption implies; where they
# disagree, the assumption is wrong for this dataset.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                  # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
source(file.path(here, "grna-sim-utils.R"))            # guide_matrix_stats(), gini(), cv()

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else
  stop("usage: measure-real-targets.R <gasperini|replogle-rd7>")

# Cells at or above this count are treated as genuinely perturbed. Higher for
# replogle because its guides carry more UMIs per cell.
THRESHOLD <- if (grepl("replogle", DATASET)) 10L else 5L

mtx_fp <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", DATASET,
                    "cleanser", "grna_matrix.mtx")
cat("\n=== real-matrix summary:", DATASET, "===\n")
cat("reading", mtx_fp, "\n")
counts <- as(Matrix::readMM(mtx_fp), "CsparseMatrix")

s <- guide_matrix_stats(counts, threshold = THRESHOLD, label = paste0(DATASET, "_real"))
G <- s$n_guides


# ---- 1. anchors the simulator consumes -------------------------------------
cat(sprintf("\n[1] anchors  (%d guides x %d cells, threshold %d)\n",
            s$n_guides, s$n_cells, THRESHOLD))
cat(sprintf("  moi                %8.3f   mean guides/cell at or above threshold\n", s$moi))
cat(sprintf("  count_per_cell     %8.0f   median guide UMIs per cell\n", s$umis_per_cell_median))
cat(sprintf("  n_cells / n_guides %8.3f\n", s$n_cells / s$n_guides))
cat(sprintf("  pert_per_guide     %8.0f   cells per guide at this size\n", s$pert_per_guide))
cat(sprintf("  zero_frac          %8.5f\n", s$zero_frac))
cat(sprintf("  nnz_per_cell       %8.2f   = moi + ambient nonzeros\n", s$nnz_per_cell))
cat(sprintf("  frac cells with no guide  %.4f\n", s$frac_cells_unpert))
# Identity, not an estimate: both sides count the same perturbed entries.
cat(sprintf("  check: moi == n_guides * pert_rate   %.6f vs %.6f\n",
            s$moi, G * s$pert_rate))


# ---- 2. how much of this depends on the threshold --------------------------
# Everything below splits entries by THRESHOLD, so the first question is whether
# that split is stable. A moi that slides steadily with the threshold means the
# perturbed/ambient boundary is not sharp in this dataset and every number that
# follows carries that ambiguity.
cat("\n[2] sensitivity of moi and snr to the threshold\n")
ct <- as(counts, "TsparseMatrix"); xv <- as.numeric(ct@x)
tot_per_cell <- sum(xv) / s$n_cells
cat("       T     moi   ambient UMIs/cell     snr\n")
for (T_ in c(2L, 3L, 5L, 10L, 20L, 50L)) {
  amb_T <- sum(xv[xv < T_]) / s$n_cells
  cat(sprintf("  %6d %7.3f %15.2f %9.2f%s\n",
              T_, sum(xv >= T_) / s$n_cells, amb_T, (tot_per_cell - amb_T) / amb_T,
              if (T_ == THRESHOLD) "   <- THRESHOLD" else ""))
}


# ---- 3. what an ambient entry looks like -----------------------------------
# The step that converts "ambient nonzeros per cell" into "ambient UMIs per cell"
# assumes ambient counts are Poisson with a Dirichlet(1) composition. Under that
# assumption the per-entry count is geometric with mean m = A/G, so
#   E[count | count > 0] = 1 + m   and   P(count == 1 | count > 0) = 1/(1 + m).
# With m near zero both say ambient entries are essentially all single UMIs.
# If the measured values below sit well away from that, the assumption does not
# hold here and the model-based snr in section 4 should not be used.
A_model <- s$amb_nnz_per_cell / (1 - s$amb_nnz_per_cell / G)
m_hat   <- A_model / G
cat("\n[3] ambient entries\n")
cat(sprintf("  ambient nonzeros/cell     %8.2f\n", s$amb_nnz_per_cell))
cat(sprintf("  mean UMIs at an ambient nonzero   measured %6.3f   implied 1 + A/G = %6.3f\n",
            s$amb_count_mean, 1 + m_hat))
cat(sprintf("  fraction exactly 1 UMI           measured %6.3f   implied 1/(1 + A/G) = %6.3f\n",
            s$amb_frac_eq1, 1 / (1 + m_hat)))
# Ambient counts are confined to 1..threshold-1 by construction, so this ratio is
# of a doubly truncated distribution and has no reference value to be tested
# against. It is here to be compared against the same statistic on a simulated
# matrix, where the identical truncation applies.
cat(sprintf("  var/mean of ambient counts       measured %6.3f   (truncated to 1..%d)\n",
            s$amb_count_var_o_mean, THRESHOLD - 1L))
cat("  ambient count distribution:")
for (k in 1:min(4L, THRESHOLD - 1L)) {
  cat(sprintf("  %d:%.3f", k, mean(xv[xv < THRESHOLD] == k)))
}
if (THRESHOLD > 5L) cat(sprintf("  >=5:%.3f", mean(xv[xv < THRESHOLD] >= 5)))
cat("\n")


# ---- 4. snr, two ways ------------------------------------------------------
# snr in simulate_guidebender2 is the ratio of signal UMIs to noise UMIs:
# frac_signal = snr / (1 + snr) is the share of the total the signal gets
# (simulate.R:409). Both estimates below target that ratio.
A_direct   <- s$amb_umis_per_cell
snr_direct <- s$snr_observed
snr_model  <- (s$umis_per_cell_mean - A_model) / A_model
cat("\n[4] snr\n")
cat(sprintf("  ambient UMIs/cell   direct sum below threshold  %8.2f\n", A_direct))
cat(sprintf("                      inverted from nonzero count %8.2f\n", A_model))
cat(sprintf("  snr                 direct                      %8.2f\n", snr_direct))
cat(sprintf("                      model-based                 %8.2f\n", snr_model))
cat("  The direct estimate assumes only the threshold split. The model-based one\n")
cat("  additionally assumes section 3's ambient model; use the direct one unless\n")
cat("  section 3 says the threshold split is what is unreliable.\n")


# ---- 5. guide-level spread -------------------------------------------------
# guide_infection_alpha and endo_shape_sum are both left at 1, i.e. a symmetric
# Dirichlet(1) over guides, for library abundance and for ambient composition
# respectively. Each weight is then Beta(1, G-1): CV = sqrt((G-1)/(G+1)) and
# Gini -> 0.5. Those are properties of that distribution; the measured values
# are what this dataset actually has.
cat("\n[5] guide-level spread\n")
cat(sprintf("  Dirichlet(1) reference at G = %d:  CV = %.3f,  Gini -> 0.500\n",
            G, sqrt((G - 1) / (G + 1))))
cat(sprintf("  ambient UMIs per guide        CV = %.3f,  Gini = %.3f\n",
            s$guide_amb_cv, s$guide_amb_gini))
cat(sprintf("  perturbed cells per guide     CV = %.3f,  Gini = %.3f\n",
            s$guide_pert_cv, s$guide_pert_gini))
cat(sprintf("\n  per-cell signal CV %.3f, per-cell ambient CV %.3f, cor(log sig, log amb) %.3f\n",
            s$cell_sig_cv, s$cell_amb_cv, s$cor_log_sig_amb))
cat(sprintf("  var/mean at perturbed entries %.2f  (mean count there %.1f)\n",
            s$pert_count_var_o_mean, s$pert_count_mean))


# ---- 6. machine-readable row -----------------------------------------------
s$snr_direct <- snr_direct
s$snr_model  <- snr_model
s$A_model    <- A_model
out_fp <- file.path(here, sprintf("%s_real_stats.csv", DATASET))
write.csv(s, out_fp, row.names = FALSE)
cat("\nstats row ->", out_fp, "\n")
cat("Compare against a simulated matrix with guide_matrix_stats() from grna-sim-utils.R.\n\n")
