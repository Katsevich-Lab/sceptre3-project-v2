#!/usr/bin/env Rscript
# =============================================================================
# Turn the real-data characterization into a GRID OF REGIMES of Model B -- one
# calibrated regime per real dataset (barnyard excluded -- doublet-heavy, so it
# lives in Model A where the real ambient is used directly).
#
# Per-guide FIT design: each sim guide IS the per-guide fit of a real guide.
# From real_perguide.rds we sample G triples (mu, theta, signal_frac) jointly
# from the dataset's per-guide fits.  Then we calibrate the SHARED ambient
# (kappa_bar, phi_noise, alpha_amb) so the simulated below-per-guide-valley
# distribution matches the real one, measured the IDENTICAL way.
#
# Calibration objective: |sim_pct3 - target_pct3| + 6 * |log1p(sim_p99) -
# log1p(target_p99)|.  p99 (not max) is used because max is sample-size
# dependent (the calibration / build C mismatch is gone now anyway, but p99
# is more stable across seeds).  Each grid point is averaged over multiple
# seeds to reduce single-seed luck.  build_modelB() args are shared between
# the calibrate call and the final build (factored into BUILD_ARGS) so the
# calibrated tuple is for the SAME DGP that gets shipped.
#
# Output: results/sim_framework/regimes.csv (params + target/sim diagnostics),
#         the regime datasets in datasets/B__regime__*/, and manifest_modelB.csv.
# =============================================================================
suppressPackageStartupMessages({library(Matrix)})
source(file.path(getwd(), "scripts", "sim_modelB.R"))   # build_modelB(), .gen_cells()
OUT <- SIMFW()
ch <- read.csv(file.path(OUT, "real_characterization.csv"))

MIN_NZ <- 30   # min nonzero cells to fit a guide (shared with sim_characterize.R)

# below-per-guide-valley distribution of a matrix (SAME estimator as the real
# measurement -- see sim_characterize.R).  Returns pct>=3 and the 99th percentile
# of below-valley counts (size-invariant proxy for the ambient tail).
below_valley <- function(counts, ng = 150, seed = 11) {
  gmr <- as(counts, "RsparseMatrix"); below <- c()
  set.seed(seed)
  gi <- if (nrow(gmr) > ng) sort(sample(nrow(gmr), ng)) else seq_len(nrow(gmr))
  for (g in gi) { cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < MIN_NZ) next
    v <- smoothed_valley_threshold(cv); if (isTRUE(v$ok)) below <- c(below, cv[cv > 0 & cv < v$t]) }
  if (!length(below)) return(c(pct3 = 0, p99 = 0))
  c(pct3 = 100 * mean(below >= 3), p99 = as.numeric(quantile(below, 0.99)))
}

# Shared build args (the SAME DGP for calibration and final shipping).
# d_sigma controls per-cell signal-capture spread; w_sigma controls per-guide
# capture spread.  We deliberately keep both small because the per-guide fit
# already encodes signal spread via (mu, theta); larger d_sigma would inject
# extra signal variance not present in the fit and bias signal histograms wide.
# CAVEAT: small d_sigma combined with independently-drawn kappa_c means the
# library margin (signal + ambient) is a worse proxy for ambient depth than
# the ambient-depth margin -- a structural advantage for the depth_fix method
# vs fishash that grows with MOI.  A follow-up sensitivity sweep over
# d_sigma in [0.05, 0.5] is recommended before the depth_fix claim is final.
# See SIMULATION_FRAMEWORK.md "Known caveats".
.d_sigma <- as.numeric(Sys.getenv("SIM_D_SIGMA", "0.05"))   # see CAVEAT above
.w_sigma <- as.numeric(Sys.getenv("SIM_W_SIGMA", "0.05"))
stopifnot(is.finite(.d_sigma), .d_sigma >= 0, is.finite(.w_sigma), .w_sigma >= 0)
BUILD_ARGS <- list(alpha_inf = 1.5, w_sigma = .w_sigma, d_sigma = .d_sigma, doublet_rate = 0)
cat(sprintf("BUILD_ARGS d_sigma=%g w_sigma=%g (override via SIM_D_SIGMA / SIM_W_SIGMA)\n",
            .d_sigma, .w_sigma))

# calibrate ambient (kappa_bar, phi_noise, alpha_amb) to a target (pct>=3, p99)
# at the SAME C, theta, mu_vec, rate_vec as the final build, averaged over
# multiple seeds to reduce single-seed variance.
calibrate <- function(mu_vec, rate_vec, target_pct3, target_p99, theta_sig, G, C,
                      seeds = c(7, 17, 23)) {
  grid <- expand.grid(kb = c(0.4, 1, 2, 4, 8, 16, 32, 64), phi = c(0, 0.7, 1.5, 3), alpha = c(0.3, 0.7))
  best <- NULL; berr <- Inf
  for (r in seq_len(nrow(grid))) {
    bvs <- vapply(seeds, function(s) {
      args <- c(list(G = G, C = C, mu_pert = mu_vec, pert_rate = rate_vec, theta_sig = theta_sig,
                     alpha_amb = grid$alpha[r], kappa_bar = grid$kb[r], phi_noise = grid$phi[r],
                     seed = s), BUILD_ARGS)
      d <- do.call(build_modelB, args); below_valley(d$counts, seed = 100 + s)
    }, numeric(2))
    bv <- rowMeans(bvs)
    err <- abs(bv["pct3"] - target_pct3) + 6 * abs(log1p(bv["p99"]) - log1p(target_p99))
    if (is.finite(err) && err < berr) { berr <- err
      best <- list(kappa_bar = grid$kb[r], phi_noise = grid$phi[r], alpha_amb = grid$alpha[r],
                   sim_pct3 = round(bv["pct3"], 1), sim_p99 = round(bv["p99"], 1)) }
  }
  best
}

pg <- readRDS(file.path(OUT, "real_perguide.rds"))   # per-dataset df: mu, theta, signal_frac
ch <- ch[!grepl("^barnyard", ch$dataset), ]          # barnyard -> Model A (real ambient handles doublets)
dir.create(DSDIR(), recursive = TRUE, showWarnings = FALSE)
G <- 150; C <- 6000

# Draw per-guide (mu, theta, signal_frac) jointly from a dataset's per-guide fits.
# DEFLATION: the real signal_frac = mean(cv >= valley) counts ambient leakage
# above the valley as signal -- worst on high-ambient regimes.  We deflate by
# the dataset's ambient-above-valley mass (1 - dataset signal_frac is the
# ambient mass; the fraction of ambient that exceeds the valley is small but
# nonzero).  We use a conservative geometric deflation: rate * (1 - ambient_leak)
# where ambient_leak is approximated by the dataset's amb_pct3 fraction.
draw_pg <- function(ds, fb_mu, fb_rate, amb_pct3) {
  leak <- pmin(pmax(amb_pct3 / 100, 0), 0.5)         # ambient mass above the "3" threshold
  deflate <- 1 - leak                                # conservative deflation factor
  d <- pg[[ds]]
  if (!is.null(d) && nrow(d) >= 5) { idx <- sample(nrow(d), G, replace = TRUE)
    list(mu = pmax(d$mu[idx], 1), theta = pmax(pmin(d$theta[idx], 1e3), 0.05),
         rate = pmin(pmax(d$signal_frac[idx] * deflate, 1e-4), 0.3)) }
  else list(mu = rep(max(fb_mu, 2), G), theta = rep(2, G),
            rate = rep(min(max(fb_rate * deflate, 1e-3), 0.3), G))
}

# Real p99 target from the same below_valley estimator on the real matrix (see
# sim_characterize.R: amb_pct3 / amb_max are stored, but we need p99 here for
# the size-invariant loss).  We derive p99 from amb_max as a fallback only.
real_p99 <- function(ds_amb_p99, ds_amb_max) if (is.finite(ds_amb_p99)) ds_amb_p99 else ds_amb_max

rows <- list(); reg <- list()
for (i in seq_len(nrow(ch))) {
  ds <- ch$dataset[i]
  set.seed(9000 + i); pgd <- draw_pg(ds, ch$right_mode[i], ch$signal_frac[i], ch$amb_pct3[i])
  mu_vec <- pgd$mu; theta_vec <- pgd$theta; rate_vec <- pgd$rate
  emoi <- round(G * mean(rate_vec), 1)                            # emergent per-cell MOI
  tgt3 <- ifelse(is.finite(ch$amb_pct3[i]), ch$amb_pct3[i], 0)
  # use amb_max as p99 surrogate if no amb_p99 in characterization (legacy)
  tgt99 <- if ("amb_p99" %in% names(ch) && is.finite(ch$amb_p99[i])) ch$amb_p99[i]
           else if (is.finite(ch$amb_max[i])) ch$amb_max[i] else 5
  cal <- calibrate(mu_vec, rate_vec, tgt3, tgt99, theta_vec, G = G, C = C)
  id <- paste0("B__regime__", ds); seed <- 9000 + i
  args <- c(list(G = G, C = C, mu_pert = mu_vec, pert_rate = rate_vec, theta_sig = theta_vec,
                 alpha_amb = cal$alpha_amb, kappa_bar = cal$kappa_bar, phi_noise = cal$phi_noise,
                 seed = seed), BUILD_ARGS)
  d <- do.call(build_modelB, args)
  meta <- list(model = "B", id = id, regime = ds, chemistry = NA_character_, purity = NA,
               mu_level = "cal", mu_pert = round(median(mu_vec)),
               mu_p10 = round(quantile(mu_vec, .1)), mu_p90 = round(quantile(mu_vec, .9)),
               signal_frac = signif(mean(rate_vec), 3), emergent_moi = emoi,
               moi_level = ifelse(emoi >= 5, "high", "low"), moi = emoi,
               separation = round(ch$separation[i], 2), kappa_bar = cal$kappa_bar,
               phi_noise = cal$phi_noise, alpha_amb = cal$alpha_amb,
               target_pct3 = tgt3, sim_pct3 = as.numeric(cal$sim_pct3),
               target_p99 = tgt99, sim_p99 = as.numeric(cal$sim_p99),
               n_guides = G, n_cells = C, seed = seed)
  save_dataset(id, d$counts, d$Z, meta)
  rows[[id]] <- as.data.frame(meta[c("model", "id", "regime", "chemistry", "purity", "mu_level",
                                     "mu_pert", "moi_level", "moi")])
  reg[[id]]  <- as.data.frame(meta[c("regime", "separation", "mu_pert", "mu_p10", "mu_p90",
                                     "signal_frac", "emergent_moi", "moi", "kappa_bar", "phi_noise",
                                     "alpha_amb", "target_pct3", "sim_pct3", "target_p99", "sim_p99")])
  cat(sprintf("%-30s mu=%4d [%4d-%4d] sig_frac=%.4f emoi=%-4.1f  amb%%>=3 t=%-4.1f s=%-4.1f  p99 t=%-4.0f s=%-4.0f\n",
              ds, round(median(mu_vec)), round(quantile(mu_vec, .1)), round(quantile(mu_vec, .9)),
              mean(rate_vec), emoi, tgt3, as.numeric(cal$sim_pct3), tgt99, as.numeric(cal$sim_p99)))
}
man <- do.call(rbind, rows); rownames(man) <- NULL
write.csv(man, file.path(OUT, "manifest_modelB.csv"), row.names = FALSE)
regtab <- do.call(rbind, reg); rownames(regtab) <- NULL
write.csv(regtab, file.path(OUT, "regimes.csv"), row.names = FALSE)
cat(sprintf("\nwrote %d data-derived Model B regimes + manifest_modelB.csv + regimes.csv\n", nrow(man)))
