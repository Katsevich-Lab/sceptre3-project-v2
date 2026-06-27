#!/usr/bin/env Rscript
# =============================================================================
# Turn the real-data characterization into a GRID OF REGIMES of Model B -- one
# calibrated regime per real dataset.  Each regime sets:
#   mu_pert  = the dataset's right (signal) mode        (signal level)
#   moi      = the dataset's MOI (capped for tractability)
#   ambient  = kappa_bar, phi_noise tuned so the SIMULATED below-per-guide-valley
#              distribution (% >= 3 and max) matches the real one -- this is the
#              "ambiguous middle" the methods must handle, and it is calibrated
#              the SAME way it was measured on the real data.
# Output: results/sim_framework/regimes.csv (params + target vs achieved),
#         17 Model B datasets, and manifest_modelB.csv (the new data-driven Model B).
# =============================================================================
suppressPackageStartupMessages({library(Matrix)})
MODELB_FUNCS_ONLY <- TRUE
source(file.path(getwd(), "scripts", "sim_modelB.R"))   # build_modelB(), .gen_cells()
OUT <- SIMFW()
ch <- read.csv(file.path(OUT, "real_characterization.csv"))

# below-per-guide-valley distribution of a matrix (SAME estimator as the real measurement)
below_valley <- function(counts, ng = 150) {
  gmr <- as(counts, "RsparseMatrix"); below <- c()
  gi <- if (nrow(gmr) > ng) sort(sample(nrow(gmr), ng)) else seq_len(nrow(gmr))
  for (g in gi) { cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < 20) next
    v <- smoothed_valley_threshold(cv); if (isTRUE(v$ok)) below <- c(below, cv[cv > 0 & cv < v$t]) }
  if (!length(below)) return(c(pct3 = 0, max = 0))
  c(pct3 = 100 * mean(below >= 3), max = max(below))
}

# calibrate ambient (kappa_bar depth, phi_noise overdispersion, alpha_amb pool
# concentration) to a target below-valley (% >= 3, max), with the signal (per-guide
# level mu_vec + per-guide perturbation rate) held fixed.
calibrate <- function(mu_vec, rate_vec, target_pct3, target_max, theta_sig, G = 150, C = 2000, seed = 7) {
  grid <- expand.grid(kb = c(0.4, 1, 2, 4, 8, 16, 32, 64), phi = c(0, 0.7, 1.5, 3), alpha = c(0.3, 0.7))
  best <- NULL; berr <- Inf
  for (r in seq_len(nrow(grid))) {
    d <- build_modelB(G = G, C = C, mu_pert = mu_vec, pert_rate = rate_vec, theta_sig = theta_sig,
                      alpha_inf = 1.5, alpha_amb = grid$alpha[r], kappa_bar = grid$kb[r],
                      phi_noise = grid$phi[r], doublet_rate = 0, seed = seed)
    bv <- below_valley(d$counts)
    err <- abs(bv["pct3"] - target_pct3) + 6 * abs(log1p(bv["max"]) - log1p(target_max))
    if (is.finite(err) && err < berr) { berr <- err
      best <- list(kappa_bar = grid$kb[r], phi_noise = grid$phi[r], alpha_amb = grid$alpha[r],
                   sim_pct3 = round(bv["pct3"],1), sim_max = bv["max"]) }
  }
  best
}

pg <- readRDS(file.path(OUT, "real_perguide.rds"))   # per-dataset df: mu, theta, signal_frac
ch <- ch[!grepl("^barnyard", ch$dataset), ]          # barnyard -> Model A (doublet-heavy; hard to fit)
dir.create(DSDIR(), recursive = TRUE, showWarnings = FALSE)
G <- 150; C <- 6000
# Each sim guide is the PARAMETRIC FIT of a real guide: its (signal mean mu, NB
# dispersion theta, perturbed fraction) are sampled jointly from the dataset's
# per-guide fits, so each guide's histogram matches a real guide's by construction.
# Minimal d/w-sigma so signal ~ NB(mu, theta) directly (the fit already has spread).
draw_pg <- function(ds, fb_mu, fb_rate) {
  d <- pg[[ds]]
  if (!is.null(d) && nrow(d) >= 5) { idx <- sample(nrow(d), G, replace = TRUE)
    list(mu = pmax(d$mu[idx], 1), theta = pmax(pmin(d$theta[idx], 1e3), 0.05),
         rate = pmin(pmax(d$signal_frac[idx], 1e-4), 0.3)) }
  else list(mu = rep(max(fb_mu, 2), G), theta = rep(2, G), rate = rep(min(max(fb_rate, 1e-3), 0.3), G))
}
rows <- list(); reg <- list()
for (i in seq_len(nrow(ch))) {
  ds <- ch$dataset[i]
  set.seed(9000 + i); pgd <- draw_pg(ds, ch$right_mode[i], ch$signal_frac[i])
  mu_vec <- pgd$mu; theta_vec <- pgd$theta; rate_vec <- pgd$rate
  emoi <- round(G * mean(rate_vec), 1)                            # emergent per-cell MOI
  tgt3 <- ifelse(is.finite(ch$amb_pct3[i]), ch$amb_pct3[i], 0)
  tgtm <- ifelse(is.finite(ch$amb_max[i]), ch$amb_max[i], 5)
  cal <- calibrate(mu_vec, rate_vec, tgt3, tgtm, theta_vec)
  id <- paste0("B__regime__", ds); seed <- 9000 + i
  d <- build_modelB(G = G, C = C, mu_pert = mu_vec, pert_rate = rate_vec, theta_sig = theta_vec,
                    alpha_inf = 1.5, w_sigma = 0.05, d_sigma = 0.05, alpha_amb = cal$alpha_amb,
                    kappa_bar = cal$kappa_bar, phi_noise = cal$phi_noise, doublet_rate = 0, seed = seed)
  meta <- list(model = "B", id = id, regime = ds, chemistry = ds, purity = NA,
               mu_level = "cal", mu_pert = round(median(mu_vec)),
               mu_p10 = round(quantile(mu_vec,.1)), mu_p90 = round(quantile(mu_vec,.9)),
               signal_frac = signif(mean(rate_vec), 3), emergent_moi = emoi,
               moi_level = ifelse(emoi >= 5, "high", "low"), moi = emoi,
               separation = round(ch$separation[i], 2), kappa_bar = cal$kappa_bar, phi_noise = cal$phi_noise,
               alpha_amb = cal$alpha_amb, target_pct3 = tgt3, sim_pct3 = as.numeric(cal$sim_pct3),
               n_guides = G, n_cells = C, seed = seed)
  save_dataset(id, d$counts, d$Z, meta)
  rows[[id]] <- as.data.frame(meta[c("model","id","regime","chemistry","purity","mu_level","mu_pert","moi_level","moi")])
  reg[[id]]  <- as.data.frame(meta[c("regime","separation","mu_pert","mu_p10","mu_p90","signal_frac","emergent_moi","kappa_bar","phi_noise","alpha_amb","target_pct3","sim_pct3")])
  cat(sprintf("%-30s mu=%4d [%4d-%4d] sig_frac=%.4f emoi=%-4.1f  amb%%>=3 target=%-4.1f sim=%-4.1f\n",
              ds, round(median(mu_vec)), round(quantile(mu_vec,.1)), round(quantile(mu_vec,.9)),
              mean(rate_vec), emoi, tgt3, as.numeric(cal$sim_pct3)))
}
man <- do.call(rbind, rows); rownames(man) <- NULL
write.csv(man, file.path(OUT, "manifest_modelB.csv"), row.names = FALSE)
regtab <- do.call(rbind, reg); rownames(regtab) <- NULL
write.csv(regtab, file.path(OUT, "regimes.csv"), row.names = FALSE)
cat("\nwrote", nrow(man), "data-derived Model B regimes + manifest_modelB.csv + regimes.csv\n")
