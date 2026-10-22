# Design and simulator call for the MOI x ambient statistical-benchmarking grid,
# shared by simulate-stat-grid.R (pilot, generation) and eda-stat-grid.R.
# Needs fishash, Matrix and SummarizedExperiment loaded by the caller.
#
# Axes:
#   MOI     - lambda, the rate of the zero-truncated Poisson for integrations per
#             cell. With signal and ambient set directly (below), lambda moves only
#             the number of guides per cell (co-occurrence).
#   signal  - s, expected signal UMIs per integration, at fixed ambient B. Lowering
#             s slides the perturbed component into the ambient tail (separation)
#             and, inherently, raises dropout (true entries with count < 2).
# Ambient is held at a level where every guide's ambient histogram keeps its mode
# at 0 and decays fast; raising B instead pushed popular guides' ambient off zero.
#
# Why simulate_guidebender, not the wrapper. simulate_guidebender2 sets the per-
# CELL signal/noise ratio from snr and divides it among the integrations, so at
# fixed snr raising lambda also weakens each integration; it also ties rho's Beta
# to snr (alpha + beta = 10), so the spread of exogenous noise across cells
# changes with the design. Here the simulator's inputs are solved directly.
#
# Derivation (R/simulate.R at fishash 9b6eabc). Per cell, before Poisson draws:
#   signal   mu_gn = J_gn d_g^guide (1 - rho_n) eps_n d_n^cell
#   endo     lam_gn = chi_g^a (1 - rho_n) eps_n d_n^drop        (chi^a sums to 1)
#   exo      summed over guides, exactly rho_n * Y_n, Y_n the cell's whole total
# With A = (1 - pi) lambda / (1 - exp(-lambda)) expected integrations per cell,
# x = 1 - frac_noise_endo, m = E[rho] = alpha / (alpha + beta), E[eps] = 1, and
# rho, eps, d^cell, d^drop, d^guide drawn independently:
#   signal per integration  (1 - m) exp(mu_c + sc^2/2 + sg^2/2) = s
#   endo per cell           (1 - m) exp(mu_d + sd^2/2)          = (1 - x) B
#   exo per cell            m E[Y] = x B, where E[Y] = (A s + (1 - x) B) / (1 - m)
#                           =>   m = x B / (A s + B)
# so d_mu_cell = log(s/(1-m)) - sc^2/2 - sg^2/2, d_mu_drop = log((1-x) B/(1-m)) - sd^2/2,
# and the expected total per cell is A s + B. Exogenous noise is a share of the
# whole library, so holding it at x B while signal grows forces m down; m is the
# same mean the wrapper uses. What is held fixed instead is rho's shape: alpha is
# fixed and beta = alpha (1 - m) / m, so for small m rho ~ m * Gamma(alpha, alpha)
# and its coefficient of variation stays between about 1.0 and 1.41 (alpha = 0.5).
# median_per_cell = NULL: no rescaling, so s and B hold in expectation by
# construction. s is exact per INTEGRATION; per true entry it is higher at large
# lambda (a guide integrated twice is one entry), and is reported.

# G and N: per-entry ambient is about B / G, and perturbed cells per guide about
# N * A / G (~100, 280, 2,500 at lambda 0.3, 3, 30; real screens ~260-480). N is
# fixed across MOI, so prevalence rising with MOI is part of what MOI means.
STAT_GRID <- list(
  n_guides               = 500L,
  n_cells                = 50000L,
  hurdle_prob            = 0.1,
  signal_per_integration = c(80, 20, 5),        # s; 20 ~ Gasperini's ~18 UMIs per called guide
  frac_noise_endo        = 0.75,                # x = 0.25 of ambient is exogenous
  rho_alpha              = 0.5,                 # fixed shape of rho (package default)
  moi                    = c(0.3, 3, 30, 10),   # lambda; 10 APPENDED (not sorted) so the
                                                # seeds of the first 27 datasets, which
                                                # follow row order, are unchanged
  ambient_per_cell       = 25,                  # B, fixed
  n_reps                 = 3L,
  base_seed              = 20262000L,           # dataset k uses base_seed + k
  methods                = c("fishash", "fishashplus", "crispat", "pertpy", "cleanser"),
  # A precision stress test: ambient overdispersed per entry (Phi_noise), at
  # replicate 1 of every (moi, s) point. Each reuses its Phi_noise = 0 counterpart's
  # seed, so guides, cells, integrations, truth and signal rates are identical (all
  # drawn before the noise) and only the ambient differs: a paired comparison.
  phi_noise_variant      = list(Phi_noise = 1, reps = 1L, prefix = "stat_grid_phi1")
)

# The remaining simulate_guidebender arguments, listed so none is implicit. All
# are package defaults except d_sigma_drop (default 0.5): at 1.0 droplets vary more
# in how much ambient they carry, which lengthens each guide's ambient tail while
# keeping its mode at 0 (bypass_inputs() corrects d_mu_drop, so E[ambient] stays
# B). Two more differ from the fishash authors' runs, which went through the
# wrapper with chunk_cells = 1000:
#   median_per_cell = NULL - no rescaling; see the derivation above.
#   chunk_cells     = NULL - one chunk. The exogenous rate of guide g in cell n is
#     R_g C_n / S * rho_n / (1 - rho_n), with R_g guide g's row total over the
#     chunk. Over 1000 cells, R_g's signal part is the few integrations of g that
#     fall in the chunk, so each guide's exogenous level jumps from chunk to chunk
#     and the ambient guide profile is not shared across cells. Over all cells it
#     averages every integration of g: one smooth, library-wide profile. The
#     helper then holds dense G x n_cells rate matrices (about 32 MB each at
#     200 x 20,000).
STAT_GRID_PARAMS <- list(
  guide_infection_alpha = 1,
  d_sigma_guide         = 0.5,
  endo_shape_sum        = 1,
  endo_shape_flat       = 0,
  d_sigma_drop          = 1.0,
  d_sigma_cell          = 0.5,
  eps_alpha             = 50,
  Phi_cell              = 1,
  Phi_noise             = 0,
  return_sparse_only    = TRUE,
  median_per_cell       = NULL,
  chunk_cells           = NULL
)

#' Expected integrations per cell, over all cells.
avg_integrations <- function(moi, hurdle_prob) (1 - hurdle_prob) * moi / (1 - exp(-moi))

#' simulate_guidebender inputs that give signal s per integration and ambient B
#' per cell, a fraction x of it exogenous, with rho's alpha fixed.
bypass_inputs <- function(moi, s, B) {
  P <- STAT_GRID_PARAMS
  x <- 1 - STAT_GRID$frac_noise_endo
  A <- avg_integrations(moi, STAT_GRID$hurdle_prob)
  m <- x * B / (A * s + B)
  stopifnot(m > 0, m < 1)
  list(
    d_mu_cell = log(s / (1 - m)) - P$d_sigma_cell^2 / 2 - P$d_sigma_guide^2 / 2,
    d_mu_drop = log((1 - x) * B / (1 - m)) - P$d_sigma_drop^2 / 2,
    rho_alpha = STAT_GRID$rho_alpha,
    rho_beta  = STAT_GRID$rho_alpha * (1 - m) / m,
    rho_mean  = m,
    expected_total_per_cell = A * s + B
  )
}

simulate_grid <- function(moi, s, B, seed, Phi_noise = STAT_GRID_PARAMS$Phi_noise) {
  bi <- bypass_inputs(moi, s, B)
  params <- STAT_GRID_PARAMS
  params$Phi_noise <- Phi_noise      # same mean, so bypass_inputs() is unaffected
  set.seed(seed)
  sim <- do.call(fishash::simulate_guidebender, c(
    list(n_guides    = STAT_GRID$n_guides,
         n_cells     = STAT_GRID$n_cells,
         moi         = moi,
         hurdle_prob = STAT_GRID$hurdle_prob,
         d_mu_drop   = bi$d_mu_drop,
         d_mu_cell   = bi$d_mu_cell,
         rho_alpha   = bi$rho_alpha,
         rho_beta    = bi$rho_beta),
    params
  ))
  list(sim = sim, inputs = bi)
}

#' Realized levels of one simulated dataset.
realized <- function(sim) {
  counts <- assay(sim, "counts"); truth <- assay(sim, "ground_truth")
  signal <- assay(sim, "counts_signal")
  totals <- colSums(counts); true_cnt <- counts[truth]
  data.frame(
    true_entries_per_cell     = sum(truth) / ncol(counts),
    signal_per_true_entry     = sum(signal) / sum(truth),
    ambient_per_cell_realized = mean(totals - colSums(signal)),
    umis_per_cell_mean        = mean(totals),
    umis_per_cell_median      = median(totals),
    frac_true_count_lt2       = mean(true_cnt < 2)     # below min_count = 2
  )
}

name_dims <- function(m) {
  colnames(m) <- paste0("CELL_", seq_len(ncol(m)))       # pipeline convention
  rownames(m) <- paste0("grna_", seq_len(nrow(m)))
  m
}
