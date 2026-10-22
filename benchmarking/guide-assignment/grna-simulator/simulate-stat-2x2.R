#!/usr/bin/env Rscript
# Generate the 2x2 statistical-benchmarking datasets: MOI x signal per integration.
#
#   sbatch run-in-image.sh simulate-stat-2x2.R
#
# The fishash authors' varyMOI scenario holds snr and count_per_cell fixed. The
# wrapper then fixes each cell's expected signal UMIs (snr/(1+snr) * count_per_cell)
# and divides it among the cell's integrations (the -log(avg_infections) term in
# simulate_guidebender2), so raising MOI also weakens every integration: about
# 85 signal UMIs per integration at MOI 0.1, about 9 at MOI 10. This design moves
# the two factors separately. Each dataset is set by
#   m - MOI (the wrapper's `moi`)
#   s - expected signal UMIs per integration
#   B - expected ambient (endogenous + exogenous) UMIs per cell, held fixed
# and the wrapper's inputs are derived from them. With
#   A(m) = (1 - hurdle_prob) * m / (1 - exp(-m)),
# the expected integrations per cell over all cells,
#   count_per_cell = B + s * A(m)      (ambient plus the signal of every integration)
#   snr            = s * A(m) / B      (expected signal over expected ambient, per cell)
# These hold in expectation before the median rescaling, which scales each
# 1000-cell chunk so its MEDIAN total equals count_per_cell. Totals are right-
# skewed, so the rescaling inflates every rate, by a factor that differs between
# design cells (about 1.2 to 1.8 when count_per_cell was set from the formula
# alone). The rescaling leaves the signal/ambient ratio, and so snr, in place,
# and every final rate is proportional to count_per_cell. So count_per_cell is
# calibrated per design cell: simulate at CAL_SEED, measure the signal per true
# entry, multiply count_per_cell by target/realized, and repeat until within
# CAL_TOL. The target is matched per true ENTRY, which is what the methods see; a
# guide integrated twice in one cell is one entry, so at MOI 10 this is a few
# percent above the signal per integration. Realized values are measured again on
# every generated dataset and written to the manifest.
#
# Everything else is the authors' high-gRNA setting (fishash_analysis
# simulate_numCells20k_medUmi100_snr4_endo75_varyNumGuides.R): 200 guides, 20,000
# cells, hurdle_prob 0.1, frac_noise_endo 0.75, chunk_cells 1000, and the package
# defaults. B = 20 is that setting's ambient (count_per_cell 100, snr 4). The s
# levels, 9 and 85, are the per-integration signal at the two ends of their MOI
# sweep.
#
# Known side effect, accepted for now: the wrapper sets rho_n ~ Beta(phi_x * 10,
# (1 - phi_x) * 10) with phi_x = (1 - frac_noise_endo) / (1 + snr). The expected
# exogenous UMIs per cell stay at (1 - frac_noise_endo) * B, but at large snr the
# Beta concentrates near 0, so exogenous noise falls on fewer cells.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/input_data/:
#   stat2x2_moi<m>_s<s>_rep<NN>/   fishash and fishashplus inputs (see
#                                  SCALING_METHODS in grna-sim-utils.R),
#                                  true_pert_matrix.rds, sim_params.rds
#   stat2x2_manifest.csv           targets, wrapper inputs and realized values

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
dp    <- normalizePath(file.path(here, "..", "data-preprocessing"))
source(file.path(here, "grna-sim-utils.R"))             # write_method_inputs()
source(file.path(dp, "convert_odm_to_h5ad.R"))           # R_to_h5ad()
source(file.path(dp, "lib_make_guide_data.R"))           # write_cleanser_method()

# ---- design ------------------------------------------------------------------
STAT2X2 <- list(
  n_guides               = 200L,
  n_cells                = 20000L,
  hurdle_prob            = 0.1,
  ambient_per_cell       = 20,          # B
  moi                    = c(0.3, 10),  # m
  signal_per_integration = c(9, 85),    # s
  n_reps                 = 10L,
  base_seed              = 20261008L,   # dataset k uses base_seed + k
  methods                = c("fishash", "fishashplus")
)

CAL_SEED       <- 11L    # calibration runs only; no dataset uses it
CAL_TOL        <- 0.01   # stop within 1% of the target signal per true entry
CAL_MAX_PASSES <- 4L

# Every other simulate_guidebender2 / simulate_guidebender argument, listed so
# none is implicit. All are package defaults except frac_noise_endo (no default;
# the authors' high-gRNA value) and chunk_cells (default NULL; the authors' 1000).
STAT2X2_PARAMS <- list(
  guide_infection_alpha = 1,
  frac_noise_endo       = 0.75,
  Phi_cell              = 1,
  Phi_noise             = 0,
  rho_sum               = 10,
  eps_alpha             = 50,
  d_sigma_guide         = 0.5,
  d_sigma_cell          = 0.5,
  d_sigma_drop          = 0.5,
  endo_shape_sum        = 1,
  endo_shape_flat       = 0,
  use_median            = TRUE,
  chunk_cells           = 1000
)

#' Expected integrations per cell, over all cells (the wrapper's avg_infections).
avg_integrations <- function(moi, hurdle_prob) (1 - hurdle_prob) * moi / (1 - exp(-moi))

#' The wrapper inputs that put s signal UMIs on each integration over B ambient.
wrapper_inputs <- function(moi, s, B, hurdle_prob) {
  A <- avg_integrations(moi, hurdle_prob)
  list(count_per_cell = B + s * A, snr = s * A / B)
}

#' The single call into the simulator, used by calibration and generation alike.
simulate_stat <- function(moi, snr, count_per_cell, seed) {
  set.seed(seed)
  do.call(fishash::simulate_guidebender2, c(
    list(n_guides       = STAT2X2$n_guides,
         n_cells        = STAT2X2$n_cells,
         moi            = moi,
         hurdle_prob    = STAT2X2$hurdle_prob,
         snr            = snr,
         count_per_cell = count_per_cell,
         return_sparse_only = TRUE),
    STAT2X2_PARAMS
  ))
}

signal_per_true_entry <- function(sim) {
  sum(assay(sim, "counts_signal")) / sum(assay(sim, "ground_truth"))
}

# ---- calibrate count_per_cell per design cell --------------------------------
design <- expand.grid(s = STAT2X2$signal_per_integration, moi = STAT2X2$moi)
design$count_per_cell_nominal <- NA_real_
design$snr <- NA_real_
design$count_per_cell <- NA_real_
for (d in seq_len(nrow(design))) {
  wi <- wrapper_inputs(design$moi[d], design$s[d], STAT2X2$ambient_per_cell,
                       STAT2X2$hurdle_prob)
  C <- wi$count_per_cell
  for (pass in seq_len(CAL_MAX_PASSES)) {
    realized <- signal_per_true_entry(simulate_stat(design$moi[d], wi$snr, C, CAL_SEED))
    cat(sprintf("calibrate moi %g, s %g, pass %d: count_per_cell %.3f -> signal per true entry %.3f\n",
                design$moi[d], design$s[d], pass, C, realized))
    if (abs(realized / design$s[d] - 1) <= CAL_TOL) break
    if (pass == CAL_MAX_PASSES)
      stop(sprintf("calibration did not reach %g%% of s = %g at moi %g",
                   100 * CAL_TOL, design$s[d], design$moi[d]))
    C <- C * design$s[d] / realized
  }
  design$count_per_cell_nominal[d] <- wi$count_per_cell
  design$snr[d]                    <- wi$snr
  design$count_per_cell[d]         <- C
}
design$cal_factor <- design$count_per_cell / design$count_per_cell_nominal
cat("\n=== calibrated design ===\n")
print(signif(design, 5), row.names = FALSE)

grid <- expand.grid(rep  = seq_len(STAT2X2$n_reps),
                    s    = STAT2X2$signal_per_integration,
                    moi  = STAT2X2$moi)
grid$dataset <- sprintf("stat2x2_moi%g_s%g_rep%02d", grid$moi, grid$s, grid$rep)
grid$seed    <- STAT2X2$base_seed + seq_len(nrow(grid))

out_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data")

# ---- generate ----------------------------------------------------------------
rows <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  g  <- grid[k, ]
  wi <- design[design$moi == g$moi & design$s == g$s, ]
  stopifnot(nrow(wi) == 1)
  cat(sprintf("\n=== %s: count_per_cell %.2f, snr %.4f, seed %d ===\n",
              g$dataset, wi$count_per_cell, wi$snr, g$seed))

  sim <- simulate_stat(g$moi, wi$snr, wi$count_per_cell, g$seed)
  counts <- assay(sim, "counts")
  truth  <- assay(sim, "ground_truth")
  signal <- assay(sim, "counts_signal")
  colnames(counts) <- colnames(truth) <- paste0("CELL_", seq_len(STAT2X2$n_cells))   # pipeline convention
  rownames(counts) <- rownames(truth) <- paste0("grna_", seq_len(STAT2X2$n_guides))

  ds_dir <- file.path(out_root, g$dataset)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  write_method_inputs(counts, ds_dir, STAT2X2$methods)
  saveRDS(truth, file.path(ds_dir, "true_pert_matrix.rds"))
  saveRDS(c(STAT2X2_PARAMS,
            list(design = "stat2x2", n_guides = STAT2X2$n_guides, n_cells = STAT2X2$n_cells,
                 hurdle_prob = STAT2X2$hurdle_prob, moi = g$moi,
                 signal_per_integration = g$s, ambient_per_cell = STAT2X2$ambient_per_cell,
                 count_per_cell = wi$count_per_cell,
                 count_per_cell_nominal = wi$count_per_cell_nominal,
                 cal_factor = wi$cal_factor, cal_seed = CAL_SEED, snr = wi$snr,
                 rep = g$rep, seed = g$seed, methods = STAT2X2$methods)),
          file.path(ds_dir, "sim_params.rds"))

  # Realized values. Signal is per integrated ENTRY (a guide integrated twice in
  # a cell is one entry), so it can sit slightly above the per-integration target.
  totals    <- colSums(counts)
  true_cnt  <- counts[truth]
  rows[[k]] <- data.frame(
    dataset                     = g$dataset,
    moi                         = g$moi,
    signal_per_integration      = g$s,
    ambient_per_cell            = STAT2X2$ambient_per_cell,
    count_per_cell_nominal      = wi$count_per_cell_nominal,
    cal_factor                  = wi$cal_factor,
    count_per_cell              = wi$count_per_cell,
    snr                         = wi$snr,
    seed                        = g$seed,
    integrations_per_cell_truth = sum(truth) / STAT2X2$n_cells,
    signal_per_true_entry       = sum(signal) / sum(truth),
    ambient_per_cell_realized   = mean(totals - colSums(signal)),
    umis_per_cell_median        = median(totals),
    umis_per_cell_mean          = mean(totals),
    frac_true_count0            = mean(true_cnt == 0),
    frac_true_count_lt2         = mean(true_cnt < 2)    # below min_count = 2 for both methods
  )
  print(rows[[k]][, -1], row.names = FALSE)
  rm(sim, counts, truth, signal); invisible(gc())
}

manifest <- do.call(rbind, rows)
mf <- file.path(out_root, "stat2x2_manifest.csv")
write.csv(manifest, mf, row.names = FALSE)
cat("\nmanifest ->", mf, "\n")

# Realized values by design cell, averaged over replicates.
agg <- aggregate(manifest[, c("integrations_per_cell_truth", "signal_per_true_entry",
                              "ambient_per_cell_realized", "umis_per_cell_median",
                              "frac_true_count0", "frac_true_count_lt2")],
                 by = manifest[, c("moi", "signal_per_integration")], FUN = mean)
print(signif(agg, 4), row.names = FALSE)
