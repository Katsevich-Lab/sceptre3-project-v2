#!/usr/bin/env Rscript
# Generate the computational-scaling datasets.
#
#   Rscript simulate-scaling-datasets.R              # every rung
#   Rscript simulate-scaling-datasets.R 1 2          # just rungs 1 and 2
#
# Needs the R_461 renv library (fishash + zellkonverter + Matrix):
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# Counts come from fishash::simulate_guidebender2. Every parameter either that
# function or simulate_guidebender exposes is named below, so the settings can be
# read without opening the package.
#
# THE RAY. Datasets step along a 1-D path in (n_guides, n_cells), not a grid,
# because most (G, N) pairs are not plausible screens. Counting the true
# (guide, cell) pairs down columns and across rows gives
#
#   perturbed cells per guide = MOI * N/G
#
# so holding MOI and that count fixed forces N/G constant. The ray's slope is the
# real dataset's N/G. MOI and N/G are different quantities; they coincide only
# when MOI = 1.
#
# The manifest records realized nnz, density, MOI and perturbed cells per guide
# for each rung. Read those rather than assuming the rungs sit on the ray.

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})
source("~/.Rprofile")                                  # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
dp    <- normalizePath(file.path(here, "..", "data-preprocessing"))
source(file.path(dp, "convert_odm_to_h5ad.R"))         # R_to_h5ad()
source(file.path(dp, "lib_make_guide_data.R"))         # write_h5ad_methods(), write_cleanser_method()

H5AD_METHODS <- c("crispat", "pertpy")                 # cleanser reads the .mtx
SEED         <- 1L


# ---- the regime -------------------------------------------------------------
# gasperini-like: high MOI, CROP-seq. The dataset id keeps "gasperini" in it
# because run_cleanser.py picks --cs/--dc by substring match, and CROP-seq (--cs)
# is the chemistry of the dataset this regime mimics.
REGIME <- "gasperini"

# Measured by measure-real-targets.R on the real gasperini matrix
# (13077 guides x 207324 cells): moi 31.64 at threshold 5, median 564 guide UMIs
# per cell, 502 perturbed cells per guide.
N_OVER_G       <- 15.854   # 207324 cells / 13077 guides
MOI_TARGET     <- 31.64    # mean DISTINCT guides per cell, held equal at every rung
COUNT_PER_CELL <- 564      # median guide UMIs per cell

SNR             <- 4        # authors' primary setting
FRAC_NOISE_ENDO <- 0.75     # authors' primary setting

# The ray. n_cells is the input and n_guides is derived from it, not the other way
# round: simulate_guidebender errors unless chunk_cells divides n_cells, so n_cells
# is the constrained one and rounding it distorts the ray. n_guides has no such
# constraint, so deriving it holds N/G to within a fraction of a percent.
N_CELLS <- c(3000L, 6000L, 12000L)


# ---- parameters held fixed at every rung ------------------------------------
PARAMS_FIXED <- list(
  # PINNED to the authors' values
  hurdle_prob           = 0.1,    # share of recovered cells carrying no guide
  chunk_cells           = 1000,   # memory bound. Also holds fixed the median
                                  # rescaling and the exogenous bulk profile,
                                  # which simulate_guidebender computes per chunk
                                  # rather than per dataset.
  Phi_cell              = 1,      # signal is NB(mean, dispersion 1)
  Phi_noise             = 0,      # no gamma mixing on ambient: the rgamma step is
                                  # skipped entirely, leaving it conditionally
                                  # Poisson (the paper's primary model)

  # PACKAGE DEFAULTS, named so they are visible rather than implicit.
  # None of these is varied anywhere in the fishash paper or in
  # grna-count-modeling; endo_shape_* move only in the varySignalNoiseCorr
  # scenario, and Phi_noise only in the overdispersed-ambient arm.
  guide_infection_alpha = 1,      # symmetric Dirichlet over library abundance
  d_sigma_guide         = 0.5,    # SD of log guide expression size factor
  d_sigma_cell          = 0.5,    # SD of log cellular depth
  d_sigma_drop          = 0.5,    # SD of log ambient depth
  eps_alpha             = 50,     # capture efficiency ~ Gamma(a, rate a), mean 1
  rho_sum               = 10,     # Beta concentration for the chimeric fraction
  endo_shape_flat       = 0,      # ambient composition tracks library composition
  endo_shape_sum        = 1,      # Dirichlet concentration on that composition
  use_median            = TRUE    # count_per_cell is a median, not a mean

  # NOT settable: simulate_guidebender2 derives d_mu_drop, d_mu_cell, rho_alpha
  # and rho_beta from snr / count_per_cell / frac_noise_endo.
)


# ---- rungs ------------------------------------------------------------------
hurdle         <- PARAMS_FIXED$hurdle_prob
pert_per_guide <- MOI_TARGET * N_OVER_G

chunk <- PARAMS_FIXED$chunk_cells
stopifnot(N_CELLS %% chunk == 0)
rungs <- data.frame(n_cells  = N_CELLS,
                    n_guides = round(N_CELLS / N_OVER_G))
rungs$dataset  <- sprintf("sim_%s_g%d_c%d", REGIME, rungs$n_guides, rungs$n_cells)
rungs$n_over_g <- round(rungs$n_cells / rungs$n_guides, 3)

# Ground truth is (infection count > 0), so a guide drawn twice into one cell
# counts once, putting the realized MOI below the lambda passed in. Guides are
# drawn from a Dirichlet(1) abundance vector rather than uniformly, so with
# p_g ~ Exp(1)/G the expected number of DISTINCT guides per cell is
#   moi = lambda / (1 + lambda/G) * (1 - hurdle_prob)
# Checked against simulation for G in [200, 3200]: within ~1%, running slightly
# high. (Assuming uniform abundance instead is 4-9% high.)
#
# Inverting it gives the lambda that hits MOI_TARGET at each G. lambda therefore
# VARIES BY RUNG -- it is only a knob, and holding the realized MOI equal is what
# makes the rungs the same experiment at different sizes. Requires
# G > MOI_TARGET/(1-hurdle_prob), or no lambda reaches the target.
stopifnot(rungs$n_guides > MOI_TARGET / (1 - hurdle))
rungs$lambda <- round(MOI_TARGET / ((1 - hurdle) - MOI_TARGET / rungs$n_guides), 3)
rungs$moi_expected <- round(
  rungs$lambda / (1 + rungs$lambda / rungs$n_guides) * (1 - hurdle), 2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) rungs <- rungs[as.integer(args), , drop = FALSE]

cat(sprintf("\nray target: n_cells/n_guides = %.3f, MOI = %.2f -> %.0f perturbed cells per guide\n\n",
            N_OVER_G, MOI_TARGET, pert_per_guide))
print(rungs, row.names = FALSE); cat("\n")


# ---- generate ---------------------------------------------------------------
out_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                      "guide_assignment", "input_data")
manifest <- list()

for (i in seq_len(nrow(rungs))) {
  G <- rungs$n_guides[i]; N <- rungs$n_cells[i]; ds_id <- rungs$dataset[i]
  ds_dir <- file.path(out_root, ds_id)
  cat(sprintf("=== %s (%d guides x %d cells) ===\n", ds_id, G, N))

  t0 <- Sys.time()
  set.seed(SEED)
  sim <- do.call(simulate_guidebender2, c(
    list(n_guides = G, n_cells = N, moi = rungs$lambda[i], snr = SNR,
         count_per_cell = COUNT_PER_CELL, frac_noise_endo = FRAC_NOISE_ENDO,
         return_sparse_only = TRUE),
    PARAMS_FIXED
  ))
  gen_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  counts <- assay(sim, "counts")
  truth  <- assay(sim, "ground_truth")
  colnames(counts) <- colnames(truth) <- paste0("CELL_", seq_len(N))   # pipeline convention
  rownames(counts) <- rownames(truth) <- paste0("grna_", seq_len(G))

  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  write_cleanser_method(counts, ds_dir)
  write_h5ad_methods(counts, ds_dir, H5AD_METHODS)
  saveRDS(truth, file.path(ds_dir, "true_pert_matrix.rds"))
  saveRDS(c(PARAMS_FIXED, list(regime = REGIME, n_guides = G, n_cells = N,
                               moi = rungs$lambda[i], moi_target = MOI_TARGET, snr = SNR,
                               count_per_cell = COUNT_PER_CELL,
                               frac_noise_endo = FRAC_NOISE_ENDO, seed = SEED)),
          file.path(ds_dir, "sim_params.rds"))

  nnz <- Matrix::nnzero(counts)
  row <- data.frame(
    dataset        = ds_id, n_guides = G, n_cells = N, nnz = nnz,
    nnz_per_cell   = round(nnz / N, 2),
    nnz_per_guide  = round(nnz / G, 1),
    zero_frac      = round(1 - nnz / (as.numeric(G) * N), 5),
    pert_per_guide = round(mean(Matrix::rowSums(truth)), 1),
    moi_realized   = round(mean(Matrix::colSums(truth)), 3),
    umis_per_cell  = median(Matrix::colSums(counts)),
    gen_min        = round(gen_min, 2)
  )
  print(row, row.names = FALSE); cat("\n")
  manifest[[length(manifest) + 1L]] <- row
  rm(sim, counts, truth); invisible(gc())
}

manifest <- do.call(rbind, manifest)
mf <- file.path(out_root, sprintf("sim_scaling_manifest_%s.csv", REGIME))
write.csv(manifest, mf, row.names = FALSE)
cat("=== done ===\nmanifest ->", mf, "\n")
print(manifest, row.names = FALSE)
