#!/usr/bin/env Rscript
# Tests for simulate-scaling-datasets.R.
#
#   Rscript tests/test-simulate-scaling.R      # from grna-simulator/
#
# Needs the R_461 renv library:
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# Covers four things:
#   1. round-trip -- the .mtx and .h5ad written for each method hold the same
#      matrix as the SummarizedExperiment they came from, in the right orientation
#   2. the ray -- generated data has the intended perturbed cells per guide,
#      and the manifest columns are mutually consistent
#   3. calibration -- realized median UMIs per cell matches count_per_cell, and
#      realized MOI matches the collision-adjusted prediction
#   4. plumbing -- every name in PARAMS_FIXED reaches the simulator rather than
#      being silently swallowed by simulate_guidebender2's `...`

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment); library(zellkonverter)
})

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
dp    <- normalizePath(file.path(here, "..", "..", "data-preprocessing"))
source(file.path(dp, "convert_odm_to_h5ad.R"))
source(file.path(dp, "lib_make_guide_data.R"))

check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}
near <- function(a, b, tol) abs(a - b) <= tol * abs(b)

# Same fixed parameters the generator uses.
PARAMS_FIXED <- list(
  hurdle_prob = 0.1, chunk_cells = 1000, Phi_cell = 1, Phi_noise = 0,
  guide_infection_alpha = 1, d_sigma_guide = 0.5, d_sigma_cell = 0.5,
  d_sigma_drop = 0.5, eps_alpha = 50, rho_sum = 10,
  endo_shape_flat = 0, endo_shape_sum = 1, use_median = TRUE
)

gen <- function(G, N, lambda, count_per_cell = 564, snr = 4, endo = 0.75,
                seed = 1L, extra = list()) {
  set.seed(seed)
  do.call(simulate_guidebender2, c(
    list(n_guides = G, n_cells = N, moi = lambda, snr = snr,
         count_per_cell = count_per_cell, frac_noise_endo = endo,
         return_sparse_only = TRUE),
    modifyList(PARAMS_FIXED, extra)
  ))
}


# ---- 1. round-trip ----------------------------------------------------------
cat("\n[1] written files match the source matrix\n")
G <- 20L; N <- 1000L
sim <- gen(G, N, lambda = 2)
counts <- assay(sim, "counts")
colnames(counts) <- paste0("CELL_", seq_len(N))
rownames(counts) <- paste0("grna_", seq_len(G))

tmp <- file.path(tempdir(), "rt"); dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
write_cleanser_method(counts, tmp)
write_h5ad_methods(counts, tmp, c("crispat", "pertpy"))

mtx <- Matrix::readMM(file.path(tmp, "cleanser", "grna_matrix.mtx"))
check(identical(dim(mtx), dim(counts)), "mtx dims are guides x cells")
check(Matrix::nnzero(mtx) == Matrix::nnzero(counts), "mtx preserves nnz")
check(all(as.matrix(mtx) == as.matrix(counts)), "mtx preserves every value")

for (m in c("crispat", "pertpy")) {
  sce <- readH5AD(file.path(tmp, m, "grna_matrix.h5ad"))
  a <- assay(sce, 1L)
  check(identical(dim(a), dim(counts)), sprintf("%s h5ad round-trips guides x cells", m))
  check(all(as.matrix(a) == as.matrix(counts)), sprintf("%s h5ad preserves every value", m))
}
h1 <- readH5AD(file.path(tmp, "crispat", "grna_matrix.h5ad"))
h2 <- readH5AD(file.path(tmp, "pertpy",  "grna_matrix.h5ad"))
check(all(as.matrix(assay(h1, 1L)) == as.matrix(assay(h2, 1L))),
      "crispat and pertpy receive identical matrices")

truth <- assay(sim, "ground_truth")
check(identical(dim(truth), dim(counts)), "ground truth matches counts dims")


# ---- 2. the ray -------------------------------------------------------------
cat("\n[2] ray arithmetic holds in generated data\n")
G <- 400L; N <- 6000L; LAMBDA <- 35.1549
sim <- gen(G, N, LAMBDA)
truth <- assay(sim, "ground_truth")

pert_per_guide <- mean(Matrix::rowSums(truth))
moi_realized   <- mean(Matrix::colSums(truth))

# Both are sum(truth) rescaled, so this must hold exactly -- it guards the
# manifest columns against being computed off the wrong margin.
check(near(pert_per_guide * G / N, moi_realized, 1e-10),
      "pert_per_guide * G / N == moi_realized")

# perturbed cells per guide = MOI * N/G, the identity the ray is built on.
check(near(pert_per_guide, moi_realized * N / G, 1e-10),
      "pert_per_guide == moi_realized * N/G")

# It is NOT N/G alone -- at this MOI the two differ by ~30x.
check(pert_per_guide > 10 * (N / G),
      "pert_per_guide is far from N/G (guards the MOI factor)")


# ---- 3. calibration ---------------------------------------------------------
cat("\n[3] realized values match what the parameters asked for\n")
counts <- assay(sim, "counts")
check(near(median(Matrix::colSums(counts)), 564, 0.02),
      "median UMIs per cell == count_per_cell (within 2%)")

# Ground truth is (infection count > 0), so a guide drawn twice into one cell
# counts once. Guides come from a Dirichlet(1) abundance vector, not a uniform
# one, so with p_g ~ Exp(1)/G the expected number of DISTINCT guides per cell is
# lambda/(1 + lambda/G). Assuming uniform abundance instead is 4-9% high.
h <- PARAMS_FIXED$hurdle_prob
moi_expected <- LAMBDA / (1 + LAMBDA / G) * (1 - h)
check(near(moi_realized, moi_expected, 0.02),
      sprintf("realized MOI %.2f matches collision formula %.2f (within 2%%)",
              moi_realized, moi_expected))
check(moi_expected < LAMBDA * (1 - h),
      "collisions put realized MOI below the requested MOI")
check(moi_expected < G * (1 - exp(-LAMBDA / G)) * (1 - h),
      "Dirichlet abundance gives fewer distinct guides than uniform would")

check(inherits(try(gen(50L, 1500L, 2, extra = list(chunk_cells = 1000)), silent = TRUE),
               "try-error"),
      "chunk_cells must divide n_cells (errors otherwise)")


# ---- 3a. lambda compensation holds realized MOI constant --------------------
cat("\n[3a] per-rung lambda keeps realized MOI equal across G\n")
MOI_TARGET <- 31.64
lambda_for <- function(G, moi = MOI_TARGET, hurdle = h) moi / ((1 - hurdle) - moi / G)

for (Gi in c(189L, 378L, 756L)) {
  lam <- lambda_for(Gi)
  obs <- mean(Matrix::colSums(assay(gen(Gi, 6000L, lam), "ground_truth")))
  check(near(obs, MOI_TARGET, 0.02),
        sprintf("G=%d: lambda %.2f gives realized MOI %.2f (target %.2f)",
                Gi, lam, obs, MOI_TARGET))
}
# Compensation must move lambda, or it is doing nothing.
check(lambda_for(189L) > lambda_for(756L) * 1.05,
      "smaller G needs a larger lambda to reach the same MOI")
# No lambda can reach a MOI above the number of guides that survive the hurdle.
check(lambda_for(400L) > 0 && lambda_for(30L) < 0,
      "target is unreachable when G < MOI/(1 - hurdle_prob)")


# ---- 3b. statistical properties ---------------------------------------------
cat("\n[3b] statistical properties of the generated data\n")
signal <- assay(sim, "counts_signal")

# Total UMIs split signal vs noise in the ratio snr/(1+snr).
check(near(sum(signal) / sum(counts), 4 / 5, 0.02),
      "signal share of all UMIs == snr/(1+snr) (within 2%)")

# hurdle_prob is the share of recovered cells that received no infection at all.
check(near(mean(Matrix::colSums(truth) == 0), h, 0.10),
      "fraction of cells with no perturbation == hurdle_prob (within 10%)")

# Perturbations per guide inherit the spread of the Dirichlet(1) abundance
# vector. With p_g ~ Exp(1)/G and rate r_g = 1-exp(-lambda*p_g), the CV works out
# near 0.92 at this lambda/G; uniform abundance would put it near 0.04, so this
# separates the two by ~20x rather than testing a fitted constant.
cv_pert <- sd(Matrix::rowSums(truth)) / mean(Matrix::rowSums(truth))
check(cv_pert > 0.8 && cv_pert < 1.1,
      sprintf("CV of perturbations per guide is %.3f, consistent with Dirichlet(1)", cv_pert))

# The assignment problem must be solvable: perturbed entries carry far more UMIs.
mean_pert   <- sum(counts * truth) / sum(truth)
mean_unpert <- sum(counts * !truth) / sum(!truth)
check(mean_pert / mean_unpert > 20,
      sprintf("perturbed entries carry %.0fx the UMIs of unperturbed ones",
              mean_pert / mean_unpert))

# count_per_cell is a median and the depth factors are lognormal, so the mean
# total sits above it.
check(mean(Matrix::colSums(counts)) > median(Matrix::colSums(counts)),
      "mean UMIs per cell exceeds the median (right-skewed depth)")


# ---- 4. plumbing ------------------------------------------------------------
cat("\n[4] PARAMS_FIXED reaches the simulator\n")
# simulate_guidebender has no `...`, so a misspelled name forwarded through
# simulate_guidebender2's `...` raises rather than being ignored.
check(inherits(try(gen(50L, 1000L, 2, extra = list(Phi_nosie = 1)), silent = TRUE),
               "try-error"),
      "a misspelled parameter name errors rather than being silently dropped")

# Each of these must change the data, or it is not being applied.
base <- assay(gen(50L, 1000L, 2), "counts")
for (p in list(list(Phi_noise = 1), list(guide_infection_alpha = 100),
               list(d_sigma_guide = 2), list(eps_alpha = 1),
               list(endo_shape_flat = 1), list(endo_shape_sum = 1e6),
               list(rho_sum = 1000), list(Phi_cell = 0.01))) {
  alt <- assay(gen(50L, 1000L, 2, extra = p), "counts")
  check(!identical(as.matrix(base), as.matrix(alt)),
        sprintf("%s changes the generated counts", names(p)))
}

cat("\nAll tests passed.\n")
