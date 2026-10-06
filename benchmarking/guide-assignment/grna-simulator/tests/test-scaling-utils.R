#!/usr/bin/env Rscript
# Tests for the "GuideBender scaling datasets" section of grna-sim-utils.R.
#
#   Rscript tests/test-scaling-utils.R      # from grna-simulator/
#
# Needs the R_461 renv library:
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment); library(zellkonverter)
})

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
sim_dir <- normalizePath(file.path(here, ".."))
dp      <- normalizePath(file.path(sim_dir, "..", "data-preprocessing"))
source(file.path(sim_dir, "grna-sim-utils.R"))
source(file.path(dp, "convert_odm_to_h5ad.R"))
source(file.path(dp, "lib_make_guide_data.R"))

check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}
near <- function(a, b, tol) abs(a - b) <= tol * abs(b)


# ---- 1. sourcing attaches nothing ------------------------------------------
cat("\n[1] sourcing grna-sim-utils.R attaches no packages\n")
out <- system2(file.path(R.home("bin"), "Rscript"),
               c("-e", shQuote(sprintf(
                 "b <- search(); source('%s'); cat(identical(b, search()))",
                 file.path(sim_dir, "grna-sim-utils.R")))),
               stdout = TRUE, stderr = FALSE)
check(tail(out, 1) == "TRUE", "search path unchanged after sourcing")


# ---- 2. alpha_from_cv ------------------------------------------------------
cat("\n[2] alpha_from_cv recovers the Dirichlet concentration\n")
set.seed(3)
G <- 5000L; mean_count <- 500
for (alpha in c(1, 7)) {
  w <- as.vector(extraDistr::rdirichlet(1, rep(alpha, G)))
  counts <- rpois(G, mean_count * G * w)
  est <- alpha_from_cv(sd(counts) / mean(counts), G, mean(counts))
  check(near(est, alpha, 0.10), sprintf("alpha %g recovered as %.3f (within 10%%)", alpha, est))
}


# ---- 3. regime_rungs -------------------------------------------------------
cat("\n[3] regime_rungs builds a valid ray\n")
for (nm in names(SCALING_REGIMES)) {
  reg <- SCALING_REGIMES[[nm]]
  n_over_g <- if (nm == "gasperini") 207324 / 13077 else 616184 / 2666
  r <- regime_rungs(reg, n_over_g, nm)
  check(all(r$n_cells %% SCALING_PARAMS_SHARED$chunk_cells == 0),
        sprintf("%s: chunk_cells divides every n_cells", nm))
  check(all(abs(r$n_cells / r$n_guides / n_over_g - 1) < 0.01),
        sprintf("%s: n_cells/n_guides within 1%% of the target at every rung", nm))
  check(all(grepl(nm, r$dataset)),
        sprintf("%s: every dataset id contains the regime name", nm))
}


# ---- 4. simulate_regime passes every parameter -----------------------------
cat("\n[4] simulate_regime is the direct simulate_guidebender2 call\n")
tg <- list(count_per_cell = 100, hurdle_prob = 0.1, alpha = 2)
a <- simulate_regime(tg, 40L, 2000L, lambda = 3, snr = 5, seed = 9)
set.seed(9)
b <- simulate_guidebender2(
  n_guides = 40L, n_cells = 2000L, moi = 3, snr = 5, count_per_cell = 100,
  hurdle_prob = 0.1, guide_infection_alpha = 2, return_sparse_only = TRUE,
  Phi_cell = 1, Phi_noise = 0, frac_noise_endo = 0.75, rho_sum = 10, eps_alpha = 50,
  d_sigma_guide = 0.5, d_sigma_cell = 0.5, d_sigma_drop = 0.5,
  endo_shape_sum = 1, endo_shape_flat = 0, use_median = TRUE, chunk_cells = 1000)
check(identical(assay(a, "counts"), assay(b, "counts")),
      "identical counts to an explicit call with every parameter named")
check(identical(assay(a, "ground_truth"), assay(b, "ground_truth")),
      "identical ground truth")


# ---- 5. cal_stats agrees with guide_matrix_stats ---------------------------
# Calibration matches the simulation to the real matrix with cal_stats; the real
# matrix is measured with guide_matrix_stats. They must define MOI and nonzeros
# per cell identically.
cat("\n[5] cal_stats and guide_matrix_stats define the targets identically\n")
counts <- assay(a, "counts")
cs <- cal_stats(counts, 5L)
gs <- guide_matrix_stats(counts, 5L)
check(isTRUE(all.equal(cs[["moi"]], gs$moi)), "thresholded MOI agrees")
check(isTRUE(all.equal(cs[["nnz_per_cell"]], gs$nnz_per_cell)), "nonzeros per cell agrees")


# ---- 6. write_method_inputs ------------------------------------------------
cat("\n[6] write_method_inputs follows SCALING_METHODS\n")
colnames(counts) <- paste0("CELL_", seq_len(ncol(counts)))
rownames(counts) <- paste0("grna_", seq_len(nrow(counts)))
tmp <- file.path(tempdir(), "wmi"); unlink(tmp, recursive = TRUE)
dir.create(tmp)
write_method_inputs(counts, tmp)
check(file.exists(file.path(tmp, "cleanser", "grna_matrix.mtx")), "cleanser .mtx written")
check(all(file.exists(file.path(tmp, c("crispat", "pertpy"), "grna_matrix.h5ad"))),
      "crispat and pertpy .h5ad written")
check(!dir.exists(file.path(tmp, "fishash")),
      "fishash gets no directory of its own (main.nf maps it to cleanser/)")
fp_rds <- readRDS(file.path(tmp, "fishashplus", "grna_matrix.rds"))
check(identical(fp_rds, counts), "fishashplus .rds is the counts matrix, dimnames included")
check(!dir.exists(file.path(tmp, "sceptre")), "nothing is written to sceptre/")
unnamed <- counts; dimnames(unnamed) <- NULL
check(inherits(try(write_method_inputs(unnamed, tmp, "fishashplus"), silent = TRUE), "try-error"),
      "the fishashplus .rds refuses a matrix without names")
check(inherits(try(write_method_inputs(counts, tmp, "not_a_method"), silent = TRUE), "try-error"),
      "an unregistered method errors")
saved <- SCALING_METHODS
SCALING_METHODS$fake <- list(input = "zarr")
check(inherits(try(write_method_inputs(counts, tmp, "fake"), silent = TRUE), "try-error"),
      "a registered method with no writer for its input type errors")
SCALING_METHODS <- saved


# ---- 7. solve_log finds known roots ----------------------------------------
# The starting bracket is x0 * exp(+-0.2); every root here lies outside it, so
# uniroot has to extend the bracket in the stated direction.
cat("\n[7] solve_log finds the root of a known monotone function\n")
cases <- list(
  list(f = function(x) x^2, target = 50,   x0 = 1, dir = "upX",   root = sqrt(50), what = "increasing, root above the start"),
  list(f = function(x) x,   target = 0.02, x0 = 1, dir = "upX",   root = 0.02,     what = "increasing, root below the start"),
  list(f = function(x) 1/x, target = 0.01, x0 = 1, dir = "downX", root = 100,      what = "decreasing, root above the start")
)
for (cs in cases) {
  est <- solve_log(cs$f, cs$target, cs$x0, cs$dir, 1e-6, 1e6, "x")
  check(near(est, cs$root, 0.005), sprintf("%s: %.4g found for %.4g (within 0.5%%)", cs$what, est, cs$root))
}
# A target out of reach inside the bounds must stop with a message, from either side.
check(inherits(try(solve_log(function(x) x, 5, 0.5, "upX", 1e-3, 1, "x"), silent = TRUE), "try-error"),
      "target above f's range on [lower, upper] errors")
check(inherits(try(solve_log(function(x) 2 + x, 1, 0.5, "upX", 1e-3, 1, "x"), silent = TRUE), "try-error"),
      "target below f's range on [lower, upper] errors (the lambda-to-0 case)")


# ---- 8. calibrate recovers known parameters --------------------------------
# The targets are simulated at (lambda0, snr0) with the same seed and n_cal that
# calibrate's evaluations use, so the targets are exactly reachable at the true
# parameters: this tests the solver, not sampling error. The starting snr is
# far from the truth.
cat("\n[8] calibrate recovers lambda and snr from targets they produced\n")
G <- 100L; lambda0 <- 4; snr0 <- 20; seed <- 11L
reg <- list(threshold = 5L, n_cal = 4000L, n_seeds = 1L, cal_tol = 0.005)
tg  <- list(count_per_cell = 200, hurdle_prob = 0.1, alpha = 2, snr_start = 5)
st0 <- cal_stats(assay(simulate_regime(tg, G, reg$n_cal, lambda0, snr0, seed), "counts"),
                 reg$threshold)
tg$moi <- st0[["moi"]]; tg$nnz_per_cell <- st0[["nnz_per_cell"]]
rungs <- data.frame(n_guides = G, dataset = "test_rung")
cal <- calibrate(reg, tg, rungs, seed, "per_rung")
r <- cal[cal$role == "rung", ]
check(nrow(r) == 1, "one row for the one rung")
check(near(r$cal_moi, tg$moi, reg$cal_tol), sprintf("MOI %.4f matches target %.4f", r$cal_moi, tg$moi))
check(near(r$cal_nnz_per_cell, tg$nnz_per_cell, reg$cal_tol),
      sprintf("nonzeros per cell %.4f matches target %.4f", r$cal_nnz_per_cell, tg$nnz_per_cell))
check(near(r$lambda, lambda0, 0.05), sprintf("lambda %.4f recovered for %g (within 5%%)", r$lambda, lambda0))
check(near(r$snr, snr0, 0.10), sprintf("snr %.3f recovered for %g (within 10%%)", r$snr, snr0))


# ---- 9. the two snr modes across two rungs ---------------------------------
# Same targets as [8], now at two numbers of guides. "per_rung" must match both
# targets at every rung; "shared" must use one snr everywhere and match MOI at
# every rung (nonzeros per cell only at the rung where snr was solved).
cat("\n[9] snr modes across two rungs\n")
rungs2 <- data.frame(n_guides = c(60L, 150L), dataset = c("test_g60", "test_g150"))
pr <- calibrate(reg, tg, rungs2, seed, "per_rung")
check(identical(pr$role, c("rung", "rung")), "per_rung: one row per rung and no snr row")
for (i in 1:2) {
  check(near(pr$cal_moi[i], tg$moi, reg$cal_tol),
        sprintf("per_rung, G = %d: MOI %.4f matches target", pr$n_guides[i], pr$cal_moi[i]))
  check(near(pr$cal_nnz_per_cell[i], tg$nnz_per_cell, reg$cal_tol),
        sprintf("per_rung, G = %d: nonzeros per cell %.4f matches target", pr$n_guides[i],
                pr$cal_nnz_per_cell[i]))
}
cat(sprintf("  (per_rung snr: %.3f at G = 60, %.3f at G = 150)\n", pr$snr[1], pr$snr[2]))
sh <- calibrate(reg, tg, rungs2, seed, "shared", cal_n_guides = 150L)
shr <- sh[sh$role == "rung", ]
check(sh$role[1] == "snr" && nrow(shr) == 2, "shared: the snr row, then one row per rung")
check(length(unique(sh$snr)) == 1, "shared: one snr on every row")
check(all(near(shr$cal_moi, tg$moi, reg$cal_tol)), "shared: MOI matches target at every rung")
check(all(sh$snr_mode == "shared") && all(pr$snr_mode == "per_rung"), "snr_mode recorded on every row")
check(inherits(try(calibrate(reg, tg, rungs2, seed, "shared"), silent = TRUE), "try-error"),
      "shared without cal_n_guides errors")

cat("\nAll tests passed.\n")
