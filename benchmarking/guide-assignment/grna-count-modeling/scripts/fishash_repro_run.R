#!/usr/bin/env Rscript
# ============================================================================
# Run Fishash+ (Poisson tail) through the fishash authors' OWN simulations.
#
# Faithful reproduction: we regenerate the authors' EXACT simulated datasets by
# replaying their generation scripts verbatim (same seed, same expand.grid order,
# same simulate_guidebender2 params -> byte-identical draws), then run the method
# panel and score with the authors' own per-entry confusion metric.
#
# Scenarios (from fishash_analysis/simulate/*.R):
#   varyMOI        : simulate_numCells20k_numGuides200_varyMOI.R
#   varyNguides    : ..._medUmi100_snr4_endo75_varyNumGuides.R  (Fig 5 high-gRNA, leq20k tiers)
#   varyNguidesLow : ..._medUmi20_snr1_endo25_varyNumGuides.R   (Fig 5 low-gRNA)
#   varyCorr       : ..._varySignalNoiseCorr.R                  (Fig 7 Simpson)
# Set FISH_PHI_NOISE=1 for the overdispersed (Geometric-noise) Appendix-C variants (-> *_od dirs).
#
# Usage:
#   Rscript scripts/fishash_repro_run.R <scenario> [n_iter] [n_cells]
#   n_iter  : FAITHFUL subset -- still replays the full RNG stream, scores only iters<=n_iter.
#   n_cells : PREVIEW -- reduced cells cannot match the authors' stream, so per-row seeds are used
#             (NOT the authors' datasets). Default (no args) = the authors' full, exact-RNG datasets.
# Data cache location: sim_data_dir() (results/fishash_repro/datasets, or $FISH_DATA).
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment); library(parallel)
})
source("scripts/fishash_repro_lib.R")
source("scripts/contingency_method.R")

args     <- commandArgs(trailingOnly = TRUE)
SCENARIO <- if (length(args) >= 1) args[1] else "varyMOI"
ITER_OV  <- if (length(args) >= 2) as.integer(args[2]) else NA
NCELL_OV <- if (length(args) >= 3) as.integer(args[3]) else NA
NCORES   <- as.integer(Sys.getenv("FISH_CORES", parallel::detectCores() - 1))
PHI      <- as.numeric(Sys.getenv("FISH_PHI_NOISE", "0"))   # >0 => overdispersed (Geometric) noise, paper Appendix C
NG_ONLY  <- Sys.getenv("FISH_NG_ONLY", "")                  # e.g. "80000" to run only that guide tier
SUFFIX   <- if (PHI > 0) "_od" else ""

DATA_DIR <- sim_data_dir(paste0(SCENARIO, SUFFIX))          # single stable path (lib helper; $FISH_DATA override)
OUT_DIR  <- file.path("results", "fishash_repro", paste0(SCENARIO, SUFFIX))
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,  recursive = TRUE, showWarnings = FALSE)
.chunk <- function(nc) if (nc %% 1000 == 0) 1000L else NULL  # robust chunking for non-1000-multiple preview sizes

# ---- scenario configs (verbatim from the authors' simulate_*.R) --------------
cfg <- switch(SCENARIO,
  varyMOI = list(
    seed = 202601060721 %% (2^31 - 1), ntimes = 10,
    grid = expand.grid(iter = 1:10, moi = c(.1, .3, .5, 1, 2, 3, 5, 10)),
    label = function(g) sprintf("moi_%f_iter_%d", g$moi, g$iter),
    sim = function(g, n_cells) simulate_guidebender2(
      n_guides = 200, n_cells = n_cells, moi = g$moi, hurdle_prob = .1, snr = 4,
      count_per_cell = 100, frac_noise_endo = .75, Phi_noise = PHI, return_sparse_only = TRUE, chunk_cells = .chunk(n_cells))),
  varyMOILow = list(   # DERIVED (not an authors' scenario): the paper's LOW-gRNA params (20 UMI, SNR 1,
                       # 75% exo) swept over the varyMOI grid -- pairs with varyMOI (25% exo) to show how the
                       # endo/exo split modulates the fishash-vs-fishash+ (depth) comparison across MOI.
    seed = 6060842, ntimes = 10,
    grid = expand.grid(iter = 1:10, moi = c(.1, .3, .5, 1, 2, 3, 5, 10)),
    label = function(g) sprintf("moi_%f_iter_%d", g$moi, g$iter),
    sim = function(g, n_cells) simulate_guidebender2(
      n_guides = 200, n_cells = n_cells, moi = g$moi, hurdle_prob = .1, snr = 1,
      count_per_cell = 20, frac_noise_endo = .25, Phi_noise = PHI, return_sparse_only = TRUE, chunk_cells = .chunk(n_cells))),
  varyNguides = list(          # Fig 5 HIGH-gRNA regime (100 UMI, SNR 4, 25% exo)
    seed = 511082, ntimes = 10,
    grid = expand.grid(iter = 1:10, nguides = c(20, 200, 2000, 20000)),  # 80k = authors' runtime/mem tier, held
    label = function(g) sprintf("nguides_%d_iter_%d", g$nguides, g$iter),
    sim = function(g, n_cells) simulate_guidebender2(
      n_guides = g$nguides, n_cells = n_cells, moi = .3, hurdle_prob = .1, snr = 4,
      count_per_cell = 100, frac_noise_endo = .75, Phi_noise = PHI, return_sparse_only = TRUE, chunk_cells = .chunk(n_cells))),
  varyNguidesLow = list(       # Fig 5 LOW-gRNA regime (20 UMI, SNR 1, 75% exo)
    seed = 37152107, ntimes = 10,
    grid = expand.grid(iter = 1:10, nguides = c(20, 200, 2000, 20000)),
    label = function(g) sprintf("nguides_%d_iter_%d", g$nguides, g$iter),
    sim = function(g, n_cells) simulate_guidebender2(
      n_guides = g$nguides, n_cells = n_cells, moi = .3, hurdle_prob = .1, snr = 1,
      count_per_cell = 20, frac_noise_endo = .25, Phi_noise = PHI, return_sparse_only = TRUE, chunk_cells = .chunk(n_cells))),
  varyCorr = list(
    seed = 202601050822 %% (2^31 - 1), ntimes = 20,
    grid = local({
      m <- data.frame(endo_exo_corr = c("high","mid","low"),
                      endo_shape_uniform = c(0,0,1), endo_shape_sum = c(1e6,1,1))
      # corr-block-major order (matches the authors' dplyr::cross_join): high x iter1..20, then mid, then low.
      # (base merge() reorders to iter-major, which misaligns the shared RNG stream -- audit HIGH finding.)
      do.call(rbind, lapply(seq_len(nrow(m)), function(r)
        data.frame(m[r, , drop = FALSE], iter = 1:20, row.names = NULL)))
    }),
    label = function(g) sprintf("corr_%s_unif_%f_sum_%g_iter_%d",
                                g$endo_exo_corr, g$endo_shape_uniform, g$endo_shape_sum, g$iter),
    sim = function(g, n_cells) simulate_guidebender2(
      n_guides = 20, n_cells = n_cells, moi = .3, hurdle_prob = .1, snr = 4, count_per_cell = 100,
      frac_noise_endo = 1, endo_shape_flat = g$endo_shape_uniform, endo_shape_sum = g$endo_shape_sum,
      return_sparse_only = TRUE, chunk_cells = .chunk(n_cells))),
  stop("unknown scenario"))

N_CELLS <- if (!is.na(NCELL_OV)) NCELL_OV else 20000L

# ---- Phase 1: regenerate the authors' exact datasets (shared RNG stream) ------
# For varyMOI/varyNguides the grid's expand.grid puts `iter` first (varies fastest),
# exactly as the authors' scripts; we replay set.seed once and loop in that order.
grid <- cfg$grid
if (!is.na(ITER_OV)) grid <- grid[grid$iter <= ITER_OV, , drop = FALSE]
if (nzchar(NG_ONLY) && "nguides" %in% names(grid)) grid <- grid[grid$nguides == as.integer(NG_ONLY), , drop = FALSE]

manifest_path <- file.path(DATA_DIR, "manifest.csv")
# PREVIEW (per-row seed, NOT the authors' datasets) is triggered ONLY by a reduced n_cells --
# that genuinely can't match the authors' stream. An iter-subset (ITER_OV alone) stays FAITHFUL:
# the exact-replay branch below still draws the full grid to keep the stream aligned, and Phase 2
# just scores the kept subset. (Previously ITER_OV also forced preview -> silently non-authors data.)
PREVIEW <- !is.na(NCELL_OV)
message(sprintf("[%s] generating %d datasets (n_cells=%d, %s) ...", SCENARIO, nrow(grid), N_CELLS,
                if (PREVIEW) "PREVIEW/per-row seed" else "exact RNG replay"))
if (PREVIEW) {
  for (i in seq_len(nrow(grid))) {                 # fast: only the kept subset, per-row seed
    g <- grid[i, , drop = FALSE]; lab <- cfg$label(g)
    p <- file.path(DATA_DIR, paste0(lab, ".Rds"))
    if (!file.exists(p)) { set.seed(cfg$seed + i); saveRDS(cfg$sim(g, N_CELLS), p) }
  }
} else {
  set.seed(cfg$seed)                               # replay FULL grid order to match authors' RNG stream
  for (i in seq_len(nrow(cfg$grid))) {
    g <- cfg$grid[i, , drop = FALSE]; lab <- cfg$label(g)
    p <- file.path(DATA_DIR, paste0(lab, ".Rds"))
    sim <- cfg$sim(g, N_CELLS)                      # ALWAYS draw (keep stream aligned)
    if (!file.exists(p)) saveRDS(sim, p)
    rm(sim)
  }
}
grid$sim_label <- vapply(seq_len(nrow(grid)), function(i) cfg$label(grid[i, , drop = FALSE]), "")
grid$path <- file.path(DATA_DIR, paste0(grid$sim_label, ".Rds"))
write.csv(grid, manifest_path, row.names = FALSE)

# ---- Phase 2: run method panel + score (parallel over datasets) --------------
run_one <- function(k) {
  lab <- grid$sim_label[k]
  done <- file.path(OUT_DIR, paste0(lab, "__confusion.csv"))
  if (file.exists(done)) return(read.csv(done))
  sim <- readRDS(grid$path[k]); counts <- assay(sim, "counts")
  methods <- list(
    fishash_refit0        = function() assay(fishash(counts, refit = 0),  "assigned"),
    fishash_refit10       = function() assay(fishash(counts, refit = 10), "assigned"),
    `fishash+_poisson`    = function() fishash_plus_poisson(counts, cell_margin = "ambient"),
    fishash_poisson_rawd  = function() fishash_plus_poisson(counts, cell_margin = "observed"),
    `fishash+_contingency`= function() contingency_assign(counts, cell_margin = "ambient",
                                                          tail = "hyper", fdr = "GS")$assigned
  )
  rows <- lapply(names(methods), function(m) {
    t0 <- Sys.time(); a <- methods[[m]]()
    s <- score_subsets(sim, a, m, lab); s$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    s
  })
  out <- do.call(rbind, rows); write.csv(out, done, row.names = FALSE); out
}

message(sprintf("[%s] running panel on %d datasets, %d cores ...", SCENARIO, nrow(grid), NCORES))
res <- mclapply(seq_len(nrow(grid)), function(k) tryCatch(run_one(k),
          error = function(e) { message(sprintf("  FAIL %s: %s", grid$sim_label[k], conditionMessage(e))); NULL }),
        mc.cores = NCORES, mc.preschedule = FALSE)
combined <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(combined, file.path(OUT_DIR, "combined_confusion.csv"), row.names = FALSE)
message(sprintf("[%s] DONE -> %s (%d rows)", SCENARIO, file.path(OUT_DIR, "combined_confusion.csv"), nrow(combined)))
