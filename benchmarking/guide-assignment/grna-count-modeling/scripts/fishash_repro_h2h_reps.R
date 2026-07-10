#!/usr/bin/env Rscript
# ============================================================================
# Higher-replicate ACCURACY-ONLY head-to-head (sceptre mixture vs Fishash+ Poisson)
# across MOI, to pin down the low-MOI precision variance (the "MOI 0.3 dip").
# Accuracy is independent of scheduling, so this is PARALLELIZED. iters 1-10 reuse
# the authors' seed-matched cached datasets; iters 11-N are fresh independent draws
# (per-(moi,iter) seed) purely to shrink the error bars. NOT for timing.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(sceptre); library(fishash); library(Matrix); library(SummarizedExperiment); library(parallel) })
source("scripts/fishash_repro_lib.R"); source("scripts/sim_lib.R")

CACHE <- sim_data_dir("varyMOI")   # authors' iters 1-10 (stable cache); iters 11-N drawn fresh below
OUT <- "results/fishash_repro/sceptre_h2h"; PER <- file.path(OUT, "reps_cache"); dir.create(PER, showWarnings = FALSE, recursive = TRUE)
N_REPS <- as.integer(Sys.getenv("FISH_REPS", "20"))
NCORES <- as.integer(Sys.getenv("FISH_CORES", "8"))
MOI <- c(.1, .3, .5, 1, 2, 3, 5, 10)

get_sim <- function(moi, iter) {
  if (iter <= 10) {
    p <- file.path(CACHE, sprintf("moi_%f_iter_%d.Rds", moi, iter))
    if (file.exists(p)) return(readRDS(p))
  }
  set.seed(700000L + as.integer(iter) * 1000L + as.integer(round(moi * 100)))  # fresh independent draw
  simulate_guidebender2(n_guides = 200, n_cells = 20000, moi = moi, hurdle_prob = .1, snr = 4,
                        count_per_cell = 100, frac_noise_endo = .75, return_sparse_only = TRUE, chunk_cells = 1000)
}
jobs <- expand.grid(moi = MOI, iter = seq_len(N_REPS))
one <- function(k) {
  moi <- jobs$moi[k]; iter <- jobs$iter[k]
  f <- file.path(PER, sprintf("moi%s_iter%d.csv", moi, iter))
  if (file.exists(f)) return(read.csv(f))
  sim <- get_sim(moi, iter); counts <- assay(sim, "counts")
  As <- tryCatch(sceptre_assign(counts, moi = "high"), error = function(e) NULL)
  Ap <- fishash_plus_poisson(counts, cell_margin = "ambient")
  mk <- function(A, m) { v <- score_entries(sim, A); data.frame(moi = moi, iter = iter, method = m,
                          F1 = v["F1"], Precision = v["Precision"], Recall = v["Recall"]) }
  out <- rbind(if (!is.null(As)) mk(As, "sceptre_mixture"), mk(Ap, "fishash+ Poisson"))
  write.csv(out, f, row.names = FALSE); out
}
message(sprintf("running %d MOI x %d reps = %d datasets on %d cores ...", length(MOI), N_REPS, nrow(jobs), NCORES))
res <- mclapply(seq_len(nrow(jobs)), function(k) tryCatch(one(k), error = function(e) { message("FAIL ", k, ": ", conditionMessage(e)); NULL }),
                mc.cores = NCORES, mc.preschedule = FALSE)
df <- do.call(rbind, Filter(Negate(is.null), res)); rownames(df) <- NULL
write.csv(df, file.path(OUT, "h2h_accuracy_reps.csv"), row.names = FALSE)
message(sprintf("DONE -> %s (%d rows, %d reps)", file.path(OUT, "h2h_accuracy_reps.csv"), nrow(df), N_REPS))
