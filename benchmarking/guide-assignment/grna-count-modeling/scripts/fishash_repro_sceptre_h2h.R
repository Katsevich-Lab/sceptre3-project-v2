#!/usr/bin/env Rscript
# ============================================================================
# Head-to-head: current SCEPTRE mixture assignment  vs  Fishash+ Poisson.
# For the "should sceptre adopt fishash+" decision. Runs BOTH on the authors'
# identical guidebender sims on the SAME machine, SERIALLY (mc.cores=1) so the
# runtime comparison is fair. Accuracy = per-entry pooled F1/precision/recall
# (the paper's metric), full subset. Resumable (skips datasets already scored).
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(sceptre); library(fishash); library(Matrix); library(SummarizedExperiment) })
source("scripts/fishash_repro_lib.R"); source("scripts/sim_lib.R")

DATA <- sim_data_dir()          # stable dataset cache (lib helper)
OUT <- "results/fishash_repro/sceptre_h2h"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# dataset list: varyMOI (5 iters x 8 MOI, 200 guides) + varyNguides (5 iters x {20,200,2000} + 1 x 20000)
jobs <- rbind(
  do.call(rbind, lapply(c(.1,.3,.5,1,2,3,5,10), function(m) data.frame(
    scenario = "varyMOI", x = m, iter = 1:5,
    path = file.path(DATA, "varyMOI", sprintf("moi_%f_iter_%d.Rds", m, 1:5))))),
  do.call(rbind, lapply(c(20,200,2000), function(g) data.frame(
    scenario = "varyNguides", x = g, iter = 1:5,
    path = file.path(DATA, "varyNguides", sprintf("nguides_%d_iter_%d.Rds", g, 1:5))))),
  data.frame(scenario = "varyNguides", x = 20000, iter = 1,
             path = file.path(DATA, "varyNguides", "nguides_20000_iter_1.Rds"))
)
jobs <- jobs[file.exists(jobs$path), ]

# scoring via the lib's score_entries() (per-entry pooled, full subset); no duplicated inline scorer.
mk <- function(v, j, m, sec) data.frame(scenario = j$scenario, x = j$x, iter = j$iter,
                          method = m, F1 = v["F1"], Precision = v["Precision"], Recall = v["Recall"], secs = sec)
rows <- list()
for (k in seq_len(nrow(jobs))) {
  j <- jobs[k, ]; tag <- sprintf("%s_x%s_i%d", j$scenario, j$x, j$iter)
  done <- file.path(OUT, paste0(tag, ".csv"))
  if (file.exists(done)) { rows[[k]] <- read.csv(done); next }
  sim <- readRDS(j$path); counts <- assay(sim, "counts")
  ts <- system.time(As <- tryCatch(sceptre_assign(counts, moi = "high"), error = function(e) NULL))["elapsed"]
  tp <- system.time(Ap <- fishash_plus_poisson(counts, cell_margin = "ambient"))["elapsed"]
  vs <- if (!is.null(As)) score_entries(sim, As) else c(F1 = NA, Precision = NA, Recall = NA)
  vp <- score_entries(sim, Ap)
  out <- rbind(if (!is.null(As)) mk(vs, j, "sceptre_mixture", ts), mk(vp, j, "fishash+ Poisson", tp))
  write.csv(out, done, row.names = FALSE); rows[[k]] <- out
  cat(sprintf("%-26s %5.0f guides | sceptre %6.1fs (F1 %.3f) | fishash+ %5.2fs (F1 %.3f)\n",
      tag, nrow(counts), ts, vs["F1"], tp, vp["F1"]))
  rm(sim, counts); gc()
}
df <- do.call(rbind, rows); rownames(df) <- NULL
write.csv(df, file.path(OUT, "sceptre_h2h.csv"), row.names = FALSE)
cat("\nDONE ->", file.path(OUT, "sceptre_h2h.csv"), nrow(df), "rows\n")
