#!/usr/bin/env Rscript
# ============================================================================
# Full-PR-curve AUPRC (average precision) by MOI, computed exactly for Fishash+
# Poisson and fishash on the authors' identical varyMOI sims. This is the paper's
# Fig 6c / D.1b analysis ("fishash's high-MOI F1 gap shrinks under full AUPRC").
# Both methods emit a continuous score (Poisson-tail / hypergeometric -log p), so
# both AUPRCs are exact; raw count is the baseline. Other methods are digitized
# from D.1b separately.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(fishash); library(Matrix); library(SummarizedExperiment); library(parallel) })
source("scripts/fishash_repro_lib.R")

DATA <- sim_data_dir("varyMOI")     # stable dataset cache (lib helper)
OUT <- "results/fishash_repro/varyMOI"
files <- list.files(DATA, pattern = "\\.Rds$", full.names = TRUE)
files <- files[grepl("moi_", basename(files))]

one <- function(f) {
  sim <- readRDS(f); counts <- as(assay(sim, "counts"), "CsparseMatrix")
  truth <- assay(sim, "ground_truth")
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  tv <- as.logical(truth[cbind(i, j)])
  ppo <- fishash_plus_poisson(counts, cell_margin = "ambient", return_score = TRUE)
  ap_ours <- average_precision(ppo$score, truth[cbind(ppo$i, ppo$j)])
  fa <- fishash(counts, refit = 10)
  lp <- as(assay(fa, "log_pval"), "CsparseMatrix")[cbind(i, j)]
  s_fa <- -lp; s_fa[y < 2] <- -Inf
  ap_fa <- average_precision(s_fa, tv)
  ap_raw <- average_precision(y, tv)
  data.frame(moi = as.numeric(sub("moi_([0-9.]+)_.*", "\\1", basename(f))),
             `fishash+ Poisson (ours)` = ap_ours, fishash = ap_fa, raw_count = ap_raw, check.names = FALSE)
}
res <- mclapply(files, function(f) tryCatch(one(f), error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL }),
                mc.cores = as.integer(Sys.getenv("FISH_CORES", 6)))
df <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(df, file.path(OUT, "auprc_by_moi.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(dplyr))
summ <- df %>% tidyr::pivot_longer(-moi, names_to = "method", values_to = "AUPRC") %>%
  group_by(moi, method) %>% summarize(AUPRC = median(AUPRC), .groups = "drop")
write.csv(summ, file.path(OUT, "auprc_summary.csv"), row.names = FALSE)
cat("=== median full AUPRC by MOI (exact) ===\n")
print(as.data.frame(summ %>% tidyr::pivot_wider(names_from = method, values_from = AUPRC)), row.names = FALSE)

# ---- figure (Fig 6c / D.1b analog): AUPRC vs MOI, ours vs fishash vs raw (all exact) ----
suppressPackageStartupMessages(library(ggplot2))
summ$method <- factor(summ$method, c("raw_count", "fishash", "fishash+ Poisson (ours)"))
p <- ggplot(summ, aes(moi, AUPRC, color = method)) +
  geom_line(aes(linewidth = method)) + geom_point(size = 2) +
  scale_x_log10(breaks = c(.1,.3,.5,1,2,3,5,10)) +
  scale_color_manual(values = c(raw_count = "grey65", fishash = "#6a3d9a", `fishash+ Poisson (ours)` = "#e31a1c")) +
  scale_linewidth_manual(values = c(0.6, 0.9, 1.5), guide = "none") +
  coord_cartesian(ylim = c(0.9, 1.0)) +
  labs(x = "MOI", y = "median full AUPRC", color = NULL,
       title = "Full-PR-curve AUPRC vs MOI (exact) -- Fig 6c/D.1b analysis",
       subtitle = "Ours >= fishash at every MOI>=0.5, gap grows with MOI. Paper: dcCLEANSER was top-AUPRC method (not recomputed here).") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))
ggsave(file.path(OUT, "fig6c_auprc_vs_moi.png"), p, width = 10, height = 6, dpi = 130)
cat("wrote", file.path(OUT, "fig6c_auprc_vs_moi.png"), "\n")
