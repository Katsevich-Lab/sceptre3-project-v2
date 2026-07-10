#!/usr/bin/env Rscript
# ============================================================================
# How the endo/exo noise split modulates the assignment comparison, across the
# paper's TWO realistic regimes -- high-gRNA (25% exogenous) and low-gRNA (75%
# exogenous) -- both swept over MOI (varyMOI vs varyMOILow). No invented settings.
#
# Compares SCEPTRE mixture, fishash, and Fishash+ Poisson. fishash / fishash+ are
# read from the existing method-panel runs; sceptre is run here (resumable cache).
#
# Mechanism (paper Eq 6): exogenous noise scales with the cell's OWN library
# (d_cell+d_drop), so fishash's raw-library draws are well-specified for it while
# depth_fix's denoised depth under-counts it -> the depth advantage shrinks as the
# exo fraction rises. (The two regimes also differ in SNR/UMIs: the paper's realistic
# pair, not a clean single-variable isolation.)
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({
  library(sceptre); library(fishash); library(Matrix); library(SummarizedExperiment); library(parallel)
  library(ggplot2); library(dplyr); library(tidyr)
})
source("scripts/fishash_repro_lib.R"); source("scripts/sim_lib.R")
OUT <- "results/fishash_repro/endo_contrast"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
scen <- c(varyMOI = "high-gRNA (25% exo)", varyMOILow = "low-gRNA (75% exo)")
MOI  <- c(.1, .3, .5, 1, 2, 3, 5, 10); ITERS <- 1:10
cols <- c(sceptre_mixture = "#1f78b4", fishash = "#6a3d9a", `fishash+` = "#e31a1c")

# ---- run SCEPTRE on both regimes' cached datasets (parallel, resumable) -------
sc_path <- file.path(OUT, "sceptre_scores.csv")
if (file.exists(sc_path)) {
  sc <- read.csv(sc_path)
} else {
  jobs <- expand.grid(scn = names(scen), moi = MOI, iter = ITERS, stringsAsFactors = FALSE)
  rows <- mclapply(seq_len(nrow(jobs)), function(k) {
    j <- jobs[k, ]
    p <- file.path(sim_data_dir(j$scn), sprintf("moi_%f_iter_%d.Rds", j$moi, j$iter))
    if (!file.exists(p)) return(NULL)
    sim <- readRDS(p); A <- tryCatch(sceptre_assign(assay(sim, "counts"), moi = "high"), error = function(e) NULL)
    if (is.null(A)) return(NULL)
    v <- score_entries(sim, A)
    data.frame(scn = j$scn, moi = j$moi, method = "sceptre_mixture",
               Precision = v["Precision"], Recall = v["Recall"], F1 = v["F1"])
  }, mc.cores = as.integer(Sys.getenv("FISH_CORES", 8)))
  sc <- do.call(rbind, Filter(Negate(is.null), rows)); write.csv(sc, sc_path, row.names = FALSE)
  message(sprintf("sceptre: scored %d datasets", nrow(sc)))
}

# ---- fishash + fishash+ from the existing method-panel combined CSVs ----------
load_fp <- function(s) {
  d <- read.csv(file.path("results/fishash_repro", s, "combined_confusion.csv"))
  d <- d[d$subset == "full", ]; d$F1 <- with(d, 2 * Precision * Recall / (Precision + Recall))
  d$moi <- as.numeric(sub("moi_([0-9.]+)_.*", "\\1", d$sim_label))
  d %>% filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
    mutate(scn = s, method = recode(method, fishash_refit10 = "fishash", `fishash+_poisson` = "fishash+")) %>%
    select(scn, moi, method, Precision, Recall, F1)
}
d <- bind_rows(sc, load_fp("varyMOI"), load_fp("varyMOILow"))
d$regime <- factor(scen[d$scn], scen); d$method <- factor(d$method, names(cols))

ag <- d %>% group_by(regime, moi, method) %>%
  summarize(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1), .groups = "drop")
long <- ag %>% pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "v") %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))

p <- ggplot(long, aes(moi, v, color = method)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
  facet_grid(stat ~ regime) + scale_x_log10(breaks = c(.1,.3,.5,1,2,3,5,10)) +
  scale_color_manual(values = cols, labels = c("sceptre (mixture)", "fishash", "fishash+")) +
  labs(x = "MOI (200 guides, 20k cells)", y = "mean (per-entry, full)", color = NULL,
       title = "SCEPTRE vs fishash+ across the paper's endo-heavy vs exo-heavy regimes",
       subtitle = "fishash+ leads sceptre in both regimes; its edge over fishash (the depth fix) is larger under endo-heavy noise (left) than exo-heavy (right).") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT, "endo_exo_contrast.png"), p, width = 11, height = 7, dpi = 130)

gap <- ag %>% select(regime, moi, method, F1) %>% pivot_wider(names_from = method, values_from = F1)
cat("=== F1: sceptre vs fishash+ gap (fishash+ - sceptre) by MOI x regime ===\n")
print(as.data.frame(gap %>% mutate(`fishash+ - sceptre` = round(`fishash+` - sceptre_mixture, 3)) %>%
  select(regime, moi, `fishash+ - sceptre`) %>%
  pivot_wider(names_from = regime, values_from = `fishash+ - sceptre`)), row.names = FALSE)
cat("\nwrote", file.path(OUT, "endo_exo_contrast.png"), "\n")
