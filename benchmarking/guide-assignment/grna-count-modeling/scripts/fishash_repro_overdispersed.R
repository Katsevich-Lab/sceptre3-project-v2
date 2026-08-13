#!/usr/bin/env Rscript
# ============================================================================
# Fig C.2 (paper Appendix C): overdispersed (Geometric) ambient noise, varying MOI.
# The robustness boundary -- the Poisson tail can't absorb non-Poisson ambient
# variance, so Fishash+ Poisson OVER-CALLS at low MOI (precision below nominal),
# while still winning F1 at high MOI where recall dominates.
#
# Generate the overdispersed data first (Phi_noise=1):
#   FISH_PHI_NOISE=1 Rscript scripts/fishash_repro_run.R varyMOI      # -> results/fishash_repro/varyMOI_od/
# Then this plots ours vs fishash (both exact) from that run.
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("scripts/fishash_repro_axis.R")
OUT <- "results/fishash_repro/varyMOI_od"

d <- read.csv(file.path(OUT, "combined_confusion.csv"))
d <- d[d$subset == "full", ]; d$F1 <- with(d, 2 * Precision * Recall / (Precision + Recall))
d$moi <- as.numeric(sub("moi_([0-9.]+)_.*", "\\1", d$sim_label))
s <- d %>% filter(method %in% c("fishash_refit10", "fishash+_poisson")) %>%
  mutate(method = recode(method, fishash_refit10 = "fishash", `fishash+_poisson` = "fishash+ Poisson (ours)")) %>%
  group_by(moi, method) %>% summarize(F1 = mean(F1), Precision = mean(Precision), Recall = mean(Recall), .groups = "drop") %>%
  pivot_longer(c(F1, Precision, Recall), names_to = "stat", values_to = "v") %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))

p <- ggplot(s, aes(moi, v, color = method)) +
  geom_hline(data = data.frame(stat = factor("Precision", c("Precision", "Recall", "F1"))),
             aes(yintercept = 0.95), lty = "dotted") +
  geom_line(linewidth = 1) + geom_point(size = 1.8) + facet_wrap(~stat, nrow = 1) +
  scale_x_fishash_guide_load() +
  scale_color_manual(values = c("fishash" = "#6a3d9a", "fishash+ Poisson (ours)" = "#e31a1c")) +
  labs(x = "Mean infection events per recovered cell", y = "mean (full subset)", color = NULL,
       title = "Overdispersed (Geometric) noise -- paper Fig C.2 setting: the robustness boundary",
       subtitle = "Poisson tail can't absorb non-Poisson ambient variance: ours OVER-CALLS at low MOI (precision 0.76 vs nominal 0.95); still wins F1 at high MOI.") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
ggsave(file.path(OUT, "figC2_overdispersed_prf.png"), p, width = 12, height = 4.8, dpi = 130)
cat("wrote", file.path(OUT, "figC2_overdispersed_prf.png"), "\n")
