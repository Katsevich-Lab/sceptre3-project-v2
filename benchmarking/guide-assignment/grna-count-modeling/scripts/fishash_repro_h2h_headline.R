#!/usr/bin/env Rscript
# Side-by-side headline: SCEPTRE mixture vs Fishash+ Poisson — accuracy (F1) and runtime.
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
source("scripts/fishash_repro_axis.R")
OUT <- "results/fishash_repro/sceptre_h2h"
cols <- c(sceptre_mixture = "#1f78b4", `fishash+ Poisson` = "#e31a1c")
labs2 <- c("sceptre (mixture)", "fishash+")

# ---- (A) accuracy: mean F1 vs MOI (20 reps) ----
acc <- read.csv(file.path(OUT, "h2h_accuracy_reps.csv"))
acc$method <- factor(acc$method, names(cols))
agA <- acc %>% group_by(moi, method) %>% summarize(F1 = mean(F1), .groups = "drop")
pA <- ggplot(agA, aes(moi, F1, color = method)) +
  geom_line(linewidth = .9) + geom_point(size = 2) +
  scale_x_fishash_guide_load(secondary = FALSE) +
  scale_color_manual(values = cols, labels = labs2) +
  coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "Mean infection events per recovered cell", y = "F1 (per-entry, full)",
       color = NULL, title = "Accuracy") +
  theme_bw(base_size = 17)

# ---- (B) runtime: median seconds vs nguides ----
rt <- read.csv(file.path(OUT, "sceptre_h2h.csv")) %>% filter(scenario == "varyNguides")
rt$method <- factor(rt$method, names(cols))
agB <- rt %>% group_by(x, method) %>% summarize(secs = median(secs), .groups = "drop")
spd <- agB %>% pivot_wider(names_from = method, values_from = secs) %>%
  mutate(lab = paste0(round(sceptre_mixture / `fishash+ Poisson`), "×"))
pB <- ggplot(agB, aes(x, secs, color = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.4) +
  geom_text(data = spd, aes(x = x, y = sceptre_mixture, label = lab), inherit.aes = FALSE, vjust = -0.9, size = 5, color = "grey30") +
  scale_x_log10(breaks = c(20,200,2000,20000)) +
  scale_y_log10(expand = expansion(mult = c(0.04, 0.18))) +   # headroom for the 1030x label
  scale_color_manual(values = cols, labels = labs2) + guides(color = "none") +
  labs(x = "number of guides (20k cells)", y = "seconds / dataset (1 CPU)", color = NULL, title = "Runtime") +
  theme_bw(base_size = 17)

combined <- (pA | pB) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 18), legend.key.size = unit(1.1, "lines"))
ggsave(file.path(OUT, "h2h_headline_sidebyside.png"), combined, width = 13.5, height = 5.6, dpi = 135)
cat("wrote", file.path(OUT, "h2h_headline_sidebyside.png"), "\n")
