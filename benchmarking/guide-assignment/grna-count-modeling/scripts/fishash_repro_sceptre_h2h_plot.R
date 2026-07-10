#!/usr/bin/env Rscript
# Head-to-head plots: current SCEPTRE mixture vs Fishash+ Poisson (accuracy + runtime).
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
OUT <- "results/fishash_repro/sceptre_h2h"
d <- read.csv(file.path(OUT, "sceptre_h2h.csv"))
d$method <- factor(d$method, c("sceptre_mixture", "fishash+ Poisson"))
cols <- c(sceptre_mixture = "#1f78b4", `fishash+ Poisson` = "#e31a1c")

ag <- d %>% group_by(scenario, x, method) %>%
  summarize(F1 = mean(F1), Precision = mean(Precision), Recall = mean(Recall),
            secs = median(secs), .groups = "drop")
write.csv(ag, file.path(OUT, "sceptre_h2h_summary.csv"), row.names = FALSE)

# ---- ACCURACY: F1 vs MOI, and F1 vs nguides -------------------------------
acc <- ag %>% pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "v") %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))
pm <- acc %>% filter(scenario == "varyMOI") %>%
  ggplot(aes(x, v, color = method)) + geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~stat, nrow = 1) + scale_x_log10(breaks = c(.1,.3,.5,1,2,3,5,10)) +
  scale_color_manual(values = cols) + coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "MOI (200 guides)", y = "mean (full subset)", color = NULL,
       title = "Accuracy: SCEPTRE mixture vs Fishash+ Poisson — varying MOI") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "h2h_accuracy_moi.png"), pm, width = 11, height = 4.2, dpi = 130)

pn <- acc %>% filter(scenario == "varyNguides") %>%
  ggplot(aes(x, v, color = method)) + geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_wrap(~stat, nrow = 1) + scale_x_log10(breaks = c(20,200,2000,20000)) +
  scale_color_manual(values = cols) + coord_cartesian(ylim = c(0.6, 1)) +
  labs(x = "number of guides (MOI 0.3)", y = "mean (full subset)", color = NULL,
       title = "Accuracy: SCEPTRE mixture vs Fishash+ Poisson — varying guide count") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "h2h_accuracy_nguides.png"), pn, width = 11, height = 4.2, dpi = 130)

# ---- RUNTIME: seconds vs nguides (log-log) — the scaling story -------------
rt <- ag %>% filter(scenario == "varyNguides")
spd <- rt %>% select(x, method, secs) %>% pivot_wider(names_from = method, values_from = secs) %>%
  mutate(speedup = round(sceptre_mixture / `fishash+ Poisson`))
pr <- ggplot(rt, aes(x, secs, color = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.4) +
  geom_text(data = spd, aes(x = x, y = sceptre_mixture, label = paste0(speedup, "× faster")),
            inherit.aes = FALSE, vjust = -0.8, size = 3.6, color = "grey30") +
  scale_x_log10(breaks = c(20,200,2000,20000)) + scale_y_log10() +
  scale_color_manual(values = cols) +
  labs(x = "number of guides", y = "seconds per dataset (20k cells, 1 CPU)", color = NULL,
       title = "Runtime: SCEPTRE mixture vs Fishash+ Poisson",
       subtitle = "Same machine, single-threaded. SCEPTRE's per-guide mixture fitting scales with library size; Fishash+ is flat.") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 9))
ggsave(file.path(OUT, "h2h_runtime_nguides.png"), pr, width = 9, height = 6, dpi = 130)

cat("=== accuracy (mean F1) ===\n")
print(as.data.frame(ag %>% select(scenario, x, method, F1) %>%
        pivot_wider(names_from = method, values_from = F1)), row.names = FALSE)
cat("\n=== runtime (median s/dataset) + speedup, by guide count ===\n")
print(as.data.frame(spd), row.names = FALSE)
cat("\nwrote h2h_accuracy_moi.png, h2h_accuracy_nguides.png, h2h_runtime_nguides.png to", OUT, "\n")
