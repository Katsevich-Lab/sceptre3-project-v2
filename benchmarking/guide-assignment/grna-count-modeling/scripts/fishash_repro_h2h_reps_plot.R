#!/usr/bin/env Rscript
# Error-bar accuracy plot for the higher-rep sceptre-vs-fishash+ head-to-head.
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
OUT <- "results/fishash_repro/sceptre_h2h"
d <- read.csv(file.path(OUT, "h2h_accuracy_reps.csv"))
d$method <- factor(d$method, c("sceptre_mixture", "fishash+ Poisson"))
cols <- c(sceptre_mixture = "#1f78b4", `fishash+ Poisson` = "#e31a1c")
nrep <- max(table(d$moi, d$method))

long <- d %>% pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "v") %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))
ag <- long %>% group_by(moi, method, stat) %>%
  summarize(mean = mean(v, na.rm = TRUE), se = sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))),
            sd = sd(v, na.rm = TRUE), n = sum(!is.na(v)), .groups = "drop")
write.csv(ag, file.path(OUT, "h2h_accuracy_reps_summary.csv"), row.names = FALSE)

p <- ggplot() +
  geom_point(data = long, aes(moi, v, color = method), position = position_jitter(width = .04, height = 0),
             alpha = .18, size = .8) +
  geom_line(data = ag, aes(moi, mean, color = method), linewidth = .8) +
  geom_errorbar(data = ag, aes(moi, ymin = mean - se, ymax = mean + se, color = method), width = .06, linewidth = .7) +
  geom_point(data = ag, aes(moi, mean, color = method), size = 1.9) +
  facet_wrap(~stat, nrow = 1) + scale_x_log10(breaks = c(.1,.3,.5,1,2,3,5,10)) +
  scale_color_manual(values = cols) + coord_cartesian(ylim = c(0.55, 1)) +
  labs(x = "MOI (200 guides, 20k cells)", y = "per-entry (full subset)", color = NULL,
       title = sprintf("SCEPTRE mixture vs Fishash+ Poisson — accuracy over %d replicates", nrep),
       subtitle = "Bars = mean +/- 1 SE; faint points = individual replicates. The MOI 0.3 sceptre-precision wobble regresses with more reps.") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8.5))
ggsave(file.path(OUT, "h2h_accuracy_moi_errorbars.png"), p, width = 12, height = 4.6, dpi = 130)

cat(sprintf("=== precision by MOI (%d reps): mean +/- SE ===\n", nrep))
ptab <- ag %>% filter(stat == "Precision") %>%
  transmute(moi, method, precision = sprintf("%.3f +/- %.3f", mean, se)) %>%
  pivot_wider(names_from = method, values_from = precision)
print(as.data.frame(ptab), row.names = FALSE)
cat("\nwrote h2h_accuracy_moi_errorbars.png\n")
