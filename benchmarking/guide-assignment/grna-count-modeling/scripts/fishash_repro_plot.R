#!/usr/bin/env Rscript
# ============================================================================
# Plot Fishash+ (Poisson) vs fishash on the authors' own simulations.
# Reproduces the authors' precision/recall/F1 figures (bin/combine_confusion.R,
# 61_plot_simulation_results.R) with our methods overlaid.
#   Rscript scripts/fishash_repro_plot.R <scenario>
# ============================================================================
setwd("/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
source("scripts/fishash_repro_axis.R")
SCEN <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "varyMOI"
OUT  <- file.path("results", "fishash_repro", SCEN)

d <- read.csv(file.path(OUT, "combined_confusion.csv"))
d$F1 <- with(d, 2 * Precision * Recall / (Precision + Recall))

# Method display order / labels. The main-text figure shows only the two
# decision-relevant methods; the full five-method decomposition is retained as
# an appendix figure.
lev <- c("fishash_refit10", "fishash+_poisson", "fishash+_contingency",
         "fishash_poisson_rawd", "fishash_refit0")
lab <- c(fishash_refit10 = "fishash (refit10)",
         `fishash+_poisson` = "fishash+ Poisson (ours)",
         `fishash+_contingency` = "fishash+ contingency (depth_fix)",
         fishash_poisson_rawd = "fishash Poisson, raw depth",
         fishash_refit0 = "fishash (refit0, uncorrected)")
cols <- c("fishash (refit10)" = "#1f77b4", "fishash+ Poisson (ours)" = "#d62728",
          "fishash+ contingency (depth_fix)" = "#ff7f0e",
          "fishash Poisson, raw depth" = "#7f7f7f",
          "fishash (refit0, uncorrected)" = "#bcbd22")
headline <- c("fishash (refit10)", "fishash+ Poisson (ours)")
d <- d %>% filter(method %in% lev) %>%
  mutate(method = factor(lab[as.character(method)], levels = lab[lev]))

# x-axis variable per scenario
if (SCEN == "varyMOI") {
  d$x <- as.numeric(sub("moi_([0-9.]+)_.*", "\\1", d$sim_label))
  d$facet_x <- factor(
    fishash_recovered_guide_labels(d$x),
    levels = fishash_recovered_guide_labels(fishash_preselection_moi)
  )
  xlab <- "Mean infection events per recovered cell"; xlog <- TRUE
} else if (SCEN == "varyNguides") {
  d$x <- as.numeric(sub("nguides_([0-9]+)_.*", "\\1", d$sim_label)); xlab <- "n guides"; xlog <- TRUE
} else {
  d$x <- factor(sub("corr_([^_]+)_.*", "\\1", d$sim_label), levels = c("low", "mid", "high"))
  xlab <- "signal-noise corr"; xlog <- FALSE
}

# ---- summary table (full subset): mean over replicates -----------------------
summ <- d %>% filter(subset == "full") %>%
  group_by(x, method) %>%
  summarize(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1),
            F1_sd = sd(F1), n = n(), .groups = "drop")
write.csv(summ, file.path(OUT, "summary_by_x.csv"), row.names = FALSE)

tab <- summ %>% select(x, method, F1) %>%
  pivot_wider(names_from = method, values_from = F1) %>% arrange(x)
cat("=== mean F1 by", xlab, "(full subset) ===\n"); print(as.data.frame(tab), row.names = FALSE)

# recall gain of ours over fishash
gain <- summ %>% filter(method %in% lab[c("fishash_refit10","fishash+_poisson")]) %>%
  select(x, method, Recall, Precision, F1) %>%
  pivot_wider(names_from = method, values_from = c(Recall, Precision, F1))
cat("\n=== ours vs fishash: recall / precision / F1 by", xlab, "===\n")
print(as.data.frame(gain), row.names = FALSE)

# ---- Figure 1: headline comparison, P / R / F1 vs x --------------------------
long <- d %>% filter(subset == "full") %>%
  pivot_longer(c(Precision, Recall, F1), names_to = "stat", values_to = "value") %>%
  mutate(stat = factor(stat, c("Precision", "Recall", "F1")))
ag <- long %>% group_by(x, method, stat) %>%
  summarize(m = mean(value), se = sd(value)/sqrt(n()), .groups = "drop")

p1 <- ggplot(ag %>% filter(method %in% headline),
             aes(x, m, color = method, group = method)) +
  { if (SCEN == "varyMOI") scale_x_fishash_guide_load()
    else if (xlog) scale_x_log10() } +
  geom_line(linewidth = 1.05) +
  geom_pointrange(aes(ymin = m - se, ymax = m + se), linewidth = 0.35, size = 0.45) +
  facet_wrap(~stat, nrow = 1) +
  scale_color_manual(values = cols, breaks = headline,
                     labels = c("fishash (refit10)", "Fishash+ Poisson")) +
  labs(x = xlab, y = "mean over replicates (full subset)", color = NULL,
       title = "Fishash+ Poisson recovers recall as guide load rises",
       subtitle = sprintf("Exact comparison with fishash on the authors' guidebender simulations (%s)", SCEN)) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
ggsave(file.path(OUT, paste0(SCEN, "_prf_vs_x.png")), p1, width = 12, height = 5, dpi = 130)

# ---- Appendix figure: full five-method ablation ------------------------------
p1a <- ggplot(ag, aes(x, m, color = method, group = method)) +
  { if (SCEN == "varyMOI") scale_x_fishash_guide_load()
    else if (xlog) scale_x_log10() } +
  geom_line(linewidth = 0.7) +
  geom_pointrange(aes(ymin = m - se, ymax = m + se), linewidth = 0.25, size = 0.3) +
  facet_wrap(~stat, nrow = 1) +
  scale_color_manual(values = cols) +
  labs(x = xlab, y = "mean over replicates (full subset)", color = NULL,
       title = sprintf("Diagnostic method decomposition on guidebender simulations (%s)", SCEN)) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))
ggsave(file.path(OUT, paste0(SCEN, "_ablation_prf_vs_x.png")), p1a,
       width = 12, height = 5, dpi = 130)

# ---- Figure 2: precision-recall scatter faceted by x (authors' style) --------
p2 <- d %>% filter(subset == "full") %>%
  ggplot(aes(Recall, Precision, color = method)) +
  geom_hline(yintercept = .95, lty = "dotted") +
  geom_point(shape = 1, alpha = .8) +
  { if (SCEN == "varyMOI") facet_wrap(~facet_x, nrow = 2,
                                      labeller = labeller(facet_x = function(z) paste(z, "infections/cell")))
    else facet_wrap(~x, labeller = label_both) } +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)) +
  labs(title = sprintf("Precision vs recall by %s (%s)", xlab, SCEN))
ggsave(file.path(OUT, paste0(SCEN, "_precrecall.png")), p2, width = 12, height = 7, dpi = 130)

cat("\nwrote figures to", OUT, "\n")
