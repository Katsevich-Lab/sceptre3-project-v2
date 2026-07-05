#!/usr/bin/env Rscript
# Cross-dataset summary of ambient DEPTH + COMPOSITION (reads
# results/ambient_intuition/survey_depth_composition.csv). One 4-panel figure.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
D <- file.path(HERE, "results", "ambient_intuition")
t <- read.csv(file.path(D, "survey_depth_composition.csv"))

t$lab <- t$dataset |>
  (\(d) sub("_GSE[0-9]+$", "", d))() |>
  (\(d) sub("_figshare.*$", "", d))() |>
  (\(d) gsub("_", " ", d))()
t$lab <- substr(t$lab, 1, 24)
t$ambient_ok <- t$n_well >= 5                     # enough well-observed guides for composition
th <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8), legend.position = "none",
        axis.text.y = element_text(size = 7.5))

# A. depth contamination: median + p90 fold by which library size overstates ambient depth
tA <- t[order(t$depth_ratio_med), ]; tA$lab <- factor(tA$lab, levels = tA$lab)
pA <- ggplot(tA, aes(y = lab)) +
  geom_segment(aes(x = depth_ratio_med, xend = depth_ratio_p90, yend = lab),
               colour = "grey70", linewidth = 1) +
  geom_point(aes(x = depth_ratio_med, colour = log10(moi)), size = 2.4) +
  geom_point(aes(x = depth_ratio_p90), colour = "grey55", size = 1.1) +
  scale_x_log10() + scale_colour_viridis_c(option = "C", end = 0.9) +
  labs(title = "A. Library size overstates ambient depth — everywhere",
       subtitle = "dot = median fold (libsize / ambient depth); grey tick = 90th-percentile cell; colour = MOI",
       x = "library size / ambient depth  (log)", y = NULL) + th

# B. the decoupling is worst at low MOI (one integration dominates the library)
pB <- ggplot(t, aes(moi, depth_cor)) +
  geom_point(aes(colour = log10(moi)), size = 2.6) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 9, colour = "grey30") +
  scale_x_log10() + scale_colour_viridis_c(option = "C", end = 0.9) +
  labs(title = "B. Library size tracks ambient depth only at high MOI",
       subtitle = "corr(log library size, log ambient depth) vs MOI; never high enough to substitute",
       x = "MOI (median assigned guides / cell, log)", y = "corr( libsize , ambient depth )") + th

# C. ambient composition is uneven everywhere (Gini); flat only as a shallow-data artifact
tC <- t[order(t$comp_gini), ]; tC$lab <- factor(tC$lab, levels = tC$lab)
pC <- ggplot(tC, aes(comp_gini, lab)) +
  geom_col(aes(fill = ambient_ok), width = 0.75) +
  scale_fill_manual(values = c(`TRUE` = "#7F77DD", `FALSE` = "grey80")) +
  labs(title = "C. The ambient gRNA pool is uneven (not flat)",
       subtitle = "Gini of ambient composition (0 = uniform, 1 = concentrated); grey = too few ambient guides to estimate",
       x = "Gini of ambient composition", y = NULL) + th

# D. CLEANSER's cheap sub-threshold depth tracks ambient depth, except deep direct-capture
pD <- ggplot(t, aes(lib_med, clean_cor)) +
  geom_point(aes(colour = log10(moi)), size = 2.6) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 9, colour = "grey30") +
  scale_x_log10() + scale_colour_viridis_c(option = "C", end = 0.9) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey55") +
  labs(title = "D. The cheap (sub-threshold) depth works — except very deep data",
       subtitle = "corr(CLEANSER <=2 depth, ambient depth) vs median library size; fixed <=2 threshold breaks when deep",
       x = "median library size (log)", y = "corr( CLEANSER depth , ambient depth )") + th

fig <- (pA | pB) / (pC | pD) + plot_annotation(
  title = "Ambient depth & composition across 17 real CRISPR-screen datasets",
  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave(file.path(D, "survey.png"), fig, width = 14, height = 10, dpi = 150)
cat("wrote", file.path(D, "survey.png"), "\n")
