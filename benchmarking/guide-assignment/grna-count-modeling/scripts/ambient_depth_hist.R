#!/usr/bin/env Rscript
# Histograms of per-cell AMBIENT DEPTH for each dataset, from the cached fishash+ fit.
# Per-cell ambient depth d_c = Cc[c] = colSums of fishash's rank-1 completion of the
# signal-masked matrix (the denoised cell exposure the depth_fix cell margin uses).
# Reads results/ambient_ceiling/fit_cache/<ds>_fit.rds only -- no refit, no matrix load.
suppressPackageStartupMessages({ library(ggplot2) })

OUT   <- "results/ambient_ceiling"; CACHE <- file.path(OUT, "fit_cache")
fits  <- sort(list.files(CACHE, pattern = "_fit\\.rds$", full.names = TRUE))
stopifnot(length(fits) > 0)

rows <- list(); summ <- list()
for (f in fits) {
  fit <- readRDS(f); nm <- fit$dataset
  d   <- as.numeric(fit$Cc)                       # per-cell ambient depth (all cells)
  pos <- d[is.finite(d) & d > 0]                  # drop zero-depth cells (no ambient) for log axis
  rows[[nm]] <- data.frame(dataset = nm, depth = pos)
  summ[[nm]] <- data.frame(dataset = nm, cells = length(d), zero_depth = sum(d == 0),
    med_depth = median(pos), p10 = quantile(pos, .10), p90 = quantile(pos, .90),
    max_depth = max(pos))
}
D  <- do.call(rbind, rows)
S  <- do.call(rbind, summ)
write.csv(S, file.path(OUT, "ambient_depth_summary.csv"), row.names = FALSE)
cat("=== per-cell ambient depth summary ===\n"); print(S, row.names = FALSE)

# order facets by median ambient depth (shallow -> deep)
lev <- S$dataset[order(S$med_depth)]
D$dataset <- factor(D$dataset, levels = lev)
med <- data.frame(dataset = factor(S$dataset, levels = lev), med = S$med_depth)

p_facet <- ggplot(D, aes(depth)) +
  geom_histogram(bins = 40, fill = "#3b6ea5", colour = "white", linewidth = 0.1) +
  geom_vline(data = med, aes(xintercept = med), colour = "#D55E00", linewidth = 0.5) +
  scale_x_log10(breaks = c(1,10,100,1000,10000), labels = c("1","10","100","1k","10k")) +
  facet_wrap(~ dataset, scales = "free", ncol = 4) +
  labs(title = "Per-cell ambient depth by dataset (fishash+ rank-1 cell factor d_c)",
       subtitle = "orange line = dataset median.  x = fitted ambient UMIs per cell (log scale)",
       x = "per-cell ambient depth  d_c", y = "cells") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "ambient_depth_hist_facet.png"), p_facet, width = 14, height = 10, dpi = 150)

# pooled (density per dataset, one overlay) for cross-dataset comparison
p_pool <- ggplot(D, aes(depth, colour = dataset)) +
  geom_density(linewidth = 0.6) +
  scale_x_log10(breaks = c(1,3,10,30,100,300,1000,3000,10000),
                labels = c("1","3","10","30","100","300","1k","3k","10k")) +
  labs(title = "Per-cell ambient depth distributions (fishash+)", x = "per-cell ambient depth  d_c (log)",
       y = "density", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT, "ambient_depth_density_overlay.png"), p_pool, width = 11, height = 6.5, dpi = 150)

cat("\nwrote:\n  ", file.path(OUT, "ambient_depth_hist_facet.png"),
    "\n  ", file.path(OUT, "ambient_depth_density_overlay.png"),
    "\n  ", file.path(OUT, "ambient_depth_summary.csv"), "\n")
