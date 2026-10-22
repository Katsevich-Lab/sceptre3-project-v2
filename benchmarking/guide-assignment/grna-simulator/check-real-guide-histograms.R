#!/usr/bin/env Rscript
# What do real guides' UMI histograms look like? Is a guide's perturbed component
# a bump away from zero, or a tail decreasing from low counts (as GuideBender's
# signal is with Phi_cell = 1, geometric-like given its mean)?
#
#   sbatch -c 4 run-in-image.sh check-real-guide-histograms.R
#
# Same plot as eda-stat-grid.R: nine guides per dataset at quantiles of their
# number of cells at or above the calling threshold t (guides with at least one),
# bins exact to K_EXACT then doubling, shared within a dataset, log10(count + 1).
# Real data has no ground truth: colour marks below vs at-or-above t, a proxy.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/eda/real_guide_histograms/:
#   hist_<dataset>.png, hist_bins.csv

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2)
})
source("~/.Rprofile")                                    # .get_config_path()

REAL <- list(gasperini = 5L, `replogle-rd7` = 10L)       # dataset = calling threshold
GUIDE_QUANTILES <- seq(0.05, 0.95, length.out = 9)
K_EXACT <- 10L
STATUS_LEVELS <- c("Below threshold", "At or above threshold")
STATUS_COLS   <- setNames(c("#00BFC4", "#F8766D"), STATUS_LEVELS)

base    <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
out_dir <- file.path(base, "eda", "real_guide_histograms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mixed_bins <- function(max_x, k_exact) {
  upper <- 0:k_exact; width <- 2L
  while (tail(upper, 1L) < max_x) { upper <- c(upper, tail(upper, 1L) + width); width <- width * 2L }
  lower <- c(0L, head(upper, -1L) + 1L)
  data.frame(bin_id = seq_along(upper), lower = lower, upper = upper,
             label = ifelse(lower == upper, as.character(upper), paste0(lower, "-", upper)))
}

logp1 <- function(x) log10(x + 1)
y_breaks <- c(0, 1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5)
all_rows <- list()

for (name in names(REAL)) {
  t  <- REAL[[name]]
  fp <- file.path(base, "input_data", name, "cleanser", "grna_matrix.mtx")
  cat("reading", fp, "\n")
  m  <- as(readMM(fp), "CsparseMatrix")                  # guides x cells
  N  <- ncol(m)
  mt <- as(m, "TsparseMatrix")
  d  <- data.table(g = mt@i + 1L, x = as.integer(mt@x))
  n_called <- tabulate(d[x >= t, g], nbins = nrow(m))
  cand <- which(n_called > 0)
  ord  <- cand[order(n_called[cand])]
  idx  <- unique(ord[pmax(1L, round(GUIDE_QUANTILES * length(ord)))])

  sel  <- d[g %in% idx]
  bins <- mixed_bins(max(sel$x), K_EXACT)
  sel[, bin_id := cut(x, breaks = c(-1, bins$upper), labels = FALSE, right = TRUE)]
  sel[, status := ifelse(x >= t, STATUS_LEVELS[2], STATUS_LEVELS[1])]
  h <- sel[, .(count = .N), by = .(g, bin_id, status)]
  zeros <- data.table(g = idx, bin_id = 1L, status = STATUS_LEVELS[1],
                      count = N - tabulate(sel$g, nbins = nrow(m))[idx])   # cells with 0 UMIs
  h <- rbind(h, zeros)
  h <- h[, .(count = sum(count)), by = .(g, bin_id, status)]
  h[, `:=`(dataset = name, threshold = t, n_called = n_called[g], bin = bins$label[bin_id])]
  all_rows[[name]] <- h

  h[, status := factor(status, levels = STATUS_LEVELS)]
  setorder(h, g, bin_id, status)
  h[, y1 := cumsum(count), by = .(g, bin_id)][, y0 := y1 - count]
  h[, panel := sprintf("guide %d (%s cells >= %d)", g, format(n_called, big.mark = ","), t)]
  h[, panel := factor(panel, levels = unique(panel[order(n_called)]))]
  x_keep <- bins$bin_id[(bins$upper <= K_EXACT & bins$upper %% 2 == 0) | bins$upper > K_EXACT]
  yb <- y_breaks[y_breaks <= max(h$y1)]
  fig <- ggplot(h) +
    geom_rect(aes(xmin = bin_id - 0.45, xmax = bin_id + 0.45,
                  ymin = logp1(y0), ymax = logp1(y1), fill = status),
              colour = "black", linewidth = 0.15) +
    facet_wrap(~ panel, nrow = 3) +
    scale_fill_manual(values = STATUS_COLS) +
    scale_x_continuous(breaks = x_keep, labels = bins$label[match(x_keep, bins$bin_id)],
                       limits = c(0.4, nrow(bins) + 0.6), expand = c(0, 0)) +
    scale_y_continuous(breaks = logp1(yb), labels = format(yb, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(title = sprintf("%s (real): %s guides x %s cells, calling threshold %d",
                         name, format(nrow(m), big.mark = ","), format(N, big.mark = ","), t),
         x = "gRNA UMI count", y = "Number of cells", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  tryCatch(ggsave(file.path(out_dir, sprintf("hist_%s.png", name)), fig, width = 10, height = 7.5, dpi = 110),
           error = function(e) cat("PNG not written for", name, ":", conditionMessage(e), "\n"))
  rm(m, mt, d, sel); invisible(gc())
}
fwrite(rbindlist(all_rows, fill = TRUE)[, .(dataset, threshold, g, n_called, bin_id, bin, status, count)],
       file.path(out_dir, "hist_bins.csv"))
cat("wrote", out_dir, "\n")
