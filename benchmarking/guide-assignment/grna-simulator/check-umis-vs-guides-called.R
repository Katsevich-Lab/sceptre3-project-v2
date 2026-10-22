#!/usr/bin/env Rscript
# Do a cell's gRNA UMIs grow with the number of guides it carries, or is there a
# fixed per-cell budget shared among them?
#
#   sbatch -c 4 run-in-image.sh check-umis-vs-guides-called.R
#
# simulate_guidebender2 (the fishash authors' wrapper) holds a cell's expected
# total fixed and divides the signal among its integrations, so UMIs per guide
# fall as about 1/k. The stat-grid design holds signal per integration fixed, so
# totals grow with k. Per cell, at the calling threshold t of the scaling work:
#   k            guides called (count >= t)
#   total        all gRNA UMIs
#   called       UMIs of the called guides (signal proxy); called / k per guide
#   below        UMIs of entries below t (ambient proxy)
# k comes from a threshold, so deeper cells can call more guides; UMIs per called
# guide is the less threshold-sensitive readout.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/eda/umis_vs_guides_called/:
#   by_k.csv, plot.png

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(ggplot2)
})
source("~/.Rprofile")                                    # .get_config_path()

base    <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
out_dir <- file.path(base, "eda", "umis_vs_guides_called")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

REAL <- list(gasperini = 5L, `replogle-rd7` = 10L)       # dataset = calling threshold

per_cell <- function(name, t) {
  fp <- file.path(base, "input_data", name, "cleanser", "grna_matrix.mtx")
  cat("reading", fp, "\n")
  m <- as(readMM(fp), "TsparseMatrix")                   # guides x cells
  d <- data.table(cell = m@j + 1L, x = m@x)
  pc <- d[, .(k = sum(x >= t), total = sum(x), called = sum(x[x >= t]),
              below = sum(x[x < t])), by = cell]
  pc <- merge(data.table(cell = seq_len(ncol(m))), pc, by = "cell", all.x = TRUE)
  for (v in c("k", "total", "called", "below")) set(pc, which(is.na(pc[[v]])), v, 0)
  pc[, dataset := name]
  pc[]
}

pc <- rbindlist(lapply(names(REAL), function(n) per_cell(n, REAL[[n]])))
by_k <- pc[, .(n_cells             = .N,
               total_median        = as.numeric(median(total)),
               total_mean          = mean(total),
               called_per_guide_median = if (k[1] > 0) as.numeric(median(called / k)) else NA_real_,
               called_per_guide_mean   = if (k[1] > 0) mean(called / k) else NA_real_,
               below_median        = as.numeric(median(below)),
               below_mean          = mean(below)),
           by = .(dataset, k)][order(dataset, k)]
fwrite(by_k, file.path(out_dir, "by_k.csv"))
print(by_k[n_cells >= 50], nrows = 200)

# log-log slopes over k >= 1 bins with at least 50 cells, weighted by cell count
fit <- by_k[k >= 1 & n_cells >= 50, {
  w <- n_cells
  list(slope_total            = coef(lm(log(total_median) ~ log(k), weights = w))[2],
       slope_called_per_guide = coef(lm(log(called_per_guide_median) ~ log(k), weights = w))[2])
}, by = dataset]
cat("\nlog-log slopes vs k (1 = totals grow in proportion to k; 0 = flat):\n"); print(fit)

long <- melt(by_k[k >= 1 & n_cells >= 50],
             id.vars = c("dataset", "k", "n_cells"),
             measure.vars = c("total_median", "called_per_guide_median", "below_median"),
             variable.name = "quantity", value.name = "umis")
long[, quantity := factor(quantity, c("total_median", "called_per_guide_median", "below_median"),
                          c("Total gRNA UMIs", "UMIs per called guide", "UMIs below threshold"))]
fig <- ggplot(long, aes(k, umis, colour = quantity)) +
  geom_line() + geom_point(aes(size = n_cells), alpha = 0.6) +
  scale_x_log10() + scale_y_log10() + scale_size_area(max_size = 3, guide = "none") +
  facet_wrap(~ dataset, scales = "free") +
  labs(x = "Guides called per cell (log)", y = "Median UMIs per cell (log)", colour = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
tryCatch(ggsave(file.path(out_dir, "plot.png"), fig, width = 9, height = 4, dpi = 110),
         error = function(e) cat("PNG not written:", conditionMessage(e), "\n"))
cat("wrote", out_dir, "\n")
