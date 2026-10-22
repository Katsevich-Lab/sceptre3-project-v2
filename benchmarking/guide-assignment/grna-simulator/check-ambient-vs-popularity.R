#!/usr/bin/env Rscript
# Does a guide's ambient level track how often it is integrated?
#
#   sbatch -c 4 run-in-image.sh check-ambient-vs-popularity.R
#
# For each guide, at the calling threshold t used for the scaling work:
#   n_called      cells with count >= t (popularity proxy)
#   amb_mean      mean count over cells with count < t, zeros included
#   amb_frac1/2   share of those cells with exactly 1 / 2 UMIs
# Cells below t include some weak true signal, so this is rough. Real matrices:
# gasperini (t = 5), replogle-rd7 (t = 10). For comparison, the same statistics on
# one simulated stat-grid dataset (moi 30, s 20, B 20) at t = 5.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/eda/ambient_vs_popularity/:
#   per_guide.csv, summary.csv, scatter.png

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
  library(data.table); library(ggplot2)
})
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "stat-grid-lib.R"))               # simulate_grid()

base    <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
out_dir <- file.path(base, "eda", "ambient_vs_popularity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#' Per-guide popularity and below-threshold ambient statistics (guides x cells).
guide_stats <- function(counts, t, label) {
  counts <- as(counts, "TsparseMatrix")
  d <- data.table(g = counts@i + 1L, x = counts@x)
  N <- ncol(counts)
  s <- d[, .(n_called = sum(x >= t),
             amb_sum  = sum(x[x < t]),
             amb_n1   = sum(x == 1),
             amb_n2   = sum(x == 2)), by = g]
  s <- merge(data.table(g = seq_len(nrow(counts))), s, by = "g", all.x = TRUE)
  for (v in c("n_called", "amb_sum", "amb_n1", "amb_n2")) set(s, which(is.na(s[[v]])), v, 0)
  n_below <- N - s$n_called
  s[, `:=`(dataset = label, threshold = t,
           amb_mean  = amb_sum / n_below,
           amb_frac1 = amb_n1 / n_below,
           amb_frac2 = amb_n2 / n_below)]
  s[]
}

read_real <- function(name) {
  fp <- file.path(base, "input_data", name, "cleanser", "grna_matrix.mtx")
  cat("reading", fp, "\n")
  as(readMM(fp), "CsparseMatrix")
}

res <- list(
  guide_stats(read_real("gasperini"), 5L, "gasperini (real)"),
  guide_stats(read_real("replogle-rd7"), 10L, "replogle-rd7 (real)")
)
invisible(gc())
sim <- simulate_grid(moi = 30, s = 20, B = 20, seed = 107L)$sim
res[[3]] <- guide_stats(assay(sim, "counts"), 5L, "simulated (moi 30, s 20, B 20)")
per_guide <- rbindlist(res)
fwrite(per_guide, file.path(out_dir, "per_guide.csv"))

summ <- per_guide[, .(n_guides          = .N,
                      n_called_median   = as.numeric(median(n_called)),
                      amb_mean_median   = as.numeric(median(amb_mean)),
                      amb_mean_max      = max(amb_mean),
                      spearman_mean     = cor(n_called, amb_mean, method = "spearman"),
                      spearman_frac1    = cor(n_called, amb_frac1, method = "spearman")),
                  by = dataset]
fwrite(summ, file.path(out_dir, "summary.csv"))
print(summ)

fig <- ggplot(per_guide[n_called > 0 & amb_mean > 0], aes(n_called, amb_mean)) +
  geom_point(alpha = 0.3, size = 0.7) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ dataset, scales = "free") +
  labs(x = "Cells at or above the calling threshold (log)",
       y = "Mean count in cells below the threshold (log)") +
  theme_bw(base_size = 10)
tryCatch(ggsave(file.path(out_dir, "scatter.png"), fig, width = 10, height = 3.6, dpi = 110),
         error = function(e) cat("PNG not written:", conditionMessage(e), "\n"))
cat("wrote", out_dir, "\n")
