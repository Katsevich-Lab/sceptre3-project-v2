#!/usr/bin/env Rscript
# What do the fishash authors' simulated guides look like? Per-guide UMI histograms
# and an ambient plausibility check for their scenarios, simulated through their
# wrapper with their arguments (fishash_analysis simulate/*.R). Writes no datasets.
#
#   sbatch -c 4 run-in-image.sh eda-authors-scenarios.R
#
# Scenarios (20,000 cells, hurdle_prob 0.1, chunk_cells 1000, one dataset each):
#   varyNumGuides high expression: snr 4, count_per_cell 100, frac_noise_endo .75,
#     moi .3, G in 20, 200, 2000, 20000
#   varyNumGuides low expression:  snr 1, count_per_cell 20,  frac_noise_endo .25,
#     moi .3, same G
#   varyMOI: high expression, G 200, moi in .1, 1, 10 (moi .3 is the G 200 point above)
# The authors' 80,000-guide tier (their runtime tier) is skipped.
#
# Plausibility, per guide, over its ambient (not truly integrated) cells:
#   mode_nonzero  some count above 0 holds more cells than count 0 does
#   rises_0to3    the counts at 0..3 UMIs increase somewhere, where the larger bin
#                 has at least RISE_MIN_CELLS cells (ignores few-cell noise)
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/eda/authors_scenarios/:
#   point_summary.csv, hist_bins.csv, hist_<point>.png

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
  library(data.table); library(ggplot2)
})
source("~/.Rprofile")                                    # .get_config_path()

hi <- list(snr = 4, count_per_cell = 100, frac_noise_endo = 0.75)
lo <- list(snr = 1, count_per_cell = 20,  frac_noise_endo = 0.25)
POINTS <- rbind(
  data.frame(scenario = "high", n_guides = c(20, 200, 2000, 20000), moi = 0.3),
  data.frame(scenario = "low",  n_guides = c(20, 200, 2000, 20000), moi = 0.3),
  data.frame(scenario = "high", n_guides = 200, moi = c(0.1, 1, 10))
)
POINTS$label <- sprintf("%s_G%d_moi%g", POINTS$scenario, POINTS$n_guides, POINTS$moi)
POINTS$seed  <- 200L + seq_len(nrow(POINTS))

GUIDE_QUANTILES <- seq(0.05, 0.95, length.out = 9)
K_EXACT         <- 10L
RISE_MIN_CELLS  <- 10L
STATUS_LEVELS   <- c("Ambient", "Perturbed")
STATUS_COLS     <- c(Ambient = "#00BFC4", Perturbed = "#F8766D")

out_dir <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "eda",
                     "authors_scenarios")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

simulate_authors <- function(p) {
  reg <- if (p$scenario == "high") hi else lo
  set.seed(p$seed)
  fishash::simulate_guidebender2(
    n_guides = p$n_guides, n_cells = 20000, moi = p$moi, hurdle_prob = 0.1,
    snr = reg$snr, count_per_cell = reg$count_per_cell,
    frac_noise_endo = reg$frac_noise_endo,
    return_sparse_only = TRUE, chunk_cells = 1000)
}

#' Per-guide ambient plausibility, from the nonzero entries.
ambient_stats <- function(counts, truth) {
  ct  <- as(counts, "TsparseMatrix")
  amb <- !truth[cbind(ct@i + 1L, ct@j + 1L)]
  d   <- data.table(g = ct@i[amb] + 1L, x = ct@x[amb])
  G   <- nrow(counts)
  n_amb <- ncol(counts) - rowSums(truth)                  # ambient cells per guide
  nz  <- d[, .N, by = g]
  zeros <- n_amb; zeros[nz$g] <- n_amb[nz$g] - nz$N        # ambient cells at count 0
  top <- d[, .N, by = .(g, x)][, .(max_n = max(N)), by = g]
  mode_nonzero <- logical(G); mode_nonzero[top$g] <- top$max_n > zeros[top$g]
  c13 <- dcast(d[x <= 3, .N, by = .(g, x)], g ~ x, value.var = "N", fill = 0)
  tab <- matrix(0, G, 4); tab[, 1] <- zeros
  for (k in 1:3) if (as.character(k) %in% names(c13)) tab[c13$g, k + 1] <- c13[[as.character(k)]]
  rises <- apply(tab, 1, function(v) any(diff(v) > 0 & v[-1] >= RISE_MIN_CELLS))
  amb_sum <- numeric(G); s <- d[, .(s = sum(x)), by = g]; amb_sum[s$g] <- s$s
  list(mode_nonzero = mode_nonzero, rises = rises, amb_mean = amb_sum / n_amb)
}

mixed_bins <- function(max_x, k_exact) {
  upper <- 0:k_exact; width <- 2L
  while (tail(upper, 1L) < max_x) { upper <- c(upper, tail(upper, 1L) + width); width <- width * 2L }
  lower <- c(0L, head(upper, -1L) + 1L)
  data.frame(bin_id = seq_along(upper), lower = lower, upper = upper,
             label = ifelse(lower == upper, as.character(upper), paste0(lower, "-", upper)))
}

picked <- list(); summaries <- list()
for (k in seq_len(nrow(POINTS))) {
  p <- POINTS[k, ]
  cat(sprintf("\n=== %s (seed %d) ===\n", p$label, p$seed))
  sim <- simulate_authors(p)
  counts <- assay(sim, "counts"); truth <- assay(sim, "ground_truth")
  signal <- assay(sim, "counts_signal")
  a <- ambient_stats(counts, truth)

  n_true <- rowSums(truth)
  ord <- order(n_true)
  idx <- unique(ord[pmax(1L, round(GUIDE_QUANTILES * length(ord)))])
  picked[[p$label]] <- lapply(idx, function(g) list(
    guide = g, n_true = n_true[g], umi = as.integer(counts[g, ]), pert = as.logical(truth[g, ])))

  totals <- colSums(counts); true_cnt <- counts[truth]
  summaries[[k]] <- cbind(p, data.frame(
    true_entries_per_cell     = sum(truth) / ncol(counts),
    perturbed_cells_per_guide = mean(n_true),
    signal_per_true_entry     = sum(signal) / sum(truth),
    ambient_per_cell          = mean(totals - colSums(signal)),
    ambient_per_entry         = mean(totals - colSums(signal)) / p$n_guides,
    max_guide_ambient_mean    = max(a$amb_mean),
    frac_guides_mode_nonzero  = mean(a$mode_nonzero),
    frac_guides_rise_0to3     = mean(a$rises),
    umis_per_cell_median      = median(totals),
    frac_true_count_lt2       = mean(true_cnt < 2)))
  print(signif(Filter(is.numeric, summaries[[k]]), 4), row.names = FALSE)
  rm(sim, counts, truth, signal); invisible(gc())
}
summary_df <- do.call(rbind, summaries)
fwrite(summary_df, file.path(out_dir, "point_summary.csv"))

# One bin structure per point (scenarios differ in scale by orders of magnitude).
hist_rows <- list()
for (lab in names(picked)) {
  bins <- mixed_bins(max(vapply(picked[[lab]], function(x) max(x$umi), numeric(1))), K_EXACT)
  for (x in picked[[lab]]) {
    b <- cut(x$umi, breaks = c(-1, bins$upper), labels = FALSE, right = TRUE)
    for (st in STATUS_LEVELS) {
      keep <- if (st == "Perturbed") x$pert else !x$pert
      hist_rows[[length(hist_rows) + 1L]] <- data.frame(
        point = lab, guide = x$guide, n_true = x$n_true, bin_id = bins$bin_id,
        bin = bins$label, status = st, count = tabulate(b[keep], nbins = nrow(bins)))
    }
  }
}
hist_df <- rbindlist(hist_rows)
fwrite(hist_df, file.path(out_dir, "hist_bins.csv"))

logp1 <- function(x) log10(x + 1)
y_breaks <- c(0, 1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4)
for (lab in names(picked)) {
  d <- as.data.frame(hist_df[point == lab])
  nb <- max(d$bin_id)
  lab_of <- unique(d[, c("bin_id", "bin")]); lab_of <- lab_of[order(lab_of$bin_id), ]
  x_keep <- lab_of$bin_id[(lab_of$bin_id <= K_EXACT + 1 & (lab_of$bin_id - 1) %% 2 == 0) | lab_of$bin_id > K_EXACT + 1]
  d$status <- factor(d$status, levels = STATUS_LEVELS)
  d <- d[order(d$guide, d$bin_id, d$status), ]
  d$y1 <- ave(d$count, d$guide, d$bin_id, FUN = cumsum); d$y0 <- d$y1 - d$count
  d$panel <- sprintf("guide %d (%d perturbed cells)", d$guide, d$n_true)
  d$panel <- factor(d$panel, levels = unique(d$panel[order(d$n_true)]))
  p <- POINTS[POINTS$label == lab, ]
  yb <- y_breaks[y_breaks <= max(d$y1)]
  fig <- ggplot(d[d$count > 0, ]) +
    geom_rect(aes(xmin = bin_id - 0.45, xmax = bin_id + 0.45,
                  ymin = logp1(y0), ymax = logp1(y1), fill = status),
              colour = "black", linewidth = 0.15) +
    facet_wrap(~ panel, nrow = 3) +
    scale_fill_manual(values = STATUS_COLS) +
    scale_x_continuous(breaks = x_keep, labels = lab_of$bin[match(x_keep, lab_of$bin_id)],
                       limits = c(0.4, nb + 0.6), expand = c(0, 0)) +
    scale_y_continuous(breaks = logp1(yb), labels = format(yb, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(title = sprintf("fishash authors' %s-expression scenario: %d guides, MOI %g",
                         p$scenario, p$n_guides, p$moi),
         x = "gRNA UMI count", y = "Number of cells", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  tryCatch(ggsave(file.path(out_dir, sprintf("hist_%s.png", lab)), fig, width = 10, height = 7.5, dpi = 110),
           error = function(e) cat("PNG not written for", lab, ":", conditionMessage(e), "\n"))
}
cat("\nwrote", out_dir, "\n")
