#!/usr/bin/env Rscript
# EDA for the MOI x ambient grid: per-guide UMI histograms with true perturbation
# status, at a few grid points. Simulates in memory; writes no datasets.
#
#   sbatch run-in-image.sh eda-stat-grid.R
#
# Each grid point gets one figure: one panel per guide, guides chosen at quantiles
# of their number of truly perturbed cells. Every panel at every point shares one
# bin structure, so panels and points compare directly:
#   - exact bins 0..K_EXACT, then bins doubling in width, out to the global max;
#   - all cells counted, zeros included;
#   - ambient (not perturbed) stacked under perturbed, per bin;
#   - y on log10(count + 1).
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/eda/stat_grid/:
#   hist_bins.csv     point, guide, bin, status, count (the histograms, as data)
#   point_summary.csv per point: design values, realized levels
#   hist_<point>.png  one figure per point, if a PNG device is available

suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
  library(ggplot2)
})
source("~/.Rprofile")                                    # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "stat-grid-lib.R"))               # STAT_GRID, simulate_grid(), realized()

# Guides and cells for this EDA; override STAT_GRID here only. At fixed per-cell
# ambient B, per-entry ambient is about B / G, so G sets how close the popular
# guides' ambient comes to losing its mode at 0; perturbed cells per guide are
# about N * A / G.
EDA_N_GUIDES <- 500L
EDA_N_CELLS  <- 50000L
STAT_GRID$n_guides <- EDA_N_GUIDES
STAT_GRID$n_cells  <- EDA_N_CELLS
# Ambient per cell, and the SD of the droplet ambient size factor (log scale).
# A larger d_sigma_drop spreads ambient across cells, which lengthens each guide's
# ambient tail; bypass_inputs() corrects d_mu_drop for it, so E[ambient] stays B.
EDA_B            <- 25
EDA_D_SIGMA_DROP <- 1.0
STAT_GRID_PARAMS$d_sigma_drop <- EDA_D_SIGMA_DROP
# Overdispersion of the ambient (endogenous + exogenous) rate per entry: with
# Phi_noise > 0 each entry's rate is Gamma(1 / Phi_noise, Phi_noise * rate), same
# mean, so the per-entry ambient tail is heavier than Poisson (geometric at 1).
EDA_PHI_NOISE <- 1
STAT_GRID_PARAMS$Phi_noise <- EDA_PHI_NOISE
RISE_MIN_CELLS <- 10L   # a rise over counts 0..3 counts only if the larger bin has this many cells

# MOI x signal per integration, at one ambient level. Ambient must stay plausible
# for every guide (mode at 0, decaying), which B = 20 met at MOI 0.3 and 30 in the
# first EDA; higher B pushed popular guides' ambient off zero.
EDA_POINTS <- expand.grid(s = c(80, 20, 5), moi = c(0.3, 3, 30))
EDA_POINTS$B     <- EDA_B
EDA_POINTS$label <- sprintf("moi%g_s%g", EDA_POINTS$moi, EDA_POINTS$s)
EDA_POINTS$seed  <- 100L + seq_len(nrow(EDA_POINTS))
PLAUSIBLE_MAX_K  <- 3L   # a guide's ambient counts must be non-increasing over 0..3 UMIs
GUIDE_QUANTILES <- seq(0.05, 0.95, length.out = 9)   # of truly perturbed cells per guide
K_EXACT <- 10L
STATUS_LEVELS <- c("Ambient", "Perturbed")
STATUS_COLS   <- c(Ambient = "#00BFC4", Perturbed = "#F8766D")

out_dir <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "eda", "stat_grid",
                     sprintf("G%d_N%d_B%g_sdrop%g_phinoise%g", EDA_N_GUIDES, EDA_N_CELLS, EDA_B,
                             EDA_D_SIGMA_DROP, EDA_PHI_NOISE))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
s <- STAT_GRID$signal_per_integration

#' Exact bins 0..k_exact, then bins doubling in width until max_x is covered.
mixed_bins <- function(max_x, k_exact) {
  upper <- 0:k_exact
  width <- 2L
  while (tail(upper, 1L) < max_x) {
    upper <- c(upper, tail(upper, 1L) + width)
    width <- width * 2L
  }
  lower <- c(0L, head(upper, -1L) + 1L)
  data.frame(bin_id = seq_along(upper), lower = lower, upper = upper,
             label = ifelse(lower == upper, as.character(upper), paste0(lower, "-", upper)))
}

# ---- simulate each point; keep the chosen guides' counts and truth -------------
picked <- list(); summaries <- list()
for (k in seq_len(nrow(EDA_POINTS))) {
  p <- EDA_POINTS[k, ]
  cat(sprintf("\n=== %s: moi %g, s %g, B %g, seed %d ===\n", p$label, p$moi, p$s, p$B, p$seed))
  r <- simulate_grid(p$moi, p$s, p$B, p$seed)
  counts <- assay(r$sim, "counts"); truth <- assay(r$sim, "ground_truth")

  # Plausibility: each guide's ambient (not perturbed) cells have counts that
  # are non-increasing over 0..PLAUSIBLE_MAX_K.
  # Computed from the nonzero entries: an ambient cell is any cell where the guide
  # is not truly integrated; its zeros are the ambient cells with no nonzero entry.
  ct  <- as(counts, "TsparseMatrix")
  amb <- !truth[cbind(ct@i + 1L, ct@j + 1L)]
  gi  <- ct@i[amb] + 1L; xv <- ct@x[amb]
  G   <- nrow(counts); n_amb_cells <- ncol(counts) - rowSums(truth)
  amb_tab <- vapply(seq_len(PLAUSIBLE_MAX_K), function(k) tabulate(gi[xv == k], nbins = G),
                    numeric(G))
  amb_tab <- cbind(n_amb_cells - tabulate(gi, nbins = G), amb_tab)   # column 1: count 0
  # Primary: the ambient mode is above 0 (some count k >= 1 has more cells than 0).
  max_k <- min(max(c(xv, 1)), 50L)
  max_above0 <- do.call(pmax, lapply(seq_len(max_k), function(k) tabulate(gi[xv == k], nbins = G)))
  n_implausible <- sum(max_above0 > amb_tab[, 1])
  # Secondary: a rise over counts 0..3 where the larger bin has RISE_MIN_CELLS cells.
  n_rise <- sum(apply(amb_tab, 1, function(x) any(diff(x) > 0 & x[-1] >= RISE_MIN_CELLS)))
  amb_sum <- vapply(split(xv, factor(gi, levels = seq_len(G))), sum, numeric(1))
  max_amb_mean  <- max(amb_sum / n_amb_cells)

  n_true <- rowSums(truth)
  ord    <- order(n_true)
  idx    <- unique(ord[pmax(1L, round(GUIDE_QUANTILES * length(ord)))])
  picked[[p$label]] <- lapply(idx, function(g) list(
    guide = g, n_true = n_true[g],
    umi = as.integer(counts[g, ]), pert = as.logical(truth[g, ])))

  summaries[[k]] <- cbind(p, data.frame(rho_mean = r$inputs$rho_mean,
                                        expected_total_per_cell = r$inputs$expected_total_per_cell,
                                        nnz_per_cell = nnzero(counts) / ncol(counts),
                                        n_guides_implausible = n_implausible,
                                        n_guides_rise_0to3 = n_rise,
                                        max_guide_ambient_mean = max_amb_mean),
                          realized(r$sim))
  print(signif(Filter(is.numeric, summaries[[k]]), 4), row.names = FALSE)
  rm(r, counts, truth); invisible(gc())
}
summary_df <- do.call(rbind, summaries)
write.csv(summary_df, file.path(out_dir, "point_summary.csv"), row.names = FALSE)

# ---- one bin structure for every panel at every point ----------------------------
max_x <- max(unlist(lapply(picked, function(pl) lapply(pl, function(x) max(x$umi)))))
bins  <- mixed_bins(max_x, K_EXACT)
cat("\nglobal max UMI over the picked guides:", max_x, "->", nrow(bins), "bins\n")

hist_rows <- list()
for (lab in names(picked)) for (x in picked[[lab]]) {
  b <- cut(x$umi, breaks = c(-1, bins$upper), labels = FALSE, right = TRUE)
  for (st in STATUS_LEVELS) {
    keep <- if (st == "Perturbed") x$pert else !x$pert
    cnt  <- tabulate(b[keep], nbins = nrow(bins))
    hist_rows[[length(hist_rows) + 1L]] <- data.frame(
      point = lab, guide = x$guide, n_true = x$n_true, bin_id = bins$bin_id,
      bin = bins$label, status = st, count = cnt)
  }
}
hist_df <- do.call(rbind, hist_rows)
write.csv(hist_df, file.path(out_dir, "hist_bins.csv"), row.names = FALSE)

# ---- figures (best effort: the image may lack a PNG device) ----------------------
logp1 <- function(x) log10(x + 1)
y_breaks <- c(0, 1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4)
x_keep <- bins$bin_id[bins$upper <= K_EXACT & bins$upper %% 2 == 0 | bins$upper > K_EXACT]

for (lab in names(picked)) {
  d <- hist_df[hist_df$point == lab, ]
  d$status <- factor(d$status, levels = STATUS_LEVELS)
  d <- d[order(d$guide, d$bin_id, d$status), ]
  d$y1 <- ave(d$count, d$guide, d$bin_id, FUN = cumsum)    # ambient first, perturbed on top
  d$y0 <- d$y1 - d$count
  d$panel <- sprintf("guide %d (%d perturbed cells)", d$guide, d$n_true)
  d$panel <- factor(d$panel, levels = unique(d$panel[order(d$n_true)]))   # panels in quantile order
  p <- EDA_POINTS[EDA_POINTS$label == lab, ]
  yb <- y_breaks[y_breaks <= max(d$y1)]
  fig <- ggplot(d[d$count > 0, ]) +
    geom_rect(aes(xmin = bin_id - 0.45, xmax = bin_id + 0.45,
                  ymin = logp1(y0), ymax = logp1(y1), fill = status),
              colour = "black", linewidth = 0.15) +
    facet_wrap(~ panel, nrow = 3) +
    scale_fill_manual(values = STATUS_COLS) +
    scale_x_continuous(breaks = x_keep, labels = bins$label[match(x_keep, bins$bin_id)],
                       limits = c(0.4, nrow(bins) + 0.6), expand = c(0, 0)) +
    scale_y_continuous(breaks = logp1(yb), labels = format(yb, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(title = sprintf("MOI %g, signal %g per integration (ambient %g UMIs/cell, sd log droplet ambient %g, Phi_noise %g, %d guides, %s cells)",
                         p$moi, p$s, p$B, EDA_D_SIGMA_DROP, EDA_PHI_NOISE, EDA_N_GUIDES,
                         format(EDA_N_CELLS, big.mark = ",")),
         x = "gRNA UMI count", y = "Number of cells", fill = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          panel.grid.minor = element_blank(), legend.position = "bottom")
  fp <- file.path(out_dir, sprintf("hist_%s.png", lab))
  tryCatch(ggsave(fp, fig, width = 10, height = 7.5, dpi = 110),
           error = function(e) cat("PNG not written for", lab, ":", conditionMessage(e), "\n"))
}
cat("\nwrote", out_dir, "\n")
