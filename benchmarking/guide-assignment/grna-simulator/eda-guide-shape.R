#!/usr/bin/env Rscript
# Per-guide count shape, across every guide, for resampled datasets against the
# real screen they came from. PROTOTYPE.
#
#   Rscript eda-guide-shape.R <out-dir-name> <source>=<dataset> [<source>=<dataset> ...]
#
# Resampling mixes each simulated guide's entries from many real guides, which
# should broaden its signal component (thesis Fig 1a, redrawn for simulations).
# These statistics measure that over all guides rather than one per quantile:
#   signal_median   median count among the guide's entries >= threshold
#   signal_iqr      IQR of log2(count) among those entries
#   valley_share    fraction of the guide's nonzeros with 2 <= count < threshold
# signal_median and signal_iqr are computed only for guides with at least
# MIN_SIGNAL signal entries; valley_share for guides with any nonzeros.
#
# Three screens per source: the full real screen, the same real cells the
# simulation drew (real guides, so the cell count matches the simulation's), and
# the simulation.
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/scaling-eda/<out-dir-name>/:
#   shape_hist.csv      histograms over guides, per source, screen and statistic
#   shape_summary.csv   quartiles over guides, per source, screen and statistic

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()
.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "resample-utils.R"))              # load_real_counts(), RESAMPLE_SOURCES

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: eda-guide-shape.R <out-dir-name> <source>=<dataset> ...")
pairs <- strsplit(args[-1], "=", fixed = TRUE)
if (any(lengths(pairs) != 2)) stop("each dataset argument must be <source>=<dataset>")
bench   <- .get_config_path("LOCAL_BENCHMARKING_DIR")
in_root <- file.path(bench, "guide_assignment", "input_data")
out_dir <- file.path(bench, "guide_assignment", "scaling-eda", args[1])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MIN_SIGNAL <- 10L
BREAKS <- list(
  signal_median = 2^seq(0, 16, by = 0.25),   # UMIs, log2-spaced
  signal_iqr    = seq(0, 8, by = 0.2),       # log2 units
  valley_share  = seq(0, 1, by = 0.02))

# Per-guide statistics. Guides are the columns of t(m), so each guide's entries
# are contiguous.
guide_shape <- function(m, threshold) {
  mt <- as(Matrix::t(m), "CsparseMatrix")
  p <- mt@p; x <- mt@x
  res <- t(vapply(seq_len(ncol(mt)), function(g) {
    v <- x[seq.int(p[g] + 1L, length.out = p[g + 1L] - p[g])]
    s <- v[v >= threshold]
    c(signal_median = if (length(s) >= MIN_SIGNAL) stats::median(s) else NA_real_,
      signal_iqr    = if (length(s) >= MIN_SIGNAL) stats::IQR(log2(s)) else NA_real_,
      valley_share  = if (length(v)) mean(v >= 2 & v < threshold) else NA_real_)
  }, numeric(3)))
  as.data.frame(res)
}

hists <- list(); summ <- list()
for (pr in pairs) {
  src <- pr[1]; ds <- pr[2]
  thr <- RESAMPLE_SOURCES[[src]]$threshold
  real <- load_real_counts(src, "rds")
  sim  <- readRDS(file.path(in_root, ds, "counts.rds"))
  cells <- readRDS(file.path(in_root, ds, "sim_params.rds"))$source_cells
  screens <- list(real = real, same_cells = real[, cells], simulated = sim)
  for (sc in names(screens)) {
    st <- guide_shape(screens[[sc]], thr)
    for (stat in names(BREAKS)) {
      v <- st[[stat]]; v <- v[!is.na(v)]
      b <- BREAKS[[stat]]
      h <- tabulate(findInterval(pmin(pmax(v, b[1]), b[length(b)]), b, rightmost.closed = TRUE),
                    nbins = length(b) - 1L)
      hists[[length(hists) + 1]] <- data.frame(source = src, screen = sc, stat = stat,
        bin_lo = head(b, -1), bin_hi = b[-1], frac = h / length(v))
      q <- stats::quantile(v, c(.25, .5, .75), names = FALSE)
      summ[[length(summ) + 1]] <- data.frame(source = src, screen = sc, stat = stat,
        n_guides = length(v), q25 = q[1], median = q[2], q75 = q[3])
    }
    cat(sprintf("%s %-10s: %d guides, %d with >= %d signal entries\n", src, sc,
                nrow(st), sum(!is.na(st$signal_median)), MIN_SIGNAL))
  }
  rm(real, sim, screens); invisible(gc())
}
summ <- do.call(rbind, summ)
write.csv(do.call(rbind, hists), file.path(out_dir, "shape_hist.csv"), row.names = FALSE)
write.csv(summ, file.path(out_dir, "shape_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 4)
cat("\nwrote", out_dir, "\n")
