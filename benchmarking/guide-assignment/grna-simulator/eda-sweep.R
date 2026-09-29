#!/usr/bin/env Rscript
# Thesis Fig 1a/1b (the guide and cell "sweeps") for resampled datasets, next to
# the real screen they came from. PROTOTYPE.
#
#   Rscript eda-sweep.R <out-dir-name> <source>=<dataset> [<source>=<dataset> ...]
#   e.g. Rscript eda-sweep.R sweep-weighted-g1000 \
#          gasperini=rs_gasperini_weighted_g1000_c15854 \
#          replogle-rd7=rs_replogle-rd7_weighted_g1000_c231127
#
# Same definitions as assignment-analysis/assignment-for-thesis.Rmd §1 with
# features = "total": the guide (cell) whose total UMIs is closest to each
# quantile q = 0.049, 0.95 (labelled 5th, 95th percentile), all its counts
# including zeros, binned with make_mixed_bin_info() -- one shared bin structure
# per figure. Base R only (the image has no tidyverse).
#
# Writes to <LOCAL_BENCHMARKING_DIR>/guide_assignment/scaling-eda/<out-dir-name>/:
#   sweep_bins.csv   figure (guide|cell), source, kind (real|simulated),
#                    quantile, bin_id, label, count
#   sweep_picks.csv  the picked guide or cell in each panel: its index, total
#                    UMIs, nonzeros, and length (cells for a guide, guides for a cell)

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()
.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- dirname(normalizePath(sub("^--file=", "", .self)))
source(file.path(here, "resample-utils.R"))              # load_real_counts()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: eda-sweep.R <out-dir-name> <source>=<dataset> ...")
pairs <- strsplit(args[-1], "=", fixed = TRUE)
if (any(lengths(pairs) != 2)) stop("each dataset argument must be <source>=<dataset>")
bench   <- .get_config_path("LOCAL_BENCHMARKING_DIR")
in_root <- file.path(bench, "guide_assignment", "input_data")
out_dir <- file.path(bench, "guide_assignment", "scaling-eda", args[1])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Q <- c(0.049, 0.95)
q_label <- function(p) { n <- round(p * 100); paste0(n, if (n %% 100 %in% 11:13) "th" else
  c("th", "st", "nd", "rd", rep("th", 6))[n %% 10 + 1], " percentile") }

# Copied from make_mixed_bin_info() in assignment-analysis/assignment-helpers.R:
# bins 0..k_exact one value each, then widths doubling from 2.
make_mixed_bin_info <- function(max_x, k_exact = 10L) {
  stopifnot(length(max_x) == 1, is.finite(max_x), max_x >= 0)
  upper <- 0:k_exact
  width <- 2L
  while (tail(upper, 1L) < max_x) {
    upper <- c(upper, tail(upper, 1L) + width)
    width <- width * 2L
  }
  lower <- c(0L, head(upper, -1L) + 1L)
  label <- ifelse(lower == upper, as.character(upper), paste0(lower, "-", upper))
  data.frame(bin_id = seq_along(label), lower = lower, upper = upper, label = label,
             stringsAsFactors = FALSE)
}
# pick_by_feature_quantile() from assignment-helpers.R
pick_by_feature_quantile <- function(feature, p)
  which.min(abs(feature - quantile(feature, p, names = FALSE, type = 7)))

panels <- list(); picks <- list()
for (pr in pairs) {
  src <- pr[1]; ds <- pr[2]
  real <- load_real_counts(src, "rds")
  sim  <- readRDS(file.path(in_root, ds, "counts.rds"))
  for (kind in c("real", "simulated")) {
    m <- if (kind == "real") real else sim
    tot_g <- Matrix::rowSums(m); tot_c <- Matrix::colSums(m)
    for (p in Q) {
      gi <- pick_by_feature_quantile(tot_g, p)
      ci <- pick_by_feature_quantile(tot_c, p)
      umi_g <- as.numeric(m[gi, ]); umi_c <- as.numeric(m[, ci])
      panels[[length(panels) + 1]] <- list(figure = "guide", source = src, kind = kind,
                                           quantile = q_label(p), umi = umi_g)
      panels[[length(panels) + 1]] <- list(figure = "cell", source = src, kind = kind,
                                           quantile = q_label(p), umi = umi_c)
      picks[[length(picks) + 1]] <- data.frame(
        figure = c("guide", "cell"), source = src, kind = kind, quantile = q_label(p),
        dataset = if (kind == "real") src else ds, index = c(gi, ci),
        total_umis = c(sum(umi_g), sum(umi_c)), nonzeros = c(sum(umi_g > 0), sum(umi_c > 0)),
        length = c(length(umi_g), length(umi_c)))
    }
  }
  rm(real, sim); invisible(gc())
}

# One shared bin structure per figure, as the thesis figures do.
bins <- do.call(rbind, lapply(c("guide", "cell"), function(fig) {
  ps <- Filter(function(x) x$figure == fig, panels)
  info <- make_mixed_bin_info(max(vapply(ps, function(x) max(x$umi), numeric(1))))
  do.call(rbind, lapply(ps, function(x) {
    cnt <- tabulate(findInterval(x$umi, c(-1, info$upper), left.open = TRUE),
                    nbins = nrow(info))
    data.frame(figure = fig, source = x$source, kind = x$kind, quantile = x$quantile,
               bin_id = info$bin_id, label = info$label, count = cnt)
  }))
}))
picks <- do.call(rbind, picks)
write.csv(bins,  file.path(out_dir, "sweep_bins.csv"),  row.names = FALSE)
write.csv(picks, file.path(out_dir, "sweep_picks.csv"), row.names = FALSE)
print(picks, row.names = FALSE)
cat("\nwrote", out_dir, "\n")
