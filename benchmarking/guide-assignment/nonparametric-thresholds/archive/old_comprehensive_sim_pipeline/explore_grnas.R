#!/usr/bin/env Rscript
# Look closely at a few randomly chosen gRNAs: does simple mode-separation on
# log(1+count) put the threshold in the gap between the background and the
# perturbed mode? Ground truth (from the simulation) is overlaid so we can see
# the two true classes directly.
#
# Usage: Rscript explore_grnas.R [dataset_id] [n_grnas] [seed]

suppressPackageStartupMessages(library(Matrix))
HERE     <- dirname(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
DATA_DIR <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(HERE, "..", "guide-assignment-pipeline", "bin", "script", "lib",
                 "threshold_methods.R"))

args       <- commandArgs(trailingOnly = TRUE)
dataset_id <- if (length(args) >= 1) args[1] else "replogle-rd7_simulated"
n_grnas    <- if (length(args) >= 2) as.integer(args[2]) else 6L
seed       <- if (length(args) >= 3) as.integer(args[3]) else 7L
outdir     <- file.path(HERE, "results", dataset_id); dir.create(outdir, FALSE, TRUE)

gm    <- as(readRDS(file.path(DATA_DIR, dataset_id, "sceptre", "grna_matrix.rds")),
            "RsparseMatrix")
gt_fp     <- file.path(DATA_DIR, dataset_id, "true_perturbation_status.rds")
has_truth <- file.exists(gt_fp)
truth     <- if (has_truth) as(readRDS(gt_fp), "RsparseMatrix") else NULL

set.seed(seed)
sel <- sort(sample(seq_len(nrow(gm)), min(n_grnas, nrow(gm))))

png(file.path(outdir, "explore_grnas.png"),
    width = 1300, height = 850, res = 115)
op <- par(mfrow = c(2, ceiling(length(sel) / 2)),
          mar = c(3.4, 3.6, 2.4, 0.8), mgp = c(2.2, 0.6, 0))

cat(sprintf("%-22s %6s %6s | %5s %5s | %5s %5s %5s\n",
            "gRNA", "otsu_t", "vall_t", "Mode1", "Mode2", "sens", "spec", "miss"))
for (g in sel) {
  counts <- as.numeric(gm[g, ])
  tr     <- if (has_truth) as.numeric(truth[g, ]) != 0 else logical(length(counts))
  z      <- log1p(counts)

  ot <- otsu_threshold_log1p(counts)
  vv <- smoothed_valley_threshold(counts)

  # Log-log histogram: x = log(1 + count), y = log10(frequency). The background
  # decays from the count=0/1 spike; the perturbed cells form a separate bump.
  br      <- seq(0, max(z) + 1e-6, length.out = 45)
  h_all   <- hist(z,      breaks = br, plot = FALSE)
  h_pt    <- hist(z[tr],  breaks = br, plot = FALSE)
  yall    <- h_all$counts; yall[yall == 0] <- NA
  ypt     <- h_pt$counts;  ypt[ypt == 0]   <- NA
  ymax    <- max(h_all$counts)
  plot(h_all$mids, yall, type = "h", log = "y", lwd = 3, col = "grey55",
       xlim = c(0, max(z) + 0.5), ylim = c(1, ymax),
       xlab = "log(1 + count)", ylab = "frequency (log scale)", main = "")
  if (has_truth)
    points(h_pt$mids, ypt, type = "h", lwd = 3, col = adjustcolor("firebrick", 0.85))
  if (is.finite(ot$t)) abline(v = log1p(ot$t), col = "darkorange2", lwd = 2)
  if (is.finite(vv$t)) abline(v = log1p(vv$t), col = "royalblue", lwd = 2, lty = 2)
  ttl <- if (has_truth) sprintf("%s  (n_pert=%d)", substr(rownames(gm)[g], 1, 18), sum(tr))
         else sprintf("%s  (nz=%d)", substr(rownames(gm)[g], 1, 18), sum(counts > 0))
  title(main = ttl, cex.main = 0.9)
  legend("topright", bty = "n", cex = 0.7,
         legend = c("all cells", if (has_truth) "perturbed (truth)", "Otsu", "valley"),
         col    = c("grey55", if (has_truth) "firebrick", "darkorange2", "royalblue"),
         lwd = 2, lty = c(1, if (has_truth) 1, 1, 2))

  n_otsu <- if (is.finite(ot$t)) sum(counts >= ot$t) else 0L
  n_vall <- if (is.finite(vv$t)) sum(counts >= vv$t) else 0L
  if (has_truth) {
    pred <- counts >= min(ot$t, Inf)   # report otsu vs truth
    sens <- sum(pred & tr) / max(sum(tr), 1)
    spec <- sum(!pred & !tr) / max(sum(!tr), 1)
    cat(sprintf("%-22s %6s %6s | %5s %5s | sens=%.3f spec=%.3f miss=%d\n",
                substr(rownames(gm)[g], 1, 22), ot$t, vv$t,
                vv$mode1 %||% NA, vv$mode2 %||% NA, sens, spec, sum(tr) - sum(pred & tr)))
  } else {
    cat(sprintf("%-22s %6s %6s | %5s %5s | n_otsu=%d n_valley=%d\n",
                substr(rownames(gm)[g], 1, 22), ot$t, vv$t,
                vv$mode1 %||% NA, vv$mode2 %||% NA, n_otsu, n_vall))
  }
}
par(op); dev.off()
cat("\nWrote", file.path(outdir, "explore_grnas.png"), "\n")
