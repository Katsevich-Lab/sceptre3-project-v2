#!/usr/bin/env Rscript
# Paired real-vs-sim per-guide histograms via the PER-GUIDE FIT, for every real
# dataset.  For a representative real guide we fit (perturbed fraction, signal
# mean mu, NB dispersion size) on its above-valley cells, then simulate a guide
# of the same #cells from those params + the guide's own ambient (real below-
# valley counts, resampled).  This is the per-guide model each Model B regime
# guide uses; here each panel pairs ONE real guide with a sim from ITS fit.
# Output: results/sim_framework/paired_histograms.png
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(patchwork)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
GA <- normalizePath(file.path(HERE, "..")); source(file.path(GA, "grna-simulator", "sims-for-paper.R"))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
SURV <- path.expand("~/data/external/perturbseq-survey"); REPRO <- file.path(HERE, "external", "repro_work")
OUT  <- SIMFW(); ch <- read.csv(file.path(OUT, "real_characterization.csv"))

loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) { f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) }) }
loaders[["gasperini"]] <- function() readRDS(file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
loaders[["replogle"]]  <- function() readRDS(file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { suppressPackageStartupMessages(library(fishash)); data(crispat_schraivogel); crispat_schraivogel }
for (s in BARN_SAMPLES) { mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (file.exists(mtx)) loaders[[paste0("barnyard_", sub("mix","",s))]] <- local({ mm <- mtx
    function() as(readMM(mm), "CsparseMatrix") }) }

fit_guide <- function(cv) {
  v <- smoothed_valley_threshold(cv); t <- if (isTRUE(v$ok)) v$t else as.numeric(quantile(cv[cv > 0], .9))
  sig <- cv[cv >= t]; amb <- cv[cv < t]; m <- mean(sig); vr <- if (length(sig) > 1) var(sig) else m
  list(frac = mean(cv >= t), mu = m, size = if (vr > m && m > 0) m^2 / (vr - m) else 1e3, amb = amb)
}
sim_guide <- function(p, C) {
  is_sig <- sample(C, round(p$frac * C))
  out <- if (length(p$amb)) sample(p$amb, C, replace = TRUE) else integer(C)   # ambient = real, resampled
  out[is_sig] <- out[is_sig] + rnbinom(length(is_sig), mu = p$mu, size = p$size)
  list(counts = out, is_sig = seq_len(C) %in% is_sig)
}
pick <- function(gm) { gm <- as(gm, "CsparseMatrix"); set.seed(1); if (ncol(gm) > 30000) gm <- gm[, sort(sample(ncol(gm), 30000))]
  gmr <- as(gm, "RsparseMatrix"); bn <- -1; best <- NULL
  for (g in sample(nrow(gmr), min(400, nrow(gmr)))) { cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < 50) next
    v <- smoothed_valley_threshold(cv); if (isTRUE(v$ok)) { n <- sum(cv >= v$t); if (n > bn && n < 0.3 * length(cv)) { bn <- n; best <- cv } } }
  best }

ord <- ch$dataset[order(ch$separation)]; plots <- list()
for (ds in ord) { if (is.null(loaders[[ds]])) next
  cv <- tryCatch(pick(loaders[[ds]]()), error = function(e) NULL); if (is.null(cv)) next
  set.seed(1); p <- fit_guide(cv); s <- sim_guide(p, length(cv))
  sep <- ch$separation[ch$dataset == ds]
  plots[[ds]] <- plot_umi_histogram_real_vs_sim(umis_real = cv, umis_sim = s$counts, is_pert = s$is_sig,
      title = sprintf("%s (sep %.1f): frac=%.3f mu=%.0f size=%.2f", ds, sep[1], p$frac, p$mu, p$size)) +
    theme(plot.title = element_text(size = 7.5), legend.position = "none")
}
ggsave(file.path(OUT, "paired_histograms.png"), wrap_plots(plots, ncol = 2),
       width = 16, height = 2.1 * ceiling(length(plots) / 2), dpi = 120, limitsize = FALSE)
cat("wrote paired_histograms.png (", length(plots), "datasets )\n")
