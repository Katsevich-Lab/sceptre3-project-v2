#!/usr/bin/env Rscript
# Validate the data-derived regimes by reproducing real-vs-sim per-guide UMI
# histograms for every regime (one panel per dataset).  For each, pick a
# representative bimodal real guide and the sim guide whose above-valley mean is
# closest to it (fair comparison).
# Output: results/sim_framework/regime_histograms_all.png
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(patchwork)})
source(file.path(getwd(), "scripts", "sim_lib.R"))
source(file.path(getwd(), "scripts", "barnyard_io.R"))
GA <- normalizePath(file.path(HERE, ".."))
source(file.path(GA, "grna-simulator", "sims-for-paper.R"))   # plot_umi_histogram_real_vs_sim
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
SURV <- path.expand("~/data/external/perturbseq-survey")
REPRO<- file.path(HERE, "external", "repro_work")
OUT  <- SIMFW()
ch <- read.csv(file.path(OUT, "real_characterization.csv"))

# real-data loaders keyed by the SAME names as the regimes (B__regime__<name>)
loaders <- list()
for (d in list.dirs(SURV, recursive = FALSE)) { f <- file.path(d, "grna_matrix.rds")
  if (file.exists(f)) loaders[[basename(d)]] <- local({ ff <- f; function() readRDS(ff) }) }
loaders[["gasperini"]] <- function() readRDS(file.path(DATA, "gasperini/sceptre/grna_matrix.rds"))
loaders[["replogle"]]  <- function() readRDS(file.path(DATA, "replogle-rd7_medium/sceptre/grna_matrix.rds"))
loaders[["schraivogel"]] <- function() { suppressPackageStartupMessages(library(fishash)); data(crispat_schraivogel); crispat_schraivogel }
for (s in BARN_SAMPLES) { mtx <- file.path(REPRO, paste0(s, "_grna_counts.mtx"))
  if (file.exists(mtx)) loaders[[paste0("barnyard_", sub("mix","",s))]] <- local({ mm <- mtx
    function() as(readMM(mm), "CsparseMatrix") }) }

pick_guide <- function(gm, max_cells = 30000, seed = 1) {
  set.seed(seed)                                       # seed BOTH samples (cells + guides) -- reproducible
  gm <- as(gm, "CsparseMatrix"); if (ncol(gm) > max_cells) gm <- gm[, sort(sample(ncol(gm), max_cells))]
  gmr <- as(gm, "RsparseMatrix"); best <- NULL; bn <- -1
  gi <- if (nrow(gmr) > 400) sort(sample(nrow(gmr), 400)) else seq_len(nrow(gmr))
  for (g in gi) { cv <- as.numeric(gmr[g, ]); if (sum(cv > 0) < 50) next
    v <- smoothed_valley_threshold(cv); if (isTRUE(v$ok)) { n <- sum(cv >= v$t)
      if (n > bn && n < 0.3 * length(cv)) { bn <- n; best <- cv } } }
  if (is.null(best)) best <- as.numeric(gmr[which.max(Matrix::rowSums(gmr > 0)), ]); best
}

ord <- ch$dataset[order(ch$separation)]      # hard -> easy
plots <- list()
for (ds in ord) {
  if (is.null(loaders[[ds]])) next
  real_umi <- tryCatch(pick_guide(loaders[[ds]]()), error = function(e) NULL); if (is.null(real_umi)) next
  rv <- smoothed_valley_threshold(real_umi); real_mode <- if (isTRUE(rv$ok)) rv$mode2 else median(real_umi[real_umi > 0])
  sim <- tryCatch(load_dataset(paste0("B__regime__", ds)), error = function(e) NULL); if (is.null(sim)) next
  smr <- as(sim$counts, "RsparseMatrix"); zr <- as(sim$Z, "RsparseMatrix")
  cand <- which(Matrix::rowSums(zr) > 15)
  if (length(cand)) {
    modes <- vapply(cand, function(g){ cv <- as.numeric(smr[g, ]); v <- smoothed_valley_threshold(cv)
      if (isTRUE(v$ok)) v$mode2 else NA_real_ }, numeric(1))
    ok <- which(is.finite(modes)); g <- if (length(ok)) cand[ok][which.min(abs(modes[ok] - real_mode))] else cand[1]
  } else g <- 1
  sep <- ch$separation[ch$dataset == ds]
  plots[[ds]] <- plot_umi_histogram_real_vs_sim(umis_real = real_umi, umis_sim = as.numeric(smr[g, ]),
      is_pert = as.numeric(zr[g, ]) > 0, title = sprintf("%s  (sep %.1f)", ds, sep[1])) +
    theme(plot.title = element_text(size = 8), legend.position = "none")
}
p <- wrap_plots(plots, ncol = 2)
ggsave(file.path(OUT, "regime_histograms_all.png"), p, width = 16, height = 2.1 * ceiling(length(plots)/2),
       dpi = 120, limitsize = FALSE)
cat("wrote regime_histograms_all.png (", length(plots), "datasets )\n")
