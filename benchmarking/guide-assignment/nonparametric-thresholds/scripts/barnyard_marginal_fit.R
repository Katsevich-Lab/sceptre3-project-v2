#!/usr/bin/env Rscript
# Fit our depth_fix ambient model on the barnyard direct-capture data and show it reproduces
# the per-guide marginal count histograms.
#
#   (A) wrong-species guides (ground-truth PURE ambient): the WHOLE marginal is ambient, so a
#       correct ambient model must reproduce it end-to-end. We overlay:
#         - observed cell counts at k = 0,1,2,3,>=4
#         - our model: depth-mixed Poisson with the denoised rank-1 rate a_g d_c  (depth_fix)
#         - baseline: a single homogeneous Poisson at the guide's mean  (ignores depth variation)
#   (B) right-species guides (real integrations): the ambient Poisson should capture the LOW mode;
#       the high mode is signal (integration), left unmodeled here.
#
# The denoised rate field is taken from our real method: run contingency_assign (depth_fix), then
# lambda_{gc} = R_g C_c / T from fishash's rank-1 completion of the signal-masked matrix.
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(fishash) })
source("scripts/contingency_method.R")

REPRO <- "external/repro_work"; OUT <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

## ---- load de-doubleted, GEX-pure barnyard direct-capture -----------------------
barnyard_qc <- function(meta, purity = 0.9) {
  sg <- meta$homo_sum_gex + meta$mus_sum_gex; fh <- meta$homo_sum_gex / sg
  qc <- (meta$mito_sum/sg < .15) & (meta$features_gex <= 6000) & (sg <= 20000) &
        (meta$features_gex >= 1500) & (sg >= 3500) & (fh < (1-purity) | fh > purity)
  qc[is.na(qc)] <- FALSE; list(qc = qc, frac_homo = fh)
}
gm     <- as(readMM(file.path(REPRO, "mix0hr_DirectCapture_grna_counts.mtx")), "CsparseMatrix")
meta   <- read.csv(file.path(REPRO, "mix0hr_DirectCapture_meta.csv"))
guides <- read.csv(file.path(REPRO, "mix0hr_DirectCapture_guides.csv"))
q  <- barnyard_qc(meta, 0.90)
counts <- gm[, q$qc, drop = FALSE]; rownames(counts) <- guides$guide
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
guide_homo <- guides$guide_type == "homo_guide"
cell_homo  <- q$frac_homo[q$qc] > 0.90
# gRNA-channel doublet removal (drop co-encapsulated wrong-species cells)
amb_mask <- as(Matrix(outer(guide_homo, cell_homo, `!=`), sparse = TRUE), "lgCMatrix")
wrong_frac <- Matrix::colSums(counts * amb_mask) / pmax(Matrix::colSums(counts), 1)
keep <- wrong_frac <= 0.005    # strict doublet cut -> clean Poisson ground-truth ambient
counts <- counts[, keep, drop = FALSE]; cell_homo <- cell_homo[keep]
cat(sprintf("cohort: %d guides x %d cells (%d human, %d mouse)\n",
            nrow(counts), ncol(counts), sum(cell_homo), sum(!cell_homo)))

## ---- fit our method (depth_fix) and extract the denoised rank-1 rate a_g d_c ----
res <- contingency_assign(counts, q = 0.05, refit = 10, min_count = 2,
                          cell_margin = "ambient", tail = "hyper", fdr = "GS")
bg  <- fishash::impute_masked_counts(counts, res$assigned)   # rank-1 Poisson completion (signal masked)
Rg  <- Matrix::rowSums(bg); Cc <- Matrix::colSums(bg); Tn <- sum(bg)
lam_of <- function(g, cells) (Rg[g] / Tn) * Cc[cells]        # lambda_{gc} = a_g d_c for guide g

KMAX <- 4L
# predicted vs observed CELL COUNTS at k = 0..KMAX-1, >=KMAX, for guide g over a cell set
marg <- function(g, cells) {
  y   <- as.numeric(counts[g, cells])
  lam <- lam_of(g, cells)
  ks  <- 0:KMAX
  obs <- vapply(ks, function(k) if (k < KMAX) sum(y == k) else sum(y >= k), numeric(1))
  # our model: depth-mixed Poisson (sum of per-cell Poisson pmfs at a_g d_c)
  pm  <- vapply(ks, function(k) if (k < KMAX) sum(dpois(k, lam)) else sum(1 - ppois(k-1, lam)), numeric(1))
  # baseline: single homogeneous Poisson at the guide's mean
  mu  <- mean(y); n <- length(cells)
  ps  <- vapply(ks, function(k) if (k < KMAX) n*dpois(k, mu) else n*(1 - ppois(k-1, mu)), numeric(1))
  data.frame(guide = g, k = ks, klab = ifelse(ks < KMAX, as.character(ks), paste0(KMAX, "+")),
             observed = obs, depthmix = pm, single = ps,
             vmr = var(y)/mean(y), mean = mu)
}

library(tidyr)
mouse_g <- rownames(counts)[!guide_homo]           # wrong-species (never-integrated) guides
host    <- which(cell_homo)                         # human host cells

## ---- Per-guide NEGATIVE-CONTROL (wrong-species) marginals ------------------------
# A wrong-species guide can never integrate into a host cell, so its ENTIRE marginal is
# ambient soup (no signal mode). A correct ambient model must reproduce the whole thing.
# We fit each guide two ways and overlay on the observed counts:
#   - our model : depth-mixed Poisson with the denoised rank-1 rate a_g d_c   (depth_fix)
#   - baseline  : a single homogeneous Poisson at the guide's mean            (ignores depth spread)
# Pick the 6 highest-ambient-share guides (richest histograms); all are clean at this cut.
amb_tot <- vapply(mouse_g, function(g) sum(as.numeric(counts[g, host])), numeric(1))
sel_g   <- mouse_g[order(-amb_tot)][1:6]
cat("selected guides (max count should be small = clean, no doublets):\n")
print(data.frame(guide = sel_g,
                 mean = round(vapply(sel_g, function(g) mean(as.numeric(counts[g,host])), numeric(1)),3),
                 maxct = vapply(sel_g, function(g) max(as.numeric(counts[g,host])), numeric(1))), row.names = FALSE)

D <- do.call(rbind, lapply(sel_g, marg, cells = host))   # marg() uses KMAX bins with a "KMAX+" top bin
D$facet <- factor(D$guide, levels = sel_g,
  labels = sprintf("%s   (mean %.3f, var/mean %.2f)",
                   sel_g, D$mean[match(sel_g, D$guide)], D$vmr[match(sel_g, D$guide)]))
klev <- c(as.character(0:(KMAX-1)), paste0(KMAX, "+"))
D$klab <- factor(D$klab, levels = klev)
Dlong <- pivot_longer(D, c(depthmix, single), names_to = "model", values_to = "pred")
Dlong$model <- factor(Dlong$model, c("depthmix","single"),
  labels = c("our model: depth-mixed Poisson (a_g d_c)", "single Poisson (guide mean)"))
Dlong$klab <- factor(Dlong$klab, levels = klev)

pA <- ggplot(D, aes(klab, observed)) +
  geom_col(fill = "grey82", width = 0.68) +
  geom_line(data = Dlong, aes(klab, pred, colour = model, group = model), linewidth = 0.7) +
  geom_point(data = Dlong, aes(klab, pred, colour = model), size = 2.6) +
  facet_wrap(~facet, ncol = 3) +
  scale_y_log10() +
  scale_colour_manual(values = c("our model: depth-mixed Poisson (a_g d_c)" = "#D55E00",
                                 "single Poisson (guide mean)" = "#0072B2")) +
  labs(title = "The fishash+ ambient model fits the per-guide wrong-species (pure-ambient) marginals",
       subtitle = "Barnyard direct-capture: wrong-species (never-integrated) guides across host cells. Bars = observed cells; lines = model.\nOnly ambient — no signal mode. The single Poisson misses the count-2 cells (cell-depth spread); the depth-mixed model captures them.",
       x = "gRNA UMI count", y = "number of cells (log)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "barnyard_marginal_wrongspecies.png"), pA, width = 10, height = 6, dpi = 130)

## ---- goodness-of-fit table (per guide) ------------------------------------------
gof <- do.call(rbind, lapply(sel_g, function(g) {
  d <- marg(g, host)
  chi_dm <- sum((d$observed - d$depthmix)^2 / pmax(d$depthmix, 1e-6))
  chi_s  <- sum((d$observed - d$single  )^2 / pmax(d$single,   1e-6))
  data.frame(guide = g, mean_ambient = round(d$mean[1], 3), marginal_vmr = round(d$vmr[1], 2),
             chisq_depthmix = round(chi_dm, 1), chisq_single = round(chi_s, 1))
}))
write.csv(gof, file.path(OUT, "barnyard_marginal_gof.csv"), row.names = FALSE)
print(gof)
cat("wrote figure + gof to", OUT, "\n")
