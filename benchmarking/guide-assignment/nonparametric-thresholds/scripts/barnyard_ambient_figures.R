#!/usr/bin/env Rscript
# Section-1 barnyard figures for the collaborator writeup (all the "ambient is Poisson" figures except
# the per-guide depth-mixing fit, which is barnyard_marginal_fit.R). Produces, into results/collaborator_writeup/:
#   1. barnyard_fig3b_repro.png     — CLEANSER's Fig-3B reproduced (direct-capture vs CROP-seq, 0.90 gate)
#   2. barnyard_frac_mouse.png      — mouse-GEX-fraction per cell on a LOG axis (resolves the human spike
#                                     and its doublet shoulder), with the 0.90 and 0.99 gates
#   3. barnyard_vmr_threshold.png   — per-guide variance vs mean at the 0.90 vs 0.99 gate
#   4. barnyard_offenders.csv       — the residual doublet/chimera cells that survive the 0.99 gate
# Data: external/repro_work/mix0hr_{DirectCapture,Cropseq}_*. Run from nonparametric-thresholds/.
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
REPRO <- "external/repro_work"; OUT <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# load a sample with BASIC QC only (mito/features/depth), keeping frac_homo so we can threshold later
load_basic <- function(sample) {
  gm     <- as(readMM(file.path(REPRO, paste0(sample, "_grna_counts.mtx"))), "CsparseMatrix")
  meta   <- read.csv(file.path(REPRO, paste0(sample, "_meta.csv")))
  guides <- read.csv(file.path(REPRO, paste0(sample, "_guides.csv")))
  sg <- meta$homo_sum_gex + meta$mus_sum_gex; fh <- meta$homo_sum_gex / sg
  qc <- (meta$mito_sum/sg < .15) & (meta$features_gex <= 6000) & (sg <= 20000) &
        (meta$features_gex >= 1500) & (sg >= 3500); qc[is.na(qc)] <- FALSE
  counts <- gm[, qc, drop = FALSE]; rownames(counts) <- guides$guide
  list(counts = counts, guide_homo = guides$guide_type == "homo_guide",
       fh = fh[qc], meta = meta[qc, ])
}
# per wrong-species (mouse) guide: mean & var of UMIs across the given host cells (incl. zeros)
perguide_vmr <- function(o, host) {
  sub <- o$counts[!o$guide_homo, host, drop = FALSE]; nC <- ncol(sub)
  s1 <- Matrix::rowSums(sub); s2 <- Matrix::rowSums(sub * sub); mu <- s1/nC; v <- (s2 - nC*mu^2)/(nC-1)
  data.frame(mean = mu, var = v, vmr = ifelse(mu > 0, v/mu, NA))
}

dc <- load_basic("mix0hr_DirectCapture")
cs <- load_basic("mix0hr_Cropseq")

## ---- 1. Fig-3B reproduction: direct-capture vs CROP-seq (0.90 gate) ---------------
pgA <- rbind(transform(perguide_vmr(dc, dc$fh > 0.90), chem = "direct-capture"),
             transform(perguide_vmr(cs, cs$fh > 0.90), chem = "CROP-seq"))
pgA <- pgA[is.finite(pgA$mean) & pgA$mean > 0, ]
medA <- tapply(pgA$vmr, pgA$chem, median, na.rm = TRUE)
p1 <- ggplot(pgA, aes(mean, var, colour = chem)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.75, size = 2.3) + scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("direct-capture" = "#D55E00", "CROP-seq" = "#0072B2"),
    labels = c(sprintf("CROP-seq (median var/mean %.2f)", medA["CROP-seq"]),
               sprintf("direct-capture (median var/mean %.1f)", medA["direct-capture"]))) +
  labs(title = "Ground-truth ambient dispersion, as CLEANSER reports it (Fig 3B)",
       subtitle = "Each dot = one wrong-species guide (variance vs mean of its UMIs across host cells, incl. zeros); dashed = Poisson.",
       x = "per-guide mean ambient UMIs (log)", y = "per-guide variance (log)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(OUT, "barnyard_fig3b_repro.png"), p1, width = 8, height = 5.2, dpi = 130)

## ---- 2. mouse-GEX-fraction histogram, LOG axis -----------------------------------
mouse_frac <- dc$meta$mus_sum_gex / (dc$meta$homo_sum_gex + dc$meta$mus_sum_gex)
FLOOR <- 3e-4; mf <- pmax(mouse_frac, FLOOR)
region <- ifelse(mouse_frac < 0.01, "human — kept at 0.99 gate",
          ifelse(mouse_frac < 0.10, "human shoulder — kept at 0.90, dropped at 0.99",
          ifelse(mouse_frac < 0.90, "mixed — dropped", "mouse")))
region <- factor(region, levels = c("human — kept at 0.99 gate",
  "human shoulder — kept at 0.90, dropped at 0.99", "mixed — dropped", "mouse"))
brks <- 10^seq(log10(FLOOR), 0, length.out = 70)
p2 <- ggplot(data.frame(mf = mf, region = region), aes(mf, fill = region)) +
  geom_histogram(breaks = brks, colour = NA) +
  geom_vline(xintercept = c(0.10, 0.01), linetype = c("dashed", "solid"),
             colour = c("grey45", "firebrick"), linewidth = c(0.5, 0.8)) +
  annotate("text", x = 0.10, y = Inf, label = "old gate 0.90\n(mouse frac 0.10)", vjust = 1.4, hjust = -0.06, size = 3.1, colour = "grey35") +
  annotate("text", x = 0.01, y = Inf, label = "new gate 0.99\n(mouse frac 0.01)", vjust = 1.4, hjust = 1.06, size = 3.1, colour = "firebrick") +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = c("human — kept at 0.99 gate" = "#0072B2",
    "human shoulder — kept at 0.90, dropped at 0.99" = "#E69F00",
    "mixed — dropped" = "grey70", "mouse" = "#D55E00")) +
  labs(title = "Fraction of GEX UMIs that are mouse, per cell (log axis)",
       subtitle = "Direct-capture 0hr. On a log axis the human cells resolve into a sharp spike at ~0.007 plus a doublet\nshoulder trailing out to 0.10. The old 0.90 gate admits the whole shoulder; the new 0.99 gate removes it.",
       x = "fraction of GEX UMIs that are mouse   mus / (homo + mus)   (log)", y = "cells (sqrt)", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(OUT, "barnyard_frac_mouse.png"), p2, width = 9.5, height = 5.2, dpi = 130)

## ---- 3. per-guide var vs mean: 0.90 gate -> 0.99 gate -> 0.99 minus each guide's top cell -----
# The third panel is a NON-CIRCULAR test of what the 0.99 residual is: an over-dispersed soup is a
# distributed heavy tail (deleting one cell barely moves it); an isolated spike collapses. It collapses.
pg_vmr <- function(o, host, drop_top = FALSE) {
  M <- as.matrix(o$counts[!o$guide_homo, host, drop = FALSE])
  r <- t(apply(M, 1, function(x) {
    if (drop_top) x <- x[-which.max(x)]
    n <- length(x); mu <- mean(x); c(mean = mu, var = (sum(x^2) - n*mu^2)/(n-1)) }))
  data.frame(mean = r[, 1], var = r[, 2], vmr = ifelse(r[, 1] > 0, r[, 2]/r[, 1], NA))
}
P <- list(`GEX gate 0.90` = pg_vmr(dc, dc$fh > 0.90),
          `GEX gate 0.99` = pg_vmr(dc, dc$fh > 0.99),
          `GEX gate 0.99, minus each guide's single largest cell` = pg_vmr(dc, dc$fh > 0.99, drop_top = TRUE))
pgV <- do.call(rbind, Map(function(df, nm) transform(df, panel = nm), P, names(P)))
pgV <- pgV[is.finite(pgV$mean) & pgV$mean > 0 & is.finite(pgV$var) & pgV$var > 0, ]
medV <- tapply(pgV$vmr, pgV$panel, median, na.rm = TRUE)
pgV$panel <- factor(pgV$panel, levels = names(P),
  labels = sprintf("%s\n(median var/mean %.1f)", names(P), medV[names(P)]))
p3 <- ggplot(pgV, aes(mean, var)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, size = 1.5, colour = "#0072B2") +
  facet_wrap(~panel, nrow = 1) + scale_x_log10() + scale_y_log10() +
  labs(title = "The 0.99-gate residual is isolated single-cell spikes, not an over-dispersed soup",
       subtitle = "Each dot = one wrong-species guide (CLEANSER's Fig-3B statistic); dashed = Poisson. Tightening 0.90 → 0.99 removes most of the excess;\ndeleting each guide's ONE largest cell collapses the rest onto the line (a distributed heavy-tailed soup would not) — so the residual is spikes.",
       x = "per-guide mean ambient UMIs (log)", y = "per-guide variance (log)") +
  theme_bw(base_size = 10.5)
ggsave(file.path(OUT, "barnyard_vmr_threshold.png"), p3, width = 12, height = 4.6, dpi = 130)

## ---- 4. residual offenders surviving the 0.99 gate -------------------------------
hc99 <- which(dc$fh > 0.99); mg <- which(!dc$guide_homo)
W <- dc$counts[mg, hc99, drop = FALSE]; top <- apply(W, 2, max); which_g <- apply(W, 2, which.max)
off <- which(top >= 10)
od <- data.frame(guide = rownames(dc$counts)[mg][which_g[off]], wrong_species_count = top[off],
  frac_homo_GEX = round(dc$fh[hc99][off], 4), mouse_mRNA_UMIs = dc$meta$mus_sum_gex[hc99][off],
  human_mRNA_UMIs = dc$meta$homo_sum_gex[hc99][off], total_gRNA = Matrix::colSums(dc$counts)[hc99][off])
od <- od[order(-od$wrong_species_count), ]
write.csv(od, file.path(OUT, "barnyard_offenders.csv"), row.names = FALSE)

## ---- reported numbers ------------------------------------------------------------
clean <- top < 3
cat("=== numbers for the writeup ===\n")
cat(sprintf("fig-3B median var/mean: direct-capture %.1f, CROP-seq %.2f\n", medA["direct-capture"], medA["CROP-seq"]))
cat(sprintf("per-guide median var/mean: 0.90 -> %.1f ; 0.99 -> %.1f\n", medV[lab90], medV[lab99]))
cat(sprintf("cells kept: 0.90 -> %d ; 0.99 -> %d (%.1f%% dropped)\n",
            sum(dc$fh > 0.9), sum(dc$fh > 0.99), 100*(1 - sum(dc$fh > 0.99)/sum(dc$fh > 0.9))))
cat(sprintf("pooled var/mean @0.99 = %.1f ; offenders (>=10) = %d\n",
            {y <- as.numeric(W); var(y)/mean(y)}, length(off)))
cat(sprintf("human frac_homo: median %.4f, max %.4f ; clean-cell mouse mRNA median %d [%d-%d]\n",
            median(dc$fh[hc99]), max(dc$fh[hc99]),
            median(dc$meta$mus_sum_gex[hc99][clean]),
            round(quantile(dc$meta$mus_sum_gex[hc99][clean], .25)),
            round(quantile(dc$meta$mus_sum_gex[hc99][clean], .75))))
cat("\noffenders:\n"); print(od, row.names = FALSE)
cat("\nwrote figures + offenders.csv to", OUT, "\n")
