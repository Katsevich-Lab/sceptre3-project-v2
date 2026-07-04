#!/usr/bin/env Rscript
# Standalone barnyard direct-capture ambient-Poisson figures for the collaborator writeup.
# Reproduces CLEANSER's Fig-3B (per-guide var-vs-mean of ground-truth ambient) and shows the
# de-doubleting collapse onto the Poisson line. Logic transcribed from doublet_overdispersion.qmd.
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })

REPRO <- "external/repro_work"
OUT   <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

barnyard_qc <- function(meta, purity = 0.9) {
  sum_gex   <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_homo <- meta$homo_sum_gex / sum_gex
  qc <- (meta$mito_sum / sum_gex < .15) & (meta$features_gex <= 6000) &
        (sum_gex <= 20000) & (meta$features_gex >= 1500) & (sum_gex >= 3500) &
        (frac_homo < (1 - purity) | frac_homo > purity)
  qc[is.na(qc)] <- FALSE
  list(qc = qc, frac_homo = frac_homo)
}
load_sample <- function(sample, purity = 0.9) {
  gm     <- as(readMM(file.path(REPRO, paste0(sample, "_grna_counts.mtx"))), "CsparseMatrix")
  meta   <- read.csv(file.path(REPRO, paste0(sample, "_meta.csv")))
  guides <- read.csv(file.path(REPRO, paste0(sample, "_guides.csv")))
  stopifnot(nrow(meta) == ncol(gm), nrow(guides) == nrow(gm))
  q <- barnyard_qc(meta, purity = purity)
  counts <- gm[, q$qc, drop = FALSE]; rownames(counts) <- guides$guide
  list(counts = counts, guide_homo = guides$guide_type == "homo_guide",
       cell_homo = q$frac_homo[q$qc] > purity,
       chemistry = if (grepl("Cropseq", sample)) "CROP-seq" else "direct-capture")
}
amb_mask_of <- function(o) as(Matrix(outer(o$guide_homo, o$cell_homo, `!=`), sparse = TRUE), "lgCMatrix")

# per wrong-species guide: mean & var of UMIs across host cells (incl. zeros)
perguide <- function(o, cell_subset = NULL) {
  gsel <- !o$guide_homo; csel <- o$cell_homo
  if (!is.null(cell_subset)) csel <- csel & cell_subset
  sub <- o$counts[gsel, csel, drop = FALSE]; nC <- ncol(sub)
  s1 <- Matrix::rowSums(sub); s2 <- Matrix::rowSums(sub * sub)
  mu <- s1 / nC; v <- (s2 - nC * mu^2) / (nC - 1)
  data.frame(mean = mu, var = v, vmr = ifelse(mu > 0, v / mu, NA))
}

dc <- load_sample("mix0hr_DirectCapture", 0.90)
cs <- load_sample("mix0hr_Cropseq", 0.90)
wrong_frac <- Matrix::colSums(dc$counts * amb_mask_of(dc)) / pmax(Matrix::colSums(dc$counts), 1)

# ---- Figure 1: Fig-3B reproduction (direct-capture vs CROP-seq) ------------------
pg_dc <- transform(perguide(dc), chem = "direct-capture")
pg_cs <- transform(perguide(cs), chem = "CROP-seq")
pgA <- rbind(pg_dc, pg_cs); pgA <- pgA[is.finite(pgA$mean) & pgA$mean > 0, ]
med <- tapply(pgA$vmr, pgA$chem, median, na.rm = TRUE)

pA <- ggplot(pgA, aes(mean, var, colour = chem)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.75, size = 2.3) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("direct-capture" = "#D55E00", "CROP-seq" = "#0072B2"),
    labels = c(sprintf("CROP-seq (median var/mean %.2f)", med["CROP-seq"]),
               sprintf("direct-capture (median var/mean %.1f)", med["direct-capture"]))) +
  labs(title = "Ground-truth ambient dispersion, as CLEANSER reports it (Fig 3B)",
       subtitle = "Per wrong-species guide: variance vs mean of UMIs across host-species cells (incl. zeros)",
       x = "per-guide mean ambient UMIs (log)", y = "per-guide variance (log)", colour = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")
ggsave(file.path(OUT, "barnyard_fig3b_repro.png"), pA, width = 8, height = 5.2, dpi = 130)

# ---- Figure 2: de-doubleting collapses direct-capture onto the Poisson line -----
rows <- lapply(c(1.01, 0.10, 0.01), function(cut) {
  keep <- wrong_frac <= cut
  pg   <- perguide(dc, if (cut > 1) NULL else keep)
  lab  <- if (cut > 1) sprintf("all GEX-pure cells (median vmr %.1f)", median(pg$vmr, na.rm = TRUE))
          else sprintf("drop gRNA wrong-frac > %.2f  (median vmr %.2f)", cut, median(pg$vmr, na.rm = TRUE))
  transform(pg[is.finite(pg$mean) & pg$mean > 0, ], filter = lab)
})
pgB <- do.call(rbind, rows)
pgB$filter <- factor(pgB$filter, levels = vapply(rows, function(d) as.character(d$filter[1]), ""))
pB <- ggplot(pgB, aes(mean, var, colour = filter)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.8, size = 2.4) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("#D55E00", "#E69F00", "#0072B2")) +
  labs(title = "Removing gRNA-channel doublets collapses it onto Poisson",
       subtitle = "Direct-capture, mouse-library guides in human cells; dashed = var = mean (Poisson)",
       x = "per-guide mean ambient UMIs (log)", y = "per-guide variance (log)", colour = "cell filter") +
  theme_bw(base_size = 12) + theme(legend.position = "top", legend.direction = "vertical")
ggsave(file.path(OUT, "barnyard_dedoublet_collapse.png"), pB, width = 8, height = 5.6, dpi = 130)

cat(sprintf("direct-capture median var/mean (raw): %.2f\n", med["direct-capture"]))
cat(sprintf("CROP-seq median var/mean: %.2f\n", med["CROP-seq"]))
cat(sprintf("de-doubleted (wrong_frac<=0.01) median var/mean: %.2f\n",
            median(perguide(dc, wrong_frac <= 0.01)$vmr, na.rm = TRUE)))
cat("wrote figures to", OUT, "\n")
