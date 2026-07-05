## Careful AGGREGATE dose-response for the TWO OTHER clean-gap cell types (EndoC T2D + CD4 T-cell),
## the exact companion of scripts/collab_dose_response_aggregate.R (Replogle) so §3's "it holds in
## other cell types" can rest on per-dataset figures instead of the old cross-dataset one.
##
## Methodology MIRRORS the Replogle aggregate script exactly:
##   - power-positive clean-gap guides only: clean gap (gap_lo>=2 from the perguide CSV), well-
##     expressed target (baseline >= 0.5), integration cells (count >= gap_hi), real knockdown
##     (integration mean/baseline < 0.7).
##   - a cell's LOW mode = count <= that guide's gap_lo. Bin the low-mode counts (1,2,3,4,5,6-8),
##     AVERAGE WITHIN EACH GUIDE FIRST, then across guides (+/- SE across guides), annotate the
##     number of contributing guides per bin.
##   - the strong-integration effect is a horizontal REFERENCE BAND (a level), tied to no x.
##   - LINEAR x-axis; the "6-8" bin's point sits at x = 6 (its lowest count).
##
## DATA differences from Replogle (survey triples; loaded exactly as dose_survey() in
## scripts/lowermode_dose_response.R):
##   - grna_matrix_aligned.rds (guides x cells), response_matrix.rds (genes x cells sparse),
##     grna_target_data_frame.csv (grna_id, grna_target = Ensembl id matching response rownames).
##   - library size = colSums(response_matrix); normalized expr = response[target,]/lib*1e4.
##   - EndoC HAS non-targeting guides -> baseline = mean normalized target expr over NT-positive
##     cells (any NT count >=30) with this guide's count == 0.
##   - CD4 has NO non-targeting guides -> complement baseline = cells with count 0 for that guide
##     (dose_survey(..., has_nt=FALSE)); genome-scale, a little slow but fine.
##
## CAP: in both survey datasets the low mode is concentrated at counts 1-2 (guide gaps are narrow:
## gap_lo mostly 2-4, unlike Replogle's up to 17). We therefore keep only bins reaching
## >= MIN_GUIDES_PER_BIN guides (the analogue of Replogle's count-8 cap) -- which for these
## datasets stops at count 2. Every retained low-mode count sits at/below its guide's gap (verified).
##
## Run from nonparametric-thresholds/. CD4 step touches the genome-scale response matrix (~a minute).
suppressMessages({library(Matrix); library(ggplot2); library(patchwork)})
HERE    <- "/Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/nonparametric-thresholds"
OUT_SRC <- file.path(HERE, "results/global_ambient_poisson")
OUT     <- file.path(HERE, "results/collaborator_writeup"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SV      <- path.expand("~/data/external/perturbseq-survey")

MIN_CELLS_PER_GUIDE <- 3     # a guide contributes to a bin only with >= this many low-mode cells there (Replogle value)
MIN_GUIDES_PER_BIN  <- 5     # a bin is shown only if >= this many guides contribute (the cap; smaller than
                             #   Replogle's 20 because these datasets have far fewer clean-gap guides)
MIN_SIG_CELLS       <- 3     # min integration cells for a guide to contribute to the strong-integration band

## Returns the per-bin aggregate + strong-integration level for one survey dataset.
compute_agg <- function(dir, has_nt, pgfile, bins, tag) {
  g  <- as(readRDS(file.path(dir, "grna_matrix_aligned.rds")), "CsparseMatrix")
  r  <- as(readRDS(file.path(dir, "response_matrix.rds")),      "CsparseMatrix")
  td <- read.csv(file.path(dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
  lib  <- as.numeric(Matrix::colSums(r)); rn <- rownames(r); tmap <- setNames(td$grna_target, td$grna_id)

  ## baseline mask: NT-positive cells if the dataset has NT guides, else all cells (v==0 then selects the complement)
  nt <- td$grna_id[td$grna_target == "non-targeting"]; ntr <- match(nt, rownames(g)); ntr <- ntr[!is.na(ntr)]
  ntmask <- if (has_nt && length(ntr)) Matrix::colSums(g[ntr, , drop = FALSE] >= 30) >= 1 else rep(TRUE, ncol(g))

  pg <- read.csv(file.path(OUT_SRC, pgfile), stringsAsFactors = FALSE); pg <- pg[pg$gap_lo >= 2, ]
  pg$target <- tmap[pg$guide]
  pg <- pg[!is.na(pg$target) & pg$target != "non-targeting" & pg$target %in% rn & pg$guide %in% rownames(g), ]

  Gsub <- as.matrix(g[pg$guide, , drop = FALSE])                 # guides are rows; slice once (small)
  utg  <- unique(pg$target); Est <- t(r[utg, , drop = FALSE]); colnames(Est) <- utg

  perbin <- setNames(vector("list", length(bins)), names(bins)); sigvec <- c(); nguide <- 0; membership_ok <- TRUE
  for (i in seq_len(nrow(pg))) {
    v  <- as.numeric(Gsub[i, ]); cp <- as.numeric(Est[, pg$target[i]]) / lib * 1e4
    base <- mean(cp[ntmask & v == 0]); gh <- pg$gap_hi[i]; lo <- pg$gap_lo[i]; st <- which(v >= gh)
    if (!is.finite(base) || base < 0.5 || !length(st) || mean(cp[st]) / base > 0.7) next   # power-positive only
    nguide <- nguide + 1; rel <- cp / base
    if (length(st) >= MIN_SIG_CELLS) sigvec <- c(sigvec, mean(rel[st]))
    for (nm in names(bins)) {                                     # LOW MODE only: count in bin AND <= this guide's gap_lo
      idx <- which(v %in% bins[[nm]] & v <= lo)
      if (length(idx) >= MIN_CELLS_PER_GUIDE) {
        perbin[[nm]] <- c(perbin[[nm]], mean(rel[idx]))
        if (any(v[idx] > lo)) membership_ok <- FALSE              # must never trigger (idx enforces v<=lo)
      }
    }
  }
  agg <- do.call(rbind, lapply(names(bins), function(nm) {
    m <- perbin[[nm]]; if (length(m) < MIN_GUIDES_PER_BIN) return(NULL)
    data.frame(bin = nm, x = min(bins[[nm]]), n_guides = length(m), mean = mean(m), se = sd(m) / sqrt(length(m)))
  }))
  sig_mean <- mean(sigvec); sig_se <- sd(sigvec) / sqrt(length(sigvec))
  cat(sprintf("[%s] clean-gap guides=%d  power-positive=%d  band from %d guides\n",
              tag, nrow(pg), nguide, length(sigvec)))
  cat(sprintf("[%s] per-bin (per-guide-averaged):\n", tag)); print(agg, row.names = FALSE)
  cat(sprintf("[%s] strong-integration level: %.3f +/- %.3f (%d guides)\n", tag, sig_mean, sig_se, length(sigvec)))
  cat(sprintf("[%s] low-mode-membership OK (every binned cell <= gap_lo): %s\n\n", tag, membership_ok))
  list(agg = agg, sig_mean = sig_mean, sig_se = sig_se, nguide = nguide, n_sig_guides = length(sigvec),
       membership_ok = membership_ok)
}

## one panel, styled exactly like the Replogle aggregate figure (linear x, band = a level, n_guides annotated).
## When only counts 1-2 survive the cap (the survey datasets), keep the x-window tight so the two-point
## descent isn't visually stretched; the "N guides" labels go BELOW the points so they never collide with
## the baseline dashed line at count 1.
make_panel <- function(res, title) {
  agg <- res$agg; sig_mean <- res$sig_mean; sig_se <- res$sig_se
  xr <- range(agg$x); ymin <- min(sig_mean - sig_se, min(agg$mean - agg$se), 0)
  xpad <- if (diff(xr) <= 1) 0.7 else 0.4
  ggplot(agg, aes(x, mean)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = sig_mean - sig_se, ymax = sig_mean + sig_se,
             fill = "#D55E00", alpha = 0.15) +
    geom_hline(yintercept = sig_mean, colour = "#D55E00", linewidth = 0.7, linetype = "dashed") +
    annotate("text", x = mean(xr), y = sig_mean, hjust = 0.5, vjust = -0.6, size = 4.0, colour = "#B5490A",
             label = sprintf("strong-integration level (%.2f)", sig_mean)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    annotate("text", x = xr[2] + xpad, y = 1, hjust = 1, vjust = 1.4, size = 3.9, colour = "grey40",
             label = "no-guide baseline") +
    geom_line(colour = "#0072B2", linewidth = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.10, colour = "#0072B2") +
    geom_point(colour = "#0072B2", size = 3) +
    geom_text(aes(label = paste0(n_guides, " guides"), y = mean - se), vjust = 1.9, size = 3.4, colour = "grey35") +
    scale_x_continuous(breaks = agg$x, labels = agg$bin) +
    coord_cartesian(xlim = c(xr[1] - xpad, xr[2] + xpad), ylim = c(ymin, 1.14)) +
    labs(title = title,
         x = "gRNA UMI count in the cell (low mode, below each guide's gap)",
         y = "target expression / no-guide baseline\n(mean over guides ± SE)") +
    theme_bw(base_size = 15)
}

bins_survey <- list(`1` = 1, `2` = 2, `3` = 3, `4` = 4, `5` = 5, `6-8` = 6:8)

cat("################## EndoC (T2D perturb-seq) ##################\n")
en <- compute_agg(file.path(SV, "endoc_t2d_perturbseq_GSE273677/sceptre"), has_nt = TRUE,
                  pgfile = "perguide_endoc_t2d.csv", bins = bins_survey, tag = "EndoC")
cat("################## CD4 T-cell (genome-scale) ##################\n")
tc <- compute_agg(file.path(SV, "tcell_cd4_perturbseq_GSE314342/sceptre"), has_nt = FALSE,
                  pgfile = "perguide_tcell_cd4.csv", bins = bins_survey, tag = "CD4")

## write the per-bin tables next to the figure for provenance
write.csv(transform(en$agg, dataset = "EndoC"), file.path(OUT, "dose_response_aggregate_endoc.csv"), row.names = FALSE)
write.csv(transform(tc$agg, dataset = "CD4"),   file.path(OUT, "dose_response_aggregate_cd4.csv"),   row.names = FALSE)

pL <- make_panel(en, sprintf("EndoC (T2D)  —  %d power-positive clean-gap guides", en$nguide))
pR <- make_panel(tc, sprintf("CD4 T-cell  —  %d power-positive clean-gap guides", tc$nguide))
fig <- pL + pR + plot_annotation(
  title = "Low-mode dose-response holds in other clean-gap cell types",
  subtitle = "Per-guide-averaged target knockdown vs low-mode gRNA count; strong-integration band = full-knockdown level",
  theme = theme(plot.title = element_text(size = 16, face = "bold"),
                plot.subtitle = element_text(size = 12, colour = "grey30")))
ggsave(file.path(OUT, "dose_response_aggregate_other.png"), fig, width = 15.5, height = 5.6, dpi = 130)
cat("wrote", file.path(OUT, "dose_response_aggregate_other.png"), "\n")
cat("wrote", file.path(OUT, "dose_response_aggregate_endoc.csv"), "and _cd4.csv\n")
