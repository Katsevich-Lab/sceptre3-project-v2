#!/usr/bin/env Rscript
# =============================================================================
# barnyard_five_cells.R  --  the CLEANER, non-circular version of the §1
# "the 0.99 residual is spikes, not soup" story for the collaborator writeup.
#
# The existing panel (barnyard_vmr_threshold.png, panel 3) deletes EACH guide's
# own single largest cell -- slightly circular ("of course removing a guide's max
# shrinks its variance"). This script instead shows the residual over-dispersion
# is driven by the SAME small shared pool of cells across guides, then removes
# just THOSE cells GLOBALLY (not per-guide) and shows the whole var-vs-mean plot
# collapses onto the Poisson line.
#
# Setting: mix0hr_DirectCapture, 0.99 GEX purity gate (doublet-suppressed),
# wrong-species (mouse, nt_101..nt_200) guides across host (human, frac_homo>0.99)
# cells. Reuses the load_basic / perguide-var conventions of barnyard_ambient_figures.R
# exactly so numbers stay consistent with the current §1 figures.
#
# Outputs into results/collaborator_writeup/:
#   barnyard_top3_marginals.png   (step 3: marginal count histograms of top-3 guides)
#   barnyard_five_cells.csv       (step 4: the shared offender cells documented)
#   barnyard_vmr_five_cells.png   (step 5: var-vs-mean with offenders removed GLOBALLY)
# Run from grna-count-modeling/.
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2) })
REPRO <- "external/repro_work"; OUT <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PAL <- c(blue = "#0072B2", orange = "#D55E00", gold = "#E69F00")

# --- loader: shared GEX-QC-only barnyard loader (same basic-QC cohort as the §1 figures) ---
source("scripts/barnyard_io.R")   # load_barnyard_basic()
load_basic <- function(sample) load_barnyard_basic(sample, REPRO)

# per wrong-species (mouse) guide: mean/var/vmr of UMIs across a given host-cell
# set (incl. zeros); SAME estimator as perguide_vmr / pg_vmr.
perguide_vmr <- function(M) {
  nC <- ncol(M); s1 <- Matrix::rowSums(M); s2 <- Matrix::rowSums(M * M)
  mu <- s1/nC; v <- (s2 - nC*mu^2)/(nC - 1)
  data.frame(guide = rownames(M), mean = mu, var = v,
             vmr = ifelse(mu > 0, v/mu, NA), row.names = NULL)
}

dc <- load_basic("mix0hr_DirectCapture")

## ============================================================================
## STEP 1: 0.99 gate, wrong-species (mouse) guides x host (human) cells
## ============================================================================
host   <- which(dc$fh > 0.99)                       # host = pure human cells
mg     <- which(!dc$guide_homo)                     # mouse (wrong-species) guides
W      <- dc$counts[mg, host, drop = FALSE]         # mouse-guide x host-cell submatrix
Wm     <- as.matrix(W)
cat(sprintf("STEP 1: %d host cells (frac_homo>0.99), %d wrong-species guides\n",
            length(host), length(mg)))

## ============================================================================
## STEP 2: per-guide variance across host cells; take the 3 highest-variance guides
## ============================================================================
pg <- perguide_vmr(W)
pg <- pg[order(-pg$var), ]
top3 <- head(pg, 3)
cat("\nSTEP 2: top-3 highest-variance wrong-species guides @0.99 gate:\n")
print(transform(top3, mean = round(mean, 4), var = round(var, 2), vmr = round(vmr, 1)),
      row.names = FALSE)

## ============================================================================
## STEP 3: marginal count histograms of the top-3 guides (host cells), highlight
## the extreme cell(s) driving each guide's variance
## ============================================================================
top3_ids <- top3$guide
hist_df <- do.call(rbind, lapply(top3_ids, function(g) {
  x <- Wm[g, ]
  data.frame(guide = g, count = x)
}))
# label each guide facet with its (mean,var,vmr) and its max count
lab_map <- setNames(sprintf("%s  (max %d UMIs; var/mean %.0f)",
                            top3$guide, apply(Wm[top3_ids, , drop = FALSE], 1, max),
                            top3$vmr), top3$guide)
hist_df$facet <- factor(lab_map[hist_df$guide], levels = lab_map[top3_ids])
# only the nonzero counts are interesting; zeros dominate and would flatten the plot.
# On a log x-axis geom_histogram bars render badly, so precompute a per-(guide,count)
# frequency table (count value k -> number of host cells with that count) and plot it
# as a lollipop (geom_segment to a floor + geom_point) on a log-log scale. All nonzero
# counts are >=1 and all cell-frequencies are >=1, so log10 is well-defined.
hz <- hist_df[hist_df$count > 0, ]
freq_df <- as.data.frame(table(guide = hz$guide, count = hz$count),
                         stringsAsFactors = FALSE)
freq_df <- freq_df[freq_df$Freq > 0, ]
freq_df$count <- as.integer(freq_df$count)
freq_df$facet <- factor(lab_map[freq_df$guide], levels = lab_map[top3_ids])
YFLOOR <- 0.7  # floor for the lollipop stems on the log y-axis (below the min freq of 1)
p_hist <- ggplot(freq_df, aes(count, Freq)) +
  geom_segment(aes(xend = count, yend = YFLOOR), colour = PAL["blue"], linewidth = 0.4) +
  geom_point(colour = PAL["blue"], size = 1.8) +
  facet_wrap(~facet, nrow = 1, scales = "free_x") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "wrong-species UMI count in a host cell (log, zeros omitted)", y = "number of host cells (log)") +
  theme_bw(base_size = 16)
ggsave(file.path(OUT, "barnyard_top3_marginals.png"), p_hist, width = 10, height = 4.2, dpi = 130)

## ============================================================================
## STEP 4: verify the "same shared pool of cells" claim honestly.
## ============================================================================
# For each guide, which host cell carries its top count? Rank guides by variance,
# look at the driving cell of the top-K guides, and see if they collapse to a
# small shared set.
cell_ids <- colnames(W)
top_cell_of_guide <- function(g) { x <- Wm[g, ]; i <- which.max(x); list(cell = cell_ids[i], val = x[i]) }

ord_guides <- pg$guide                                   # guides ordered by var desc
K <- 10
drivers <- data.frame(rank = seq_len(K), guide = ord_guides[seq_len(K)],
                      driver_cell = NA_character_, driver_count = NA_integer_)
for (k in seq_len(K)) { tc <- top_cell_of_guide(ord_guides[k])
  drivers$driver_cell[k] <- tc$cell; drivers$driver_count[k] <- tc$val }
cat("\nSTEP 4: driver cell of the top-10 highest-variance guides:\n")
print(drivers, row.names = FALSE)
cat("\n  distinct driver cells among top-3 guides: ",
    length(unique(drivers$driver_cell[1:3])), "\n")
cat("  distinct driver cells among top-10 guides:",
    length(unique(drivers$driver_cell[1:10])), "\n")

# Candidate shared offender pool: cells that carry a loud wrong-species spike
# (max wrong-species count >= 10) in ANY guide -- matches barnyard_offenders.csv's
# threshold. These are the cells to consider removing globally.
cell_max <- apply(Wm, 2, max)                            # per host cell: loudest wrong-species UMI
off_cells <- cell_ids[cell_max >= 10]
cat(sprintf("\n  host cells with a wrong-species spike >=10 UMIs: %d\n", length(off_cells)))

# Cross-reference barnyard_offenders.csv (the documented 5 residual offenders).
off_csv <- file.path(OUT, "barnyard_offenders.csv")
ref <- if (file.exists(off_csv)) read.csv(off_csv, stringsAsFactors = FALSE) else NULL

# Per offender cell: which wrong-species guides does it carry, and how loudly?
# (Are the 5 cells one-spike-per-guide, or does each carry multiple guides?)
build_cell_record <- function(cid) {
  ci <- match(cid, cell_ids); x <- Wm[, ci]
  carried <- which(x > 0); ord <- carried[order(-x[carried])]
  top_g <- rownames(Wm)[ord[1]]
  meta_i <- host[ci]
  data.frame(
    cell = cid,
    top_mouse_guide = top_g,
    top_wrong_species_count = x[ord[1]],
    n_mouse_guides_carried = length(carried),
    n_mouse_guides_ge5 = sum(x >= 5),
    frac_homo_GEX = round(dc$fh[meta_i], 4),
    mouse_mRNA_UMIs = dc$meta$mus_sum_gex[meta_i],
    human_mRNA_UMIs = dc$meta$homo_sum_gex[meta_i],
    total_gRNA = Matrix::colSums(dc$counts)[meta_i],
    stringsAsFactors = FALSE)
}
five <- do.call(rbind, lapply(off_cells, build_cell_record))
five <- five[order(-five$top_wrong_species_count), ]

# How much of the TOTAL 0.99-gate variance excess is attributable to these cells?
# Total excess = sum over guides of (var - mean) [Poisson has var=mean; excess>0 is over-disp].
# Recompute per-guide var after removing the offender cells globally, and compare.
excess_of <- function(M) { p <- perguide_vmr(M); e <- p$var - p$mean; sum(e[is.finite(e)]) }
keep_after <- setdiff(cell_ids, off_cells)
W_after <- W[, keep_after, drop = FALSE]
exc_before <- excess_of(W);  exc_after <- excess_of(W_after)
frac_explained <- 1 - exc_after / exc_before
cat(sprintf("\n  per-guide summed variance-excess (var-mean): before %.1f, after removing %d cells %.1f (%.1f%% removed)\n",
            exc_before, length(off_cells), exc_after, 100*frac_explained))

# also report pooled var/mean before/after (the headline "collapse" number)
pooled_vmr <- function(M) { y <- as.numeric(as.matrix(M)); var(y)/mean(y) }
cat(sprintf("  pooled var/mean: before %.2f, after %.2f\n",
            pooled_vmr(W), pooled_vmr(W_after)))

write.csv(five, file.path(OUT, "barnyard_five_cells.csv"), row.names = FALSE)
cat("\n  wrote", file.path(OUT, "barnyard_five_cells.csv"), "\n")
print(five, row.names = FALSE)

## ============================================================================
## STEP 5: var-vs-mean scatter with the offender cells removed GLOBALLY, as a
## multi-panel comparison in the barnyard_vmr_threshold.png style.
## ============================================================================
NOFF <- length(off_cells)
mk_panel <- function(M, nm) transform(perguide_vmr(M), panel = nm)
P <- list(
  `GEX gate 0.99`                = mk_panel(W, "GEX gate 0.99"),
  `minus shared offender cells`  = mk_panel(W_after,
                                    sprintf("GEX gate 0.99, minus the %d shared offender cells (global)", NOFF)))
lev <- c("GEX gate 0.99",
         sprintf("GEX gate 0.99, minus the %d shared offender cells (global)", NOFF))
pgV <- do.call(rbind, P)
pgV <- pgV[is.finite(pgV$mean) & pgV$mean > 0 & is.finite(pgV$var) & pgV$var > 0, ]
medV <- tapply(pgV$vmr, pgV$panel, median, na.rm = TRUE)
pgV$panel <- factor(pgV$panel, levels = lev,
  labels = sprintf("%s\n(median var/mean %.1f)", lev, medV[lev]))
p_vmr <- ggplot(pgV, aes(mean, var)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6, size = 1.5, colour = PAL["blue"]) +
  facet_wrap(~panel, nrow = 1) + scale_x_log10() + scale_y_log10() +
  labs(x = "per-guide mean ambient UMIs (log)", y = "per-guide variance (log)") +
  theme_bw(base_size = 16) + theme(strip.text = element_text(size = 12))
ggsave(file.path(OUT, "barnyard_vmr_five_cells.png"), p_vmr, width = 10, height = 4.6, dpi = 130)

## ============================================================================
## STEP 6: histogram of the per-guide var/mean ratios that REMAIN after the
## 0.99 gate + global removal of the five shared offender cells. Summary
## companion to the var-vs-mean scatter: once the debris cells are gone, the
## distribution of per-guide var/mean piles up tightly at ~1 (Poisson) --
## there is no residual over-dispersed soup.
## ============================================================================
vmr_after <- perguide_vmr(W_after)$vmr
vmr_after <- vmr_after[is.finite(vmr_after) & vmr_after > 0]     # estimable guides (mean>0)
med_after <- median(vmr_after)
frac_le12 <- mean(vmr_after <= 1.2)
# pre-removal 0.99-gate median, for the light contrast annotation
vmr_before <- perguide_vmr(W)$vmr
vmr_before <- vmr_before[is.finite(vmr_before) & vmr_before > 0]
med_before <- median(vmr_before)
cat(sprintf("\nSTEP 6: residual per-guide var/mean after 0.99 gate + global 5-cell removal:\n"))
cat(sprintf("  estimable guides: %d; median var/mean = %.2f; fraction <= 1.2 = %.1f%%\n",
            length(vmr_after), med_after, 100*frac_le12))
cat(sprintf("  (0.99-gate median before removal = %.2f)\n", med_before))

hist_after <- data.frame(vmr = vmr_after)
p_vmr_hist <- ggplot(hist_after, aes(vmr)) +
  geom_histogram(binwidth = 0.05, fill = PAL["blue"], colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = med_after, colour = PAL["orange"], linewidth = 0.7) +
  annotate("text", x = Inf, y = Inf, vjust = 1.6, hjust = 1.05,
           label = sprintf("median var/mean = %.2f  (Poisson = 1)", med_after),
           colour = PAL["orange"], size = 5) +
  labs(x = "per-guide variance / mean", y = "number of wrong-species guides") +
  theme_bw(base_size = 16)
ggsave(file.path(OUT, "barnyard_vmr_histogram.png"), p_vmr_hist, width = 9, height = 5, dpi = 130)

cat("\n=== SUMMARY ===\n")
cat(sprintf("top-3 guides: %s\n", paste(top3$guide, collapse = ", ")))
cat(sprintf("shared offender cells (spike>=10): %d\n", NOFF))
cat(sprintf("variance-excess explained by them: %.1f%%\n", 100*frac_explained))
cat(sprintf("median var/mean: 0.99 gate %.1f -> minus %d cells %.1f\n",
            medV[lev[1]], NOFF, medV[lev[2]]))
cat("\nwrote figures + barnyard_five_cells.csv to", OUT, "\n")
