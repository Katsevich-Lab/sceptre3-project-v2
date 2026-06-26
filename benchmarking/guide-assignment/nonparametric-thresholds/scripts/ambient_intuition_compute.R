#!/usr/bin/env Rscript
# ============================================================================
# Hands-on look at AMBIENT DEPTH and AMBIENT COMPOSITION on a real CROP-seq
# dataset (fishash's crispat_schraivogel: 86 guides x ~22k cells).
#
# Estimate, for every cell, its background depth kappa_c, and for every guide,
# its background composition gamma_g, three ways:
#   naive     : kappa = library size (colSums) ; gamma = global share (rowSums/N)
#   CLEANSER  : kappa = sum of sub-threshold (<=2) gRNA UMIs ; gamma from same
#   fishash   : kappa = colSums(noise) ; gamma = rowSums(noise)/sum, where noise
#               is fishash's rank-1 Poisson completion of the assigned-masked matrix
# Saves per-cell + per-guide tables and a couple of highlighted slices for plots.
# ============================================================================
suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
  library(sparseMatrixStats)
})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)))), error = function(e) ".")
HERE <- dirname(HERE)
OUT  <- file.path(HERE, "results", "ambient_intuition"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

data(crispat_schraivogel)
counts <- as(crispat_schraivogel, "CsparseMatrix")    # guides x cells
G <- nrow(counts); C <- ncol(counts)
cat(sprintf("counts: %d guides x %d cells, %d UMIs\n", G, C, sum(counts)))

cat("running fishash (refit=10, default GS)...\n")
res      <- fishash(counts)
assigned <- assay(res, "assigned")                    # logical guides x cells
noise    <- impute_masked_counts(counts, assigned)    # fishash's N_hat_noise
cat(sprintf("fishash assigned %d (g,c) pairs; %.1f%% of cells get >=1 guide\n",
            sum(assigned), 100*mean(colSums(assigned) > 0)))

THR <- 2                                               # CLEANSER sub-threshold
sub <- counts * (counts <= THR)                        # keep only counts <= 2

# ---- per-cell depth -------------------------------------------------------
libsize     <- colSums(counts)
depth_clean <- colSums(sub)
depth_fish  <- colSums(noise)
topcount    <- colMaxs(counts)                         # biggest single guide in the cell
n_assigned  <- colSums(assigned)
signal_umi  <- colSums(counts * assigned)              # UMIs sitting in assigned guides
cells <- data.frame(
  cell = colnames(counts), libsize = libsize,
  depth_clean = depth_clean, depth_fish = depth_fish,
  topcount = topcount, top_frac = topcount / pmax(libsize, 1),
  n_assigned = n_assigned, signal_umi = signal_umi,
  ambient_umi = libsize - signal_umi)
write.csv(cells, file.path(OUT, "per_cell.csv"), row.names = FALSE)

# ---- per-guide composition ------------------------------------------------
gshare_naive <- rowSums(counts) / sum(counts)
gc_clean     <- rowSums(sub);   gcomp_clean <- gc_clean / sum(gc_clean)
gc_fish      <- rowSums(noise); gcomp_fish  <- gc_fish / sum(gc_fish)
guides <- data.frame(
  guide = rownames(counts),
  total = rowSums(counts), n_assigned = rowSums(assigned),
  share_naive = gshare_naive, comp_clean = gcomp_clean, comp_fish = gcomp_fish,
  # how much the naive global share OVER-states the true ambient rate:
  inflation = gshare_naive / pmax(gcomp_fish, .Machine$double.eps))
write.csv(guides, file.path(OUT, "per_guide.csv"), row.names = FALSE)

# ---- pick highlights ------------------------------------------------------
# a signal-dominated cell (biggest library dominated by one guide) vs an ambient-only cell
hi_sig <- which.max(ifelse(cells$top_frac > 0.8, cells$libsize, 0))
amb_only <- which(cells$n_assigned == 0 & cells$libsize > quantile(libsize, .6))
amb_only <- amb_only[which.min(abs(libsize[amb_only] - median(libsize)))]
# the most signal-inflated guide (global share >> ambient rate) and a near-honest one
g_inflated <- which.max(guides$inflation * (guides$n_assigned > 50))
g_honest   <- order(abs(log(guides$inflation)))[1]

hl <- list(cell_signal = colnames(counts)[hi_sig],
           cell_ambient = colnames(counts)[amb_only],
           guide_inflated = rownames(counts)[g_inflated],
           guide_honest   = rownames(counts)[g_honest])
saveRDS(hl, file.path(OUT, "highlights.rds"))
cat("\nHighlights:\n"); str(hl)

# full guide profile for the two highlighted cells (for the per-cell floor plot)
prof <- data.frame(
  guide = rownames(counts),
  count_signalcell  = as.numeric(counts[, hi_sig]),
  amb_signalcell    = gcomp_fish * depth_fish[hi_sig],
  assigned_signalcell = as.logical(assigned[, hi_sig]),
  count_ambientcell = as.numeric(counts[, amb_only]),
  amb_ambientcell   = gcomp_fish * depth_fish[amb_only])
write.csv(prof, file.path(OUT, "cell_profiles.csv"), row.names = FALSE)

# one guide's counts across all cells vs depth (for the per-guide scatter)
gi <- g_inflated
gscat <- data.frame(
  cell = colnames(counts),
  depth_fish = depth_fish, libsize = libsize,
  count = as.numeric(counts[gi, ]),
  assigned = as.logical(assigned[gi, ]))
write.csv(gscat, file.path(OUT, "guide_scatter.csv"), row.names = FALSE)
cat(sprintf("\nhighlighted guide for scatter: %s (assigned in %d cells)\n",
            rownames(counts)[gi], guides$n_assigned[gi]))
cat(sprintf("  its naive global share = %.4f ; fishash ambient rate = %.4f ; inflation = %.1fx\n",
            guides$share_naive[gi], guides$comp_fish[gi], guides$inflation[gi]))
cat("wrote tables to", OUT, "\n")
