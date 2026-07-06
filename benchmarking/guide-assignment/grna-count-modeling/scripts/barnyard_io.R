#!/usr/bin/env Rscript
# =============================================================================
# barnyard_io.R  --  canonical loader for the 4 Liu-2025 species-mix samples,
# using the ORIGINAL AUTHORS' QC choices verbatim (from fishash_analysis
# 82_plot_barnyard_results.R, reproduced in external/repro_work/score_table2.R).
# We do NOT invent our own purity sweep -- the authors' QC defines the cohort.
#
# A QC-passing cell is species-PURE (frac_homo < .1 or > .9) with standard GEX
# QC.  For a pure human cell the mouse guides (nt_101..200) are provably ambient;
# for a pure mouse cell the human guides (nt_1..100) are.  That wrong-species
# submatrix is our real-ambient ground truth (whatever its dispersion turns out
# to be -- CROP-seq near-Poisson, direct-capture heavy-tailed).
# =============================================================================
suppressPackageStartupMessages({ library(Matrix) })

BARN_SAMPLES <- c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")

# Authors' QC mask (verbatim from score_table2.R lines 39-47), with the species
# PURITY cut exposed as a single knob.  purity = 0.9 reproduces the authors' QC
# exactly; purity = 0.99 is the strict (doublet-suppressed) setting -- ONLY the
# purity threshold changes; every other cap (mito, features, sum_gex) is held
# identical, so 0.9 vs 0.99 is a clean one-parameter sensitivity.
.barnyard_qc <- function(meta, purity = 0.9) {
  sum_gex   <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_homo <- meta$homo_sum_gex / sum_gex
  qc <- (meta$mito_sum / sum_gex < .15) & (meta$features_gex <= 6000) & (sum_gex <= 20000) &
        (meta$features_gex >= 1500) & (sum_gex >= 3500) &
        (frac_homo < (1 - purity) | frac_homo > purity)
  qc[is.na(qc)] <- FALSE
  list(qc = qc, frac_homo = frac_homo)
}

# Returns, for one sample (after authors' QC):
#   counts        guides x QC-cells sparse count matrix (rownames nt_1..nt_200)
#   guide_homo    logical, TRUE for human (homo) guides
#   cell_homo     logical, TRUE for pure human cells (frac_homo > .9)
#   ambient_mask  guides x cells logical: TRUE where guide species != cell species
#                 (provably-ambient wrong-species entries)
#   chemistry     "CROP-seq" or "direct-capture"
load_barnyard <- function(sample, purity = 0.9, repro_dir = NULL) {
  if (is.null(repro_dir)) {
    h <- if (file.exists(file.path(getwd(), "external", "repro_work"))) getwd()
         else file.path(Sys.getenv("SIM_HERE", getwd()))
    repro_dir <- file.path(h, "external", "repro_work")
  }
  gm     <- as(readMM(file.path(repro_dir, paste0(sample, "_grna_counts.mtx"))), "CsparseMatrix")
  meta   <- read.csv(file.path(repro_dir, paste0(sample, "_meta.csv")))
  guides <- read.csv(file.path(repro_dir, paste0(sample, "_guides.csv")))
  stopifnot(nrow(meta) == ncol(gm), nrow(guides) == nrow(gm))
  q <- .barnyard_qc(meta, purity = purity)
  counts <- gm[, q$qc, drop = FALSE]
  rownames(counts) <- guides$guide
  colnames(counts) <- paste0(sample, "_", which(q$qc))
  guide_homo <- guides$guide_type == "homo_guide"
  cell_homo  <- q$frac_homo[q$qc] > purity
  amb_mask <- as(Matrix(outer(guide_homo, cell_homo, function(g, c) g != c), sparse = TRUE), "lgCMatrix")
  dimnames(amb_mask) <- dimnames(counts)
  list(counts = counts, guide_homo = guide_homo, cell_homo = cell_homo, purity = purity,
       ambient_mask = amb_mask, chemistry = if (grepl("Cropseq", sample)) "CROP-seq" else "direct-capture",
       sample = sample)
}

# GEX-QC-only loader (mito/features/depth caps, NO species-purity gate).  Callers
# apply the purity cut per-analysis via `fh` (e.g. host = fh > 0.90), and often
# need BOTH species, so this is deliberately looser than load_barnyard().  Used by
# the wrong-species VMR / marginal-fit collaborator figures.  Returns:
#   counts     guides x GEX-QC-cells sparse matrix (rownames = guides, colnames sample_<idx>)
#   guide_homo logical, TRUE for human (homo) guides
#   fh         frac_homo of the kept cells
#   meta       the kept cells' metadata rows
load_barnyard_basic <- function(sample, repro_dir = NULL) {
  if (is.null(repro_dir)) {
    h <- if (file.exists(file.path(getwd(), "external", "repro_work"))) getwd()
         else file.path(Sys.getenv("SIM_HERE", getwd()))
    repro_dir <- file.path(h, "external", "repro_work")
  }
  gm     <- as(readMM(file.path(repro_dir, paste0(sample, "_grna_counts.mtx"))), "CsparseMatrix")
  meta   <- read.csv(file.path(repro_dir, paste0(sample, "_meta.csv")))
  guides <- read.csv(file.path(repro_dir, paste0(sample, "_guides.csv")))
  sg <- meta$homo_sum_gex + meta$mus_sum_gex; fh <- meta$homo_sum_gex / sg
  qc <- (meta$mito_sum/sg < .15) & (meta$features_gex <= 6000) & (sg <= 20000) &
        (meta$features_gex >= 1500) & (sg >= 3500); qc[is.na(qc)] <- FALSE
  counts <- gm[, qc, drop = FALSE]; rownames(counts) <- guides$guide
  colnames(counts) <- paste0(sample, "_", which(qc))
  list(counts = counts, guide_homo = guides$guide_type == "homo_guide",
       fh = fh[qc], meta = meta[qc, ])
}
