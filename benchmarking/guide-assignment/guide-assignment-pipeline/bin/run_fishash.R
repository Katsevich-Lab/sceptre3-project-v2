#!/usr/bin/env Rscript

# fishash guide assignment.
#
#   Rscript run_fishash.R <grna_matrix.mtx> <dataset_id>
#
# Reads the same Matrix Market file cleanser does: the guides x cells count
# matrix under <dataset>/cleanser/. fishash() takes a bare count matrix, so no
# separate input format is needed. Matrix::writeMM drops dimnames, so guides and
# cells are identified by 1-based row and column index -- the same convention
# assignments_cleanser_*.csv already uses, since it reads this same file.
#
# Writes assignments_fishash.csv with columns cell_id, grna_id: one row per
# assigned (guide, cell) pair, matching the crispat/cleanser/pertpy format.

suppressPackageStartupMessages({
  library(fishash)
  library(Matrix)
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: run_fishash.R <grna_matrix.mtx> <dataset_id>")
}
mtx_fp     <- args[1]
dataset_id <- args[2]

# Every fishash() argument is named here so the settings can be read without
# opening the package.
#
# refit is the iterative Simpson's-paradox correction: the test is rerun up to
# this many times, or until the assignments stop changing. fishash_analysis runs
# refit in {0, 10} (main.nf:391); its plotting helper labels refit=10 as
# "fishash" and does not show refit=0, so 10 is the method as published and 0
# isolates the correction. The realized iteration count is reported below.
REFIT         <- 10      # fishash_analysis bin/run_fishash.R, run at 10
PADJ_CUTOFF   <- 0.05    # package default
EXCLUDE_EMPTY <- TRUE    # package default, set explicitly upstream too

# PACKAGE DEFAULTS, named so they are visible rather than implicit. Neither
# fishash_analysis nor grna-count-modeling varies any of these.
PADJ_METHOD <- "GS"      # Guo & Sarkar 2020 block-dependent FDR correction
MIN_COUNT   <- 2         # minimum UMIs to call a guide present
MIN_FRAC    <- 0         # minimum within-cell count share to call a guide present

cat("Loading:", mtx_fp, "\n")
counts <- as(Matrix::readMM(mtx_fp), "CsparseMatrix")
cat("  ", nrow(counts), "gRNAs x", ncol(counts), "cells;",
    Matrix::nnzero(counts), "nonzero entries\n")

# fishash() requires rownames: fishash.R:259-265 indexes rownames(counts) to build
# the per-cell `assignment` string, and NULL rownames error out in tapply. writeMM
# does not store dimnames, so supply them here, in the same order as the file.
# These names never leave this script -- the CSV below carries row and column
# indices -- but they follow the naming the dataset builders use.
rownames(counts) <- paste0("grna_", seq_len(nrow(counts)))
colnames(counts) <- paste0("CELL_", seq_len(ncol(counts)))

cat("Running fishash (refit =", REFIT, ", padj_cutoff =", PADJ_CUTOFF,
    ", padj_method =", PADJ_METHOD, ", min_count =", MIN_COUNT,
    ", min_frac =", MIN_FRAC, ", exclude_empty =", EXCLUDE_EMPTY, ")\n")
res <- fishash(
  counts,
  refit         = REFIT,
  padj_cutoff   = PADJ_CUTOFF,
  padj_method   = PADJ_METHOD,
  min_count     = MIN_COUNT,
  min_frac      = MIN_FRAC,
  exclude_empty = EXCLUDE_EMPTY
)
cat("  stopped after", metadata(res)$num_iter, "iterations\n")
print(table(colData(res)$demux_type))

# The "assigned" assay is a logical sparse matrix, guides x cells. Go through the
# triplet form to get the (i, j) of each TRUE, then emit j (cell) and i (guide).
# The indices are only meaningful against the input file's row/column order, so
# check that fishash returned the matrix on the same grid it was given.
assigned <- assay(res, "assigned")
stopifnot(identical(dim(assigned), dim(counts)),
          identical(rownames(assigned), rownames(counts)),
          identical(colnames(assigned), colnames(counts)))
trip <- Matrix::summary(as(assigned, "TsparseMatrix"))
if (!is.null(trip$x)) trip <- trip[as.logical(trip$x), , drop = FALSE]

out <- data.frame(cell_id = trip$j, grna_id = trip$i)
out <- out[order(out$cell_id, out$grna_id), , drop = FALSE]

write.csv(out, "assignments_fishash.csv", row.names = FALSE)
cat("Wrote assignments_fishash.csv:", nrow(out), "assignments over",
    length(unique(out$cell_id)), "cells and",
    length(unique(out$grna_id)), "gRNAs\n")
