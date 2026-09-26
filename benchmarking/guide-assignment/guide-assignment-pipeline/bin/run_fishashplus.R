#!/usr/bin/env Rscript

# fishash+ guide assignment.
#
#   Rscript run_fishashplus.R <grna_matrix.rds> <dataset_id>
#
# Reads the same guides x cells gRNA count matrix sceptre does
# (<dataset>/sceptre/grna_matrix.rds), via fishashplus::read_grna_matrix(). The
# .rds carries real dimnames (guide IDs, cell barcodes), which fishash_plus()
# passes through, so assignments are written by NAME -- unlike fishash, whose
# .mtx input has no names and is written by row/column index.
#
# Writes assignments_fishashplus.csv with columns cell_id, grna_id: one row per
# assigned (guide, cell) pair, matching the crispat/cleanser/pertpy format.

suppressPackageStartupMessages({
  library(fishashplus)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: run_fishashplus.R <grna_matrix.rds> <dataset_id>")
}
rds_fp     <- args[1]
dataset_id <- args[2]

# Every fishash_plus() argument is named here so the settings can be read without
# opening the package. All are the PACKAGE DEFAULTS (fishashplus 0.1.0).
Q         <- 0.05   # target false discovery rate (Guo-Sarkar)
REFIT     <- 10     # ambient-field refits after the first pass; stops early once
                    # the assignment stops changing (at most REFIT + 1 iterations)
MIN_COUNT <- 2      # minimum count for an entry to be assigned

cat("Loading:", rds_fp, "\n")
counts <- read_grna_matrix(rds_fp, format = "rds")
cat("  ", nrow(counts), "gRNAs x", ncol(counts), "cells;",
    Matrix::nnzero(counts), "nonzero entries\n")

# Assignments are reported by name, so refuse to run on an input without them
# rather than silently writing blank or made-up IDs.
if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
    anyNA(rownames(counts)) || anyNA(colnames(counts))) {
  stop("input matrix must have guide rownames and cell colnames: ", rds_fp)
}
if (anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
  stop("input matrix has duplicated guide or cell names: ", rds_fp)
}

cat("Running fishash_plus (q =", Q, ", refit =", REFIT, ", min_count =", MIN_COUNT, ")\n")
res <- fishash_plus(
  counts,
  q              = Q,
  refit          = REFIT,
  min_count      = MIN_COUNT,
  return_details = TRUE,
  verbose        = TRUE
)
cat("  ran", res$n_iter, "iterations; converged:", res$converged, "\n")

# The "assigned" matrix is a logical sparse matrix, guides x cells, with the
# input's dimnames. Check it came back on the same grid, then emit the names of
# each TRUE entry.
assigned <- res$assigned
stopifnot(identical(dim(assigned), dim(counts)),
          identical(dimnames(assigned), dimnames(counts)))
trip <- Matrix::summary(as(assigned, "TsparseMatrix"))
if (!is.null(trip$x)) trip <- trip[as.logical(trip$x), , drop = FALSE]

out <- data.frame(cell_id = colnames(assigned)[trip$j],
                  grna_id = rownames(assigned)[trip$i])
out <- out[order(out$cell_id, out$grna_id), , drop = FALSE]

write.csv(out, "assignments_fishashplus.csv", row.names = FALSE)
cat("Wrote assignments_fishashplus.csv:", nrow(out), "assignments over",
    length(unique(out$cell_id)), "cells and",
    length(unique(out$grna_id)), "gRNAs\n")
