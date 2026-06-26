#!/usr/bin/env Rscript
# Same-cohort scoring: my plain ambient test (geomux-core) vs the REAL geomux
# (default lor=10, and lor=0) on the identical exported barnyard cohort, with the
# identical Fishash Table-2 metric. Loads directly from results/barnyard_cohort_export/.
suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
EXP  <- file.path(HERE, "results", "barnyard_cohort_export")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

# Table-2 from a per-cell list of assigned guides
t2_from_pairs <- function(cell_of, guide_of, cells, guides) {
  sp <- setNames(cells$species, cells$barcode); nat <- setNames(guides$native, guides$guide)
  correct <- nat[guide_of] == sp[cell_of]
  byc <- split(correct, factor(cell_of, levels = cells$barcode))
  ncor <- vapply(byc, function(x) sum(x, na.rm=TRUE), numeric(1))
  nwro <- vapply(byc, function(x) sum(!x, na.rm=TRUE), numeric(1))
  ncor <- ncor[cells$barcode]; nwro <- nwro[cells$barcode]; ncor[is.na(ncor)] <- 0; nwro[is.na(nwro)] <- 0
  amb <- if (length(unlist(byc))) mean(!unlist(byc), na.rm=TRUE) else NA
  list(acc = mean(ncor >= 1 & nwro == 0), amb_fdr = amb, n_calls = length(cell_of),
       cells_assigned = sum(ncor + nwro > 0))
}
t2_from_matrix <- function(A, cells, guides) {       # A: guides x cells logical
  A <- as(A, "TsparseMatrix")
  t2_from_pairs(cells$barcode[A@j + 1L], guides$guide[A@i + 1L], cells, guides)
}

rows <- list()
for (s in list.dirs(EXP, recursive = FALSE)) {
  cells <- read.csv(file.path(s, "cells.csv")); guides <- read.csv(file.path(s, "guides.csv"))
  gm <- as(Matrix::readMM(file.path(s, "guide_counts.mtx")), "CsparseMatrix")
  rownames(gm) <- guides$guide; colnames(gm) <- cells$barcode
  sn <- basename(s)
  # my plain ambient test (= geomux core, no lor threshold)
  A <- ambient_test_assign(gm, q = 0.05, model = "hypergeometric", n_iter = 1)$assignment_matrix
  r <- t2_from_matrix(A, cells, guides)
  rows[[length(rows)+1]] <- data.frame(sample = sn, method = "ours (plain BH-hyper)",
    acc = round(r$acc,4), amb_fdr = round(r$amb_fdr,4), cells_assigned = r$cells_assigned)
  # real geomux (default lor=10) and (lor=0)
  for (tag in c("default","nolor")) {
    fp <- file.path(s, paste0("geomux_", tag, "_calls.csv"))
    if (file.exists(fp)) {
      g <- read.csv(fp)
      r <- t2_from_pairs(as.character(g$cell_barcode), as.character(g$guide), cells, guides)
      rows[[length(rows)+1]] <- data.frame(sample = sn, method = paste0("geomux ", tag),
        acc = round(r$acc,4), amb_fdr = round(r$amb_fdr,4), cells_assigned = r$cells_assigned)
    }
  }
}
df <- do.call(rbind, rows)
cat("===== SAME-COHORT: ours vs real geomux (Fishash Table-2 accuracy) =====\n")
print(df, row.names = FALSE)
write.csv(df, file.path(HERE, "results", "same_cohort_geomux.csv"), row.names = FALSE)
