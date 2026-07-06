#!/usr/bin/env Rscript
# make_mccutcheon.R -- build the benchmarking mccutcheon guide matrix from the RAW import-mccutcheon-2023
# output. QC applied here (project repo, NOT the import repo): keep only guide-bearing cells
# (colSums > 0), 24,244 -> 22,100. Reproduces the committed grna_matrix.rds byte-identically.
# Companion to make_barnyard.R / make_gasperini.R / make_replogle-rd7.R.
suppressMessages(library(Matrix))

base <- .get_config_path("LOCAL_MCCUTCHEON_2023_DATA_DIR")
grna <- as(readRDS(paste0(base, "processed/grna_matrix.rds")), "CsparseMatrix")   # 40 x all cells (raw)
keep <- Matrix::colSums(grna) > 0                                                 # guide-bearing-cell QC
out  <- as(grna[, keep, drop = FALSE], "CsparseMatrix")
saveRDS(out, paste0(base, "grna_matrix.rds"))                                     # registry location
cat(sprintf("mccutcheon: %d -> %d guide-bearing cells\n", ncol(grna), ncol(out)))
