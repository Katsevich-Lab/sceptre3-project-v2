#!/usr/bin/env Rscript
# Export gRNA matrices (Matrix Market) for the smaller datasets so CLEANSER can be run on them in
# batch for the method-comparison writeup. Skips the huge-guide datasets (CLEANSER MCMC infeasible).
suppressPackageStartupMessages(library(Matrix))
source("scripts/datasets.R")
BATCH <- "results/ambient_ceiling/cleanser_batch"
dir.create(BATCH, showWarnings=FALSE, recursive=TRUE)

# smallest guide-count first, so we get dataset coverage fastest. (Skips gasperini/replogle/ipsc/cd4:
# too many guides for CLEANSER MCMC.)
paths <- dataset_paths()[c(
  "invivo_cortex","mccutcheon","cd8_tcell","erythroid_multiome","gastric_organoid",
  "dctap_k562_lowmoi","dctap_k562_highmoi","barnyard_mch2_72hr","barnyard_mch2_0hr",
  "barnyard_lrb100_0hr","barnyard_lrb100_72hr","endoc_t2d","a549")]

for (nm in names(paths)) {
  m <- as(readRDS(paths[[nm]]),"CsparseMatrix")
  gn <- rownames(m); if (is.null(gn)) gn <- as.character(seq_len(nrow(m)))
  cn <- colnames(m); if (is.null(cn)) cn <- paste0("cell",seq_len(ncol(m)))
  d <- file.path(BATCH,nm); dir.create(d,showWarnings=FALSE)
  writeMM(m, file.path(d,"m.mtx")); writeLines(gn,file.path(d,"guides.txt")); writeLines(cn,file.path(d,"cells.txt"))
  cat(sprintf("%-20s %d guides x %d cells\n", nm, nrow(m), ncol(m)))
}
writeLines(names(paths), file.path(BATCH,"order.txt"))
