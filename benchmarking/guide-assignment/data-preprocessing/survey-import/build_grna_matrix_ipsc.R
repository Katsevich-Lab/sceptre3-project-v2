#!/usr/bin/env Rscript
# Build sparse gRNA-by-cell dgCMatrix from the HipSci genome-wide CRISPRi
# figshare "Guide-UMI-Counts.csv.gz": DENSE csv, guides as rows (col1 = guide id),
# cells as columns (header = cell barcodes). 6824 guides x 322746 cells.
# Uses data.table::fread (fast C parser), then extracts nonzeros row-by-row.
suppressMessages({library(Matrix); library(data.table)})

dir <- "/Users/ekatsevi/data/external/perturbseq-survey/ipsc_crispri_hipsci_figshare27989294"
inf <- file.path(dir, "GenomeWideScreen_FitnessGenes_Guide-UMI-Counts.csv.gz")

cat("reading with fread...\n")
DT <- fread(cmd = paste("gzcat", shQuote(inf)), header = TRUE, sep = ",",
            colClasses = "integer", showProgress = FALSE)
# first column holds guide ids (its name is "V1" or "")
guide_ids <- as.character(DT[[1]])
cells <- colnames(DT)[-1]
DT[, (1) := NULL]                  # drop guide-id column
cat("guides:", length(guide_ids), " cells:", length(cells), "\n")

ng <- length(guide_ids); nc <- length(cells)
ii <- integer(0); jj <- integer(0); xx <- integer(0)
for (g in seq_len(ng)) {
  v <- DT[[1L]]                    # placeholder; we index rows below
  break
}
# Extract nonzeros row-wise. DT is cells-in-columns; row g = guide g.
# Pull each guide row as a vector (data.table stores by column, so build via
# a single matrix coercion in chunks of cells to limit memory).
chunk <- 20000L
for (start in seq(1L, nc, by = chunk)) {
  end <- min(start + chunk - 1L, nc)
  sub <- as.matrix(DT[, start:end, with = FALSE])   # guides x chunkcells
  idx <- which(sub != 0L, arr.ind = TRUE)
  if (nrow(idx)) {
    ii <- c(ii, idx[, 1])
    jj <- c(jj, idx[, 2] + (start - 1L))
    xx <- c(xx, sub[idx])
  }
}
cat("nnz:", length(xx), "\n")

m <- sparseMatrix(i = ii, j = jj, x = as.numeric(xx),
                  dims = c(ng, nc), dimnames = list(guide_ids, cells))
m <- as(m, "CsparseMatrix"); storage.mode(m@x) <- "double"

cat("gRNA matrix dim (guides x cells):", nrow(m), "x", ncol(m), "\n")
cat("total guide UMIs:", sum(m@x), " nnz:", length(m@x), "\n")
cat("mean guides-per-cell (>0):", round(mean(colSums(m > 0)), 3), "\n")

saveRDS(m, file.path(dir, "grna_matrix.rds"))
write.csv(data.frame(id = guide_ids, name = guide_ids),
          file.path(dir, "guide_features.csv"), row.names = FALSE)
cat("Saved:", file.path(dir, "grna_matrix.rds"), "\n")
