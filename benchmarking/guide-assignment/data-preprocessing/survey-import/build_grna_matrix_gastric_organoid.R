#!/usr/bin/env Rscript
# Build sparse gRNA-by-cell dgCMatrix from the CROP-seq cell_identities.csv.gz
# (GSE280506). Each row = one cell with its ASSIGNED guide and that guide's
# per-cell UMI_count. CROP-seq low-MOI: one guide per assigned cell. Cells with
# guide_identity == "*" are unassigned (no guide) and are dropped. The result is
# a guides x cells integer UMI matrix (one nonzero per assigned cell).
suppressMessages(library(Matrix))

dir <- "/Users/ekatsevi/data/external/perturbseq-survey/gastric_organoid_cropseq_GSE280506"
inf <- file.path(dir, "GSE280506_cell_identities.csv.gz")

df <- read.csv(gzfile(inf), stringsAsFactors = FALSE)
cat("rows:", nrow(df), " cols:", paste(colnames(df), collapse=","), "\n")

# keep assigned cells only
keep <- df$guide_identity != "*" & !is.na(df$UMI_count)
df <- df[keep, ]
cat("assigned cells:", nrow(df), "\n")

guides <- sort(unique(df$guide_identity))
cells  <- df$cell_barcode
gi <- match(df$guide_identity, guides)
ci <- seq_len(nrow(df))

m <- sparseMatrix(i = gi, j = ci, x = as.numeric(df$UMI_count),
                  dims = c(length(guides), nrow(df)),
                  dimnames = list(guides, cells))
m <- as(m, "CsparseMatrix")
storage.mode(m@x) <- "double"

cat("gRNA matrix dim (guides x cells):", nrow(m), "x", ncol(m), "\n")
cat("total guide UMIs:", sum(m@x), " nnz:", length(m@x), "\n")
cat("mean guides-per-cell (>0):", round(mean(colSums(m > 0)), 3), "\n")

saveRDS(m, file.path(dir, "grna_matrix.rds"))
write.csv(data.frame(id = guides, name = guides),
          file.path(dir, "guide_features.csv"), row.names = FALSE)
cat("Saved:", file.path(dir, "grna_matrix.rds"), "\n")
