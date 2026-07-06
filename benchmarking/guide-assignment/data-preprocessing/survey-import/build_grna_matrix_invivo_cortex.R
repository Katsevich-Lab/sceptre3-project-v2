#!/usr/bin/env Rscript
# Build gRNA-by-cell dgCMatrix from GSE249416_Perturb_sg.qs.gz (in vivo cortex
# Perturb-seq, Jin et al. Cell 2024). Legacy `qs` Seurat object; the "Crispr"
# assay holds raw per-cell guide UMI counts. REQUIRES R 4.4 (legacy qs won't
# compile on R>=4.6). Run with the 4.4 Rscript binary.
suppressMessages({library(qs); library(Matrix); library(SeuratObject)})

dir <- "/Users/ekatsevi/data/external/perturbseq-survey/invivo_cortex_perturbseq_GSE249416"
gz  <- file.path(dir, "GSE249416_Perturb_sg.qs.gz")
qsf <- file.path(dir, "GSE249416_Perturb_sg.qs")

if (!file.exists(qsf)) {
  cat("gunzip-ing...\n")
  R.utils::gunzip(gz, destname = qsf, remove = FALSE, overwrite = TRUE)
}

cat("qreading Seurat object...\n")
obj <- qs::qread(qsf)
cat("assays:", paste(names(obj@assays), collapse=", "), "\n")

a <- if ("Crispr" %in% names(obj@assays)) "Crispr" else
     grep("crispr|guide|grna|sgrna", names(obj@assays), ignore.case=TRUE, value=TRUE)[1]
stopifnot(!is.na(a))
cat("using assay:", a, "\n")

m <- tryCatch(SeuratObject::GetAssayData(obj, assay = a, layer = "counts"),
              error = function(e) SeuratObject::GetAssayData(obj, assay = a, slot = "counts"))
m <- as(m, "CsparseMatrix"); storage.mode(m@x) <- "double"

cat("gRNA matrix dim (guides x cells):", nrow(m), "x", ncol(m), "\n")
cat("total guide UMIs:", sum(m@x), " nnz:", length(m@x), "\n")
cat("mean guides-per-cell (>0):", round(mean(colSums(m > 0)), 3), "\n")

saveRDS(m, file.path(dir, "grna_matrix.rds"))
write.csv(data.frame(id = rownames(m), name = rownames(m)),
          file.path(dir, "guide_features.csv"), row.names = FALSE)
cat("Saved:", file.path(dir, "grna_matrix.rds"), "\n")
