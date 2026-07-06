#!/usr/bin/env Rscript
# Parse a 10x mtx triplet (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz),
# extract CRISPR Guide Capture features, save sparse gRNA-by-cell dgCMatrix.
#
# Usage: Rscript parse_10x_mtx_guides.R <mtx_dir> <out_dir> [feature_type_regex]
suppressMessages({library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
mtx_dir <- args[1]
out_dir <- args[2]
ft_regex <- if (length(args) >= 3) args[3] else "CRISPR Guide Capture|Guide Capture|CRISPR|gRNA|guide"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

f <- function(...) {
  cands <- file.path(mtx_dir, c(...))
  cands[file.exists(cands)][1]
}
mtx_f  <- f("matrix.mtx.gz", "matrix.mtx")
feat_f <- f("features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv")
bc_f   <- f("barcodes.tsv.gz", "barcodes.tsv")
stopifnot(!is.na(mtx_f), !is.na(feat_f), !is.na(bc_f))

m <- readMM(gzfile(mtx_f))               # features x cells
feat <- read.delim(gzfile(feat_f), header = FALSE, stringsAsFactors = FALSE)
bc   <- readLines(gzfile(bc_f))
cat("matrix dim:", dim(m), " features:", nrow(feat), " barcodes:", length(bc), "\n")

# feature_type usually col 3
ftype <- if (ncol(feat) >= 3) feat[[3]] else rep("Unknown", nrow(feat))
cat("feature_type table:\n"); print(table(ftype))

is_guide <- grepl(ft_regex, ftype)
cat("n guide features matched:", sum(is_guide), "\n")
if (sum(is_guide) == 0) stop("No guide features matched: ", ft_regex)

m <- as(m, "CsparseMatrix")
rownames(m) <- make.unique(as.character(feat[[1]]))
colnames(m) <- bc
grna <- m[is_guide, , drop = FALSE]
# prefer guide name (col2) if present
gid <- if (ncol(feat) >= 2) feat[[2]][is_guide] else feat[[1]][is_guide]
rownames(grna) <- make.unique(as.character(gid))
grna <- as(grna, "CsparseMatrix")
storage.mode(grna@x) <- "double"

cat("gRNA matrix dim (guides x cells):", nrow(grna), "x", ncol(grna), "\n")
cat("total guide UMIs:", sum(grna@x), " nnz:", length(grna@x), "\n")
cat("mean guides-per-cell (>0):", round(mean(colSums(grna > 0)), 3), "\n")

saveRDS(grna, file.path(out_dir, "grna_matrix.rds"))
write.csv(data.frame(id = feat[[1]][is_guide],
                     name = if (ncol(feat) >= 2) feat[[2]][is_guide] else feat[[1]][is_guide],
                     feature_type = ftype[is_guide]),
          file.path(out_dir, "guide_features.csv"), row.names = FALSE)
cat("Saved:", file.path(out_dir, "grna_matrix.rds"), "\n")
