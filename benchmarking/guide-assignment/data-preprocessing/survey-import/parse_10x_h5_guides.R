#!/usr/bin/env Rscript
# Parse a 10x CellRanger feature-barcode .h5, extract the CRISPR Guide Capture
# modality, and save a sparse gRNA-by-cell integer dgCMatrix as grna_matrix.rds.
#
# Usage: Rscript parse_10x_h5_guides.R <input.h5> <output_dir> [guide_feature_type_regex]
#
suppressMessages({library(rhdf5); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
in_h5  <- args[1]
out_dir <- args[2]
ft_regex <- if (length(args) >= 3) args[3] else "CRISPR Guide Capture|Guide Capture|CRISPR"

stopifnot(file.exists(in_h5))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Inspecting:", in_h5, "\n")
ls <- h5ls(in_h5)
print(ls[ls$group %in% c("/", "/matrix", "/matrix/features"), c("group","name","dim")])

# CellRanger v3 layout: /matrix/{data,indices,indptr,shape,barcodes} and
#   /matrix/features/{id,name,feature_type}
grp <- if (any(ls$name == "matrix" & ls$group == "/")) "/matrix" else "/"

read1 <- function(p) as.vector(h5read(in_h5, p))
data    <- read1(paste0(grp, "/data"))
indices <- read1(paste0(grp, "/indices"))
indptr  <- read1(paste0(grp, "/indptr"))
shape   <- read1(paste0(grp, "/shape"))
barcodes<- read1(paste0(grp, "/barcodes"))

# features
fpath <- paste0(grp, "/features")
fnames <- ls[ls$group == fpath, "name"]
feat_id   <- if ("id"   %in% fnames) read1(paste0(fpath, "/id"))   else read1(paste0(grp,"/genes"))
feat_name <- if ("name" %in% fnames) read1(paste0(fpath, "/name")) else feat_id
feat_type <- if ("feature_type" %in% fnames) read1(paste0(fpath, "/feature_type")) else rep("Unknown", length(feat_id))

cat("\nfeature_type table:\n"); print(table(feat_type))

# Build full sparse matrix (features x cells). CellRanger stores CSC: data/
# indices/indptr describe a dgCMatrix of dimension shape (features x cells).
m <- new("dgCMatrix", i = as.integer(indices), p = as.integer(indptr),
         x = as.numeric(data), Dim = as.integer(shape))
rownames(m) <- make.unique(as.character(feat_name))
colnames(m) <- as.character(barcodes)

# Select guide features
is_guide <- grepl(ft_regex, feat_type)
cat("\nn guide features matched:", sum(is_guide), "of", length(feat_type), "\n")
if (sum(is_guide) == 0) stop("No guide features matched regex: ", ft_regex)

grna <- m[is_guide, , drop = FALSE]
rownames(grna) <- make.unique(as.character(feat_id[is_guide]))
grna <- as(grna, "CsparseMatrix")
storage.mode(grna@x) <- "double"

cat("\ngRNA matrix dim (guides x cells):", nrow(grna), "x", ncol(grna), "\n")
cat("total guide UMIs:", sum(grna@x), " nnz:", length(grna@x), "\n")
moi <- mean(colSums(grna > 0))
cat("mean guides-per-cell (>0):", round(moi, 3), "\n")

saveRDS(grna, file.path(out_dir, "grna_matrix.rds"))
cat("Saved:", file.path(out_dir, "grna_matrix.rds"), "\n")

# also save guide feature metadata
write.csv(data.frame(id = feat_id[is_guide], name = feat_name[is_guide],
                     feature_type = feat_type[is_guide]),
          file.path(out_dir, "guide_features.csv"), row.names = FALSE)
