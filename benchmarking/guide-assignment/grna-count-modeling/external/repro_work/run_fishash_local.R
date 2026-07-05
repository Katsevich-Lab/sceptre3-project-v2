#!/usr/bin/env Rscript
# Faithful to fishash_analysis bin/run_fishash.R (revision): refit=10, padj_cutoff=0.05,
# exclude_empty=TRUE. Input is the guides x cells grna count matrix (.mtx).
suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
})
args <- commandArgs(TRUE)
in_mtx <- args[1]; out_mtx <- args[2]; guides_csv <- args[3]; meta_csv <- args[4]
counts <- as(readMM(in_mtx), "CsparseMatrix")  # guides x cells
rownames(counts) <- read.csv(guides_csv)$guide
colnames(counts) <- read.csv(meta_csv)$barcode
res <- fishash(counts, refit = 10, padj_cutoff = 0.05, exclude_empty = TRUE)
writeMM(assay(res, "assigned"), out_mtx)
cat(sprintf("[%s] assigned pairs = %d\n", in_mtx, sum(assay(res,"assigned"))))
