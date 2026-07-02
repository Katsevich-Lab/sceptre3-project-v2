#!/usr/bin/env Rscript
# Prepare a benchmark-sized subset of the comprehensive sim (so geomux/fishash
# are tractable), export it for the external tools, and run our R methods
# (ambient, otsu, valley) on the identical subset.
suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
EXP <- file.path(HERE, "results", "comprehensive_bench"); dir.create(EXP, recursive=TRUE, showWarnings=FALSE)

gm <- readRDS(file.path(DATA,"sims_comprehensive/sceptre/grna_matrix.rds"))
truth <- readRDS(file.path(DATA,"sims_comprehensive/true_pert_matrix.rds"))
gn <- rownames(gm)
group <- vapply(strsplit(gn,"_"), function(p) paste(p[2],p[3],p[4],sep="_"), character(1))

set.seed(2)
# stratified: 40 guides per group; 12000 cells
gsel <- unlist(lapply(split(seq_along(gn), group), function(ix) sample(ix, min(40, length(ix)))))
gsel <- sort(gsel)
csel <- sort(sample(ncol(gm), 12000))
gm2 <- as(gm[gsel, csel], "CsparseMatrix"); tr2 <- as(truth[gsel, csel], "CsparseMatrix")
guides <- data.frame(guide = rownames(gm2), group = group[gsel])
cells  <- data.frame(barcode = colnames(gm2))
cat("benchmark subset:", nrow(gm2), "guides x", ncol(gm2), "cells\n")

Matrix::writeMM(gm2, file.path(EXP,"guide_counts.mtx"))           # guides x cells
write.csv(cells,  file.path(EXP,"cells.csv"),  row.names=FALSE)
write.csv(guides, file.path(EXP,"guides.csv"), row.names=FALSE)
saveRDS(tr2, file.path(EXP,"truth.rds")); saveRDS(gm2, file.path(EXP,"gm.rds"))

# our R methods
saveRDS(ambient_test_assign(gm2, q=0.05, model="hypergeometric", n_iter=1)$assignment_matrix, file.path(EXP,"assign_ambient.rds"))
saveRDS(assign_by_threshold(gm2, otsu_threshold_log1p)$assignment_matrix,       file.path(EXP,"assign_otsu.rds"))
saveRDS(assign_by_threshold(gm2, smoothed_valley_threshold)$assignment_matrix,  file.path(EXP,"assign_valley.rds"))
cat("exported to", EXP, "and ran ambient/otsu/valley\n")
