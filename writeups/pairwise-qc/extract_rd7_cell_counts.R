#!/usr/bin/env Rscript
# Recover the post-QC treatment and non-targeting cell counts used by extract_rd7.R.
suppressMessages({ library(sceptre); library(Matrix) })

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
outdir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
rd7_dir <- "/Users/ekatsevi/data/projects/sceptre3/nf_pipelines/rd7_trans"

so <- readRDS(file.path(rd7_dir, "sceptre_object.rds"))
ciu <- tryCatch(so@cells_in_use, error = function(e) NULL)
am <- readRDS(file.path(rd7_dir, "sceptre_outputs", "grna_assignment_matrix.rds"))
gtdf <- so@grna_target_data_frame
if (length(ciu) > 0 && max(ciu) <= ncol(am)) am <- am[, ciu, drop = FALSE]

target_of <- gtdf$grna_target[match(rownames(am), gtdf$grna_id)]
targets <- setdiff(unique(target_of), "non-targeting")
n_trt <- sapply(targets, function(target) {
  sum(Matrix::colSums(am[which(target_of == target), , drop = FALSE]) > 0)
})
n_nt <- sum(Matrix::colSums(am[which(target_of == "non-targeting"), , drop = FALSE]) > 0)

saveRDS(
  list(n_trt = setNames(as.integer(n_trt), targets), n_nt = as.integer(n_nt)),
  file.path(outdir, "rd7_ntrt.rds")
)
cat("Saved", file.path(outdir, "rd7_ntrt.rds"), "\n")
