#!/usr/bin/env Rscript
# Run the retained R Poisson Fishash+ implementation on a complete dataset.
# Usage: Rscript run_full_fishashplus.R <dataset-directory> <new-output-directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
input <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
if (dir.exists(output)) stop("Output directory already exists: ", output)
dir.create(output, recursive = TRUE)
self <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
implementation <- file.path(dirname(normalizePath(self)), "fishash_repro_lib.R")
source(implementation)
suppressPackageStartupMessages(library(ondisc))
started <- Sys.time()
cat("Started:", format(started), "\n")
counts_fp <- file.path(input, "sceptre", "grna_matrix.rds")
counts <- as(readRDS(counts_fp), "CsparseMatrix")
odm_fp <- file.path(input, "sceptre-pipeline", "grna.odm")
odm <- initialize_odm_from_backing_file(odm_fp)
stopifnot(identical(dim(counts), dim(odm)),
          identical(rownames(counts), rownames(odm)),
          !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
          all(is.finite(counts@x)), all(counts@x > 0),
          all(counts@x == floor(counts@x)))
# Exhaustive value/order comparison with the full Matrix Market benchmark input.
mtx_fp <- file.path(input, "cleanser", "grna_matrix.mtx")
mtx <- as(readMM(mtx_fp), "CsparseMatrix")
stopifnot(identical(dim(counts), dim(mtx)),
          identical(counts@p, mtx@p), identical(counts@i, mtx@i),
          identical(counts@x, mtx@x))
rm(mtx)
# Independent checks against ODM rows and pipeline cell covariates.
probe <- unique(round(seq(1, nrow(counts), length.out = 25)))
for (i in probe) stopifnot(all(as.numeric(counts[i, ]) == odm[i, ]))
pipeline <- readRDS(file.path(input, "sceptre-pipeline", "sceptre_object.rds"))
stopifnot(all(Matrix::colSums(counts) == pipeline@covariate_data_frame$grna_n_umis))
rm(pipeline, odm)
invisible(gc())
cat("Validated full input:", nrow(counts), "guides x", ncol(counts),
    "cells;", length(counts@x), "nonzeros\n")
cat("Settings: Poisson tail; ambient depth; GS q=0.05; min_count=2; refit=10; exclude_empty=TRUE\n")
# Trace only the existing factor function; do not modify the assignment algorithm.
.fit_iteration <- 0L
trace("rank1_factors", tracer = quote({
  .GlobalEnv$.fit_iteration <- .GlobalEnv$.fit_iteration + 1L
  cat("Fit iteration", .GlobalEnv$.fit_iteration, "at", format(Sys.time()), "\n")
}), print = FALSE)
fit_start <- Sys.time()
assigned <- fishash_plus_poisson(counts, q = 0.05, refit = 10,
                                min_count = 2, cell_margin = "ambient",
                                exclude_empty = TRUE)
fit_seconds <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
untrace("rank1_factors")
assigned <- drop0(as(assigned, "CsparseMatrix"))
stopifnot(identical(dim(assigned), dim(counts)),
          identical(dimnames(assigned), dimnames(counts)),
          all(assigned@x), all(counts[assigned] >= 2))
saveRDS(assigned, file.path(output, "grna_assignment_matrix.rds"))
trip <- summary(assigned)
calls <- data.frame(cell_id = trip$j, grna_id = trip$i)
write.csv(calls, file.path(output, "assignments_fishashplus.csv"), row.names = FALSE)
write.csv(data.frame(grna_id = seq_len(nrow(counts)), guide = rownames(counts)),
          file.path(output, "guide_index.csv"), row.names = FALSE)
barcodes <- read.csv(file.path(input, "sceptre", "covariate_data_frame_with_real_cell_names.csv"))
stopifnot(nrow(barcodes) == ncol(counts), !anyDuplicated(barcodes$X))
write.csv(data.frame(cell_id = seq_len(ncol(counts)), cell = colnames(counts),
                     source_barcode = barcodes$X),
          file.path(output, "cell_index.csv"), row.names = FALSE)
per_cell <- Matrix::colSums(assigned)
per_guide <- Matrix::rowSums(assigned)
threshold <- counts >= 5
overlap <- sum(assigned & threshold)
stats <- list(n_guides = nrow(counts), n_cells = ncol(counts),
              input_nonzeros = length(counts@x), n_calls = sum(assigned),
              mean_calls_per_cell = mean(per_cell), median_calls_per_cell = median(per_cell),
              cells_with_calls = sum(per_cell > 0), guides_with_calls = sum(per_guide > 0),
              threshold5_calls = sum(threshold), shared_with_threshold5 = overlap,
              fishashplus_only = sum(assigned) - overlap,
              threshold5_only = sum(threshold) - overlap)
write.csv(data.frame(grna_id = seq_along(per_guide), guide = rownames(counts),
                     n_assigned_cells = as.numeric(per_guide)),
          file.path(output, "per_guide_summary.csv"), row.names = FALSE)
write.csv(data.frame(cell_id = seq_along(per_cell), cell = colnames(counts),
                     n_assigned_guides = as.numeric(per_cell)),
          file.path(output, "per_cell_summary.csv"), row.names = FALSE)
metadata <- list(method = "fishash+_poisson", q = 0.05, refit = 10,
                 min_count = 2, cell_margin = "ambient", exclude_empty = TRUE,
                 iterations = .fit_iteration, fit_seconds = fit_seconds,
                 started = started, completed = Sys.time(), input = input,
                 input_md5 = tools::md5sum(c(counts_fp, mtx_fp, odm_fp)),
                 implementation = implementation,
                 implementation_md5 = tools::md5sum(implementation),
                 validation = "Full MTX values/order identical; 25 ODM rows identical; all cell UMI totals identical to pipeline covariates",
                 stats = stats, session = sessionInfo())
saveRDS(metadata, file.path(output, "assignment_metadata.rds"))
capture.output(str(metadata), file = file.path(output, "run_summary.txt"))
# Confirm the delivered matrix survives serialization without any changes.
stopifnot(identical(readRDS(file.path(output, "grna_assignment_matrix.rds")), assigned))
print(stats)
cat("Fit seconds:", fit_seconds, "; completed:", format(Sys.time()), "\n")
