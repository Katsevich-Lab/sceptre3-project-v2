#!/usr/bin/env Rscript
# Add fishashplus's input, <dataset>/fishashplus/grna_matrix.rds, to existing
# real guide-assignment datasets. One-off: the datasets predate fishashplus
# having its own input directory.
#
#   Rscript add_fishashplus_rds.R                 # gasperini and replogle-rd7
#   Rscript add_fishashplus_rds.R replogle-rd7    # one
#
# fishashplus needs guide and cell names, which the canonical input,
# cleanser/grna_matrix.mtx, does not store. So the matrix is taken from the
# dataset's existing sceptre/grna_matrix.rds, which carries the real guide IDs
# and cell barcodes -- but only after checking that it holds exactly the .mtx's
# matrix (same dimensions and every entry), so fishashplus sees what cleanser and
# fishash see. Nothing else reads sceptre/. An existing fishashplus/ file is
# never overwritten.

suppressPackageStartupMessages(library(Matrix))
source("~/.Rprofile")                                    # .get_config_path()

input_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data")
picked <- commandArgs(trailingOnly = TRUE)
datasets <- if (length(picked)) picked else c("gasperini", "replogle-rd7")

as_dgc <- function(m) drop0(as(as(m, "CsparseMatrix"), "generalMatrix"))

failed <- character(0)
for (ds in datasets) {
  cat(sprintf("\n=== %s ===\n", ds))
  src <- file.path(input_root, ds, "sceptre", "grna_matrix.rds")
  mtx <- file.path(input_root, ds, "cleanser", "grna_matrix.mtx")
  out <- file.path(input_root, ds, "fishashplus", "grna_matrix.rds")
  if (file.exists(out)) { cat("  already exists, left alone:", out, "\n"); next }

  cat("  reading", src, "\n"); a <- as_dgc(readRDS(src))
  cat("  reading", mtx, "\n"); b <- as_dgc(readMM(mtx))
  if (is.null(rownames(a)) || is.null(colnames(a)) ||
      anyNA(rownames(a)) || anyNA(colnames(a)) ||
      anyDuplicated(rownames(a)) || anyDuplicated(colnames(a))) {
    cat("  FAIL: the .rds lacks complete, unique guide and cell names\n")
    failed <- c(failed, ds); next
  }
  same <- identical(dim(a), dim(b)) && identical(a@p, b@p) && identical(a@i, b@i) &&
          identical(as.numeric(a@x), as.numeric(b@x))
  cat(sprintf("  %d guides x %d cells, %d nonzeros: .rds and .mtx %s\n",
              nrow(a), ncol(a), length(a@x), if (same) "IDENTICAL" else "DIFFER"))
  if (!same) { failed <- c(failed, ds); next }

  dir.create(dirname(out), showWarnings = FALSE)
  saveRDS(a, out)
  cat("  wrote", out, "\n")
}
if (length(failed)) stop("not written: ", paste(failed, collapse = ", "))
