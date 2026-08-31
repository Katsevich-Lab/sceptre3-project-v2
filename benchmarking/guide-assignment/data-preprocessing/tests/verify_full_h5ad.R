#!/usr/bin/env Rscript
# Verify that the FULL-matrix grna_matrix.h5ad written by build_full_h5ad.R is
# actually what we claim: read it back OFF DISK and check it against the ODM.
#
#   Rscript tests/verify_full_h5ad.R gasperini
#   Rscript tests/verify_full_h5ad.R replogle-rd7
#
# build_full_h5ad.R already cross-checks dims + total nnz against the cleanser
# .mtx, but that validates the IN-MEMORY matrix. This validates the FILE the
# pipeline will actually hand to crispat and pertpy -- covering the write step,
# the SCE->AnnData orientation flip, and name preservation.
#
# Run under R 4.6 + the R_461 renv library (ondisc + zellkonverter).

suppressPackageStartupMessages({
  library(ondisc); library(Matrix); library(zellkonverter); library(SummarizedExperiment)
})
source("~/.Rprofile")

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else stop("usage: verify_full_h5ad.R <dataset>")
N_PROBE <- if (length(args) >= 2) as.integer(args[2]) else 12

check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}

root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
ds   <- file.path(root, DATASET)

grna_odm <- initialize_odm_from_backing_file(file.path(ds, "sceptre-pipeline", "grna.odm"))
guides   <- rownames(grna_odm); n_cells <- ncol(grna_odm)
cat("\n=== verifying", DATASET, "h5ad against the ODM ===\n")
cat("ODM: ", length(guides), " guides x ", n_cells, " cells\n", sep = "")

for (method in c("crispat", "pertpy")) {
  fp <- file.path(ds, method, "grna_matrix.h5ad")
  cat("-", method, "\n")
  check(file.exists(fp), paste0(method, "/grna_matrix.h5ad exists"))

  h5 <- assay(readH5AD(fp), 1)     # comes back guides x cells
  check(nrow(h5) == length(guides), "h5ad has EVERY guide")
  check(ncol(h5) == n_cells,        "h5ad has EVERY cell")
  check(identical(rownames(h5), guides), "guide names + order match the ODM exactly")
  check(identical(colnames(h5), paste0("CELL_", seq_len(n_cells))), "cell names are CELL_<i>")

  # Values: spot-check random guides straight from the ODM (independent path).
  set.seed(42)
  idx <- sort(sample(seq_along(guides), min(N_PROBE, length(guides))))
  bad <- 0L
  for (i in idx) {
    odm_row <- as.numeric(grna_odm[i, ])
    h5_row  <- as.numeric(h5[i, ])
    if (!identical(odm_row, h5_row)) bad <- bad + 1L
  }
  check(bad == 0L, paste0("all ", length(idx), " probed guides match the ODM value-for-value"))

  # Totals across the WHOLE matrix, not just the probes.
  check(sum(h5) == sum(vapply(seq_along(guides), function(i) sum(grna_odm[i, ]), numeric(1))),
        "grand total of all counts matches the ODM")
}

cat("\nALL FULL-MATRIX H5AD CHECKS PASSED for", DATASET, "\n")
