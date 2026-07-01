## For the Replogle benchmarking datasets, I use max assignment, but
## that can call junk cells with no large UMIs perturbed.
## This script computes and writes a list of cell position idx for the cells that
## have max < `MIN_UMI_TO_KEEP`, so that those cells can be considered
## unperturbed in benchmarking dataset creation.

rm(list=ls())
library(sceptre)
source("~/.Rprofile")

dataset <- "replogle-rd7"
assignment_method = "maximum"

MIN_UMI_TO_KEEP = 5  # the max must be >= this or we call the cell unperturbed

# shared normalizer for the cell-universe fingerprint (see fingerprint contract).
# as.integer so an ODM-derived double total and the stored integer grna_n_umis
# hash identically -- only the values-in-order matter, not the R type.
stopifnot(requireNamespace("digest", quietly = TRUE))
fingerprint_cells <- function(x) digest::digest(as.integer(x))

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data",
  dataset,
  "sceptre-pipeline"
)

path_to_outputs <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs",
  dataset,
  assignment_method
)

grna_odm <- file.path(path_to_data, "grna.odm") |> ondisc::initialize_odm_from_backing_file()

cells_too_small <- rep(TRUE, ncol(grna_odm))
cell_sums <- numeric(ncol(grna_odm))
for(i in 1:nrow(grna_odm)) {
  y <- grna_odm[i,]
  cells_too_small[y >= MIN_UMI_TO_KEEP] <- FALSE
  cell_sums <- cell_sums + y
}

cat(sum(cells_too_small), " cells do NOT pass our min UMI threshold of ", MIN_UMI_TO_KEEP, ".\n", sep="")
cat(sum(!cells_too_small), "cells remain.\n")



## verify the cell universe matches the saved assignment before trusting these idx ~~~~
# cell_sums is the per-cell total guide UMI computed straight from the ODM; its
# fingerprint must equal the one the assignment script saved (which it computed from
# the covariate frame). A match means the ODM columns these idx point at are the SAME
# universe, in the same order, that the assignment matrix uses. (Replogle has no cell
# names, so this fingerprint stands in for a cell-name check.)
meta_fp <- file.path(path_to_outputs, "assignment_metadata.rds")
if(!file.exists(meta_fp)) {
  stop("assignment_metadata.rds not found -- run make_benchmarking_assns_repl.R first.")
}
meta <- readRDS(meta_fp)

fingerprint_now <- fingerprint_cells(cell_sums)
if(!identical(meta$cell_fingerprint, fingerprint_now)) {
  stop("cell-universe fingerprint mismatch: this grna ODM does not match the saved ",
       "assignment universe. Do NOT use these indices.")
}
stopifnot(meta$n_cells == ncol(grna_odm))
cat("Cell-universe fingerprint matches the saved assignment.\n")

## write the exclusion list, stamped with the same fingerprint ~~~~~~~~~~~~~~~~~~~~~~~~~
excluded_idx <- which(cells_too_small)
out_fp <- file.path(path_to_outputs, "excluded_cells.rds")
saveRDS(
  list(excluded_idx     = excluded_idx,
       n_excluded       = length(excluded_idx),
       criterion        = paste0("max guide UMI < ", MIN_UMI_TO_KEEP),
       min_umi_to_keep  = MIN_UMI_TO_KEEP,
       n_cells          = ncol(grna_odm),
       cell_fingerprint = fingerprint_now,
       date             = Sys.time()),
  out_fp
)
cat("Wrote", length(excluded_idx), "excluded cell idx to", out_fp, "\n")


