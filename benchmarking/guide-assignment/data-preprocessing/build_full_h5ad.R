#!/usr/bin/env Rscript
# Build the crispat/pertpy input (grna_matrix.h5ad) for the FULL grna matrix of
# one source dataset -- every guide, every cell, no subsetting.
#
#   Rscript build_full_h5ad.R gasperini
#   Rscript build_full_h5ad.R replogle-rd7
#
# Requires R 4.6 + the R_461 renv library (ondisc + zellkonverter):
#   R_LIBS_USER=~/katsevich-lab/R_461/renv/library/linux-ubuntu-jammy/R-4.6/x86_64-pc-linux-gnu
#
# WHY THIS EXISTS: the full datasets already ship cleanser/grna_matrix.mtx, but
# crispat/ and pertpy/ have no input at all, so 4 of the 6 rows in
# full_datasets_config.csv would fail on staging rather than on method capability.
#
# WHY WE DO NOT REBUILD THE .mtx: the existing one was verified against the ODM
# (exact dims, plus 15 random probe rows matching on both nonzero count and sum,
# positionally, for both datasets). Its rows are in ODM order with faithful
# values, so building the h5ad from the same ODM in the same order gives a
# consistent trio. This script re-checks that agreement exhaustively on the
# TOTAL nonzero count -- free, since we have the full matrix in hand anyway.

suppressPackageStartupMessages({ library(ondisc); library(Matrix) })
source("~/.Rprofile")                                  # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
source(file.path(here, "convert_odm_to_h5ad.R"))       # R_to_h5ad()
source(file.path(here, "lib_make_guide_data.R"))       # extract_guide_counts(), validate_guide_dataset()

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else stop("usage: build_full_h5ad.R <gasperini|replogle-rd7>")
METHODS <- c("crispat", "pertpy")                      # cleanser already has its .mtx

root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
ds   <- file.path(root, DATASET)

cat("\n=== full-matrix h5ad build:", DATASET, "===\n")

# --- source ODM -------------------------------------------------------------
grna_odm <- initialize_odm_from_backing_file(file.path(ds, "sceptre-pipeline", "grna.odm"))
guides     <- rownames(grna_odm)
n_cells    <- ncol(grna_odm)
cell_names <- paste0("CELL_", seq_len(n_cells))
cat("Source ODM: ", length(guides), " guides x ", n_cells, " cells\n", sep = "")

# --- extract every guide ----------------------------------------------------
cat("Extracting all guides (this is the slow part)...\n")
t0 <- Sys.time()
counts <- extract_guide_counts(grna_odm, guides, cell_names)
cat("  done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min; ",
    length(counts@x), " nonzeros\n", sep = "")

# allow_empty_rows: the FULL matrix keeps every guide, including any never
# detected in any cell. See the note on validate_guide_dataset().
validate_guide_dataset(counts, cell_names, paste0(DATASET, "_full"),
                       allow_empty_rows = TRUE)
n_empty <- sum(Matrix::rowSums(counts > 0) == 0)
if (n_empty > 0) {
  cat(sprintf("  NOTE: %d guide(s) have zero counts in every cell. KEPT, because the\n", n_empty))
  cat( "        cleanser .mtx contains them too -- dropping them here would give\n")
  cat( "        crispat/pertpy a different guide set than cleanser. If a method\n")
  cat( "        later crashes on this dataset, these are the first suspects.\n")
}

# --- cross-check against the existing cleanser .mtx --------------------------
# Exhaustive on dims + total nnz, complementing the earlier per-row probe.
mtx_fp <- file.path(ds, "cleanser", "grna_matrix.mtx")
if (file.exists(mtx_fp)) {
  hdr <- scan(mtx_fp, what = integer(), skip = 1, nmax = 3, quiet = TRUE)
  ok  <- identical(as.numeric(hdr), as.numeric(c(nrow(counts), ncol(counts), length(counts@x))))
  cat(sprintf("  [xcheck] cleanser .mtx says %d x %d, nnz=%d -> %s\n",
              hdr[1], hdr[2], hdr[3],
              if (ok) "MATCHES the matrix we just built" else "*** MISMATCH ***"))
  if (!ok) stop("Built matrix disagrees with the existing cleanser .mtx; ",
                "the three methods would NOT be receiving the same data.")
} else {
  cat("  [xcheck] no cleanser .mtx to compare against; skipping\n")
}

# --- write the per-method inputs --------------------------------------------
write_h5ad_methods(counts, ds, METHODS)

cat("\n=== Done:", DATASET, "->", ds, "===\n")
for (m in METHODS) {
  fp <- file.path(ds, m, "grna_matrix.h5ad")
  cat(sprintf("  %-8s %s (%.1f MB)\n", m, fp, file.size(fp) / 1048576))
}
