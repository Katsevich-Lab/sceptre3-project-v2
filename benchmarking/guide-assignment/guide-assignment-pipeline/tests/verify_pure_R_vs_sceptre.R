#!/usr/bin/env Rscript
# Verify that script_baseline (pure-R reimplementation of sceptre's mixture-
# model gRNA assignment) produces identical assignments to the sceptre
# package, on small unit-test datasets.
#
# Run after:
#   nextflow run main.nf --run_id script_unit_test_pure_R_sceptre

suppressPackageStartupMessages(library(Matrix))

source("~/.Rprofile")
OUTPUTS_DIR <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/script_unit_test_pure_R_sceptre"
)

DATASETS <- c(
  "gasperini_assign_nguides=100_ncells=100k",
  "replogle-rd7_assign_nguides=100_ncells=100k"
)

assn_path <- function(method, dataset) {
  file.path(OUTPUTS_DIR, sprintf("assignment_matrix_%s_%s.rds", method, dataset))
}
extras_path <- function(method, dataset) {
  file.path(OUTPUTS_DIR, sprintf("assignment_extras_%s_%s.rds", method, dataset))
}

# Jaccard similarity over two index sets; returns 1 if both are empty.
jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0L) 1 else length(intersect(a, b)) / u
}

report_mismatches <- function(jacs, view_label) {
  bad <- which(jacs < 1)
  if (length(bad) == 0L) return(invisible())
  stop(sprintf(
    "%s mismatch on %d/%d guides (e.g. %s)",
    view_label, length(bad), length(jacs),
    paste(head(names(jacs)[bad], 5L), collapse = ", ")
  ))
}

verify_dataset <- function(dataset) {
  cat(sprintf("== %s ==\n", dataset))

  assns_R    <- readRDS(assn_path("script_baseline", dataset))
  assns_scep <- readRDS(assn_path("sceptre",         dataset))
  extras_R   <- readRDS(extras_path("script_baseline", dataset))

  guides <- rownames(assns_scep)

  # 1. Guide sets match across all three outputs.
  stopifnot(
    "assignment-matrix row names differ between pure-R and sceptre" =
      setequal(rownames(assns_R), guides),
    "extras$per_guide names differ from sceptre's row names" =
      setequal(names(extras_R$per_guide), guides)
  )

  # 2. Per-guide assignments identical (matrix view).
  jacs_mat <- vapply(guides, function(g) {
    jaccard(which(assns_R[g, ]), which(assns_scep[g, ]))
  }, numeric(1))
  report_mismatches(jacs_mat, "matrix-view")

  # 3. Per-guide assignments identical (extras view -- guards the matrix-
  # building code in baseline.R against drifting from the per-guide list).
  jacs_extras <- vapply(guides, function(g) {
    jaccard(extras_R$per_guide[[g]]$assignments, which(assns_scep[g, ]))
  }, numeric(1))
  report_mismatches(jacs_extras, "extras-view")

  cat(sprintf("  PASS: all %d guides match (matrix + extras views)\n",
              length(guides)))
}

for (ds in DATASETS) verify_dataset(ds)
cat("\nAll datasets verified.\n")
