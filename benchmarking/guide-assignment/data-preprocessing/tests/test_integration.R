#!/usr/bin/env Rscript
# ============================================================================
# DORMANT (paused 2026-08-30). Part of the NESTED GUIDE-SUBSET benchmark:
# datasets of 100/200/400/800 guides selected by perturbed sample size, built to
# measure how method memory and runtime SCALE with guide count.
#
# That is a different question from the one currently being run, which is "which
# methods can handle the ENTIRE grna matrix at all" -- see build_full_h5ad.R.
# Nothing outside this directory references these files.
#
# Paused, not abandoned: this code was unit- and integration-tested and works.
# ============================================================================
# Integration test: build a SMALL real dataset into a temp dir and verify the
# written per-method files against the source ODM -- i.e. that the datasets are
# what we say they are. Needs the ODMs + Matrix + zellkonverter (assumed present
# locally). Runs a full stat pass over the chosen dataset, so it is not instant;
# defaults to Replogle (fewer guides).
#
#   Rscript tests/test_integration.R            # replogle-rd7
#   Rscript tests/test_integration.R gasperini  # (slower)

suppressPackageStartupMessages({
  library(Matrix); library(ondisc); library(zellkonverter); library(SummarizedExperiment)
})
source("~/.Rprofile")                          # .get_config_path()

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
source(file.path(here, "..", "convert_odm_to_h5ad.R"))
source(file.path(here, "..", "lib_make_guide_data.R"))

args    <- commandArgs(trailingOnly = TRUE)
DATASET <- if (length(args) >= 1) args[1] else "replogle-rd7"
T_THRESH <- 5L; K_FLOOR <- 10L; SIZES <- c(25L, 50L)   # 50 %% 25 == 0

check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}

# Independent ground truth: read the raw ODM rows for `guides` (does NOT go
# through extract_guide_counts), assemble guides x cells in the given order.
truth_from_odm <- function(grna_odm, guides, n_cells) {
  parts <- lapply(seq_along(guides), function(a) {
    row <- grna_odm[guides[a], ]
    nz  <- which(row > 0)
    if (length(nz) == 0) return(NULL)
    list(i = rep.int(a, length(nz)), j = nz, x = row[nz])
  })
  parts <- Filter(Negate(is.null), parts)
  sparseMatrix(i = unlist(lapply(parts, `[[`, "i")),
               j = unlist(lapply(parts, `[[`, "j")),
               x = unlist(lapply(parts, `[[`, "x")),
               dims = c(length(guides), n_cells))
}

# --- build into a throwaway dir (real source, temp output) ------------------
input_root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
tmp <- file.path(tempdir(), paste0("ga_itest_", DATASET))
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

cat("=== integration build:", DATASET, "-> ", tmp, "===\n")
sel <- build_guide_assignment_datasets(
  source_name = DATASET, out_prefix = "zzz_itest",
  threshold = T_THRESH, K = K_FLOOR, sizes = SIZES,
  methods = c("crispat", "cleanser"),           # pertpy h5ad == crispat h5ad
  input_root = input_root, output_root = tmp
)

# --- independent handle on the source ODM -----------------------------------
grna_odm <- initialize_odm_from_backing_file(
  file.path(input_root, DATASET, "sceptre-pipeline", "grna.odm"))
n_cells    <- ncol(grna_odm)
cell_names <- paste0("CELL_", seq_len(n_cells))

cat("\n=== verifying", length(SIZES), "sizes against the ODM ===\n")
for (m in SIZES) {
  cat("- size", m, "\n")
  guides_m <- sel$guide_sets[[as.character(m)]]
  ds_dir   <- file.path(tmp, paste0("zzz_itest_T", T_THRESH, "_K", K_FLOOR, "_", m))

  # (1) selection claim: right count, clears K floor, spans a real range
  check(length(guides_m) == m, "selected guide count == size")
  stats_m <- vapply(guides_m, function(g) sum(grna_odm[g, ] >= T_THRESH), numeric(1))
  check(all(stats_m >= K_FLOOR), "every selected guide has >= K cells at count >= T (recomputed from ODM)")
  check(min(stats_m) < max(stats_m), "selected guides span a range of perturbed sample sizes")

  # (2) ground truth counts for these guides, straight from the ODM
  truth <- truth_from_odm(grna_odm, guides_m, n_cells)

  # (3) crispat/pertpy h5ad: cells x guides on disk, round-trips to guides x cells
  sce <- readH5AD(file.path(ds_dir, "crispat", "grna_matrix.h5ad"))
  h5  <- assay(sce, 1)
  check(identical(rownames(h5), guides_m),   "h5ad var order == selected guides")
  check(ncol(h5) == n_cells,                 "h5ad has ALL cells")
  check(identical(colnames(h5), cell_names), "h5ad obs names == CELL_<i>")
  check(sum(abs(h5 - truth)) == 0,           "h5ad counts == ODM counts (faithful, correctly oriented)")

  # (4) cleanser mtx: guides x cells, positional order == selection order
  mtx <- readMM(file.path(ds_dir, "cleanser", "grna_matrix.mtx"))
  check(all(dim(mtx) == c(m, n_cells)),      "mtx is guides x cells")
  check(sum(abs(mtx - truth)) == 0,          "mtx counts == ODM counts (== h5ad, cross-encoding consistent)")
}

# (5) nesting realized on disk (from the h5ad var names, not just from sel)
small <- rownames(assay(readH5AD(file.path(tmp, paste0("zzz_itest_T", T_THRESH, "_K", K_FLOOR, "_", min(SIZES)), "crispat", "grna_matrix.h5ad")), 1))
big   <- rownames(assay(readH5AD(file.path(tmp, paste0("zzz_itest_T", T_THRESH, "_K", K_FLOOR, "_", max(SIZES)), "crispat", "grna_matrix.h5ad")), 1))
check(all(small %in% big), "on-disk nesting: smaller size's guides subset of larger's")

cat("\nALL INTEGRATION TESTS PASSED\n")
