#!/usr/bin/env Rscript
# Unit tests for the pure logic in lib_make_guide_data.R. Data-free (base R +
# Matrix only) -- run instantly, no ODMs / no qsub:
#
#   Rscript tests/test_unit.R
#
# lib_make_guide_data.R only DEFINES functions at source time, so sourcing it
# pulls in the functions under test without running any build.

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
source(file.path(here, "..", "lib_make_guide_data.R"))

# --- tiny harness -----------------------------------------------------------
check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}
expect_error <- function(expr, label) {
  errored <- FALSE
  tryCatch(force(expr), error = function(e) errored <<- TRUE)
  if (!errored) stop("FAIL (expected an error): ", label, call. = FALSE)
  cat("  ok (errored as expected):", label, "\n")
}
mk_stat <- function(vals) setNames(as.numeric(vals), paste0("g", seq_along(vals)))
suf     <- function(ids) as.integer(sub("^g", "", ids))   # guide id suffix == rank here


# ===========================================================================
cat("select_nested_guides: counts, nesting, anchor convention\n")
# stat = 1..100 so guide g<r> has stat r -> rank r; nothing ties.
stat <- mk_stat(1:100)
sel  <- select_nested_guides(stat, K = 1, sizes = c(25, 50))
gs   <- sel$guide_sets

check(sel$pool_size   == 100, "pool = all 100 (K=1)")
check(sel$anchor_size == 50,  "anchor == max(sizes) == 50 (final dataset)")
check(length(gs[["25"]]) == 25 && length(gs[["50"]]) == 50, "each size has its exact count")
check(all(gs[["25"]] %in% gs[["50"]]),                      "size 25 nested in size 50")

cat("select_nested_guides: EVEN COVERAGE across the eligible range\n")
# The anchor is exactly the evenly-spaced rank grid, and it reaches both ends.
expected_anchor_ranks <- as.integer(round(seq(1, 100, length.out = 50)))
check(identical(sort(suf(gs[["50"]])), sort(expected_anchor_ranks)),
      "anchor guides = evenly-spaced rank grid over the pool")
check(min(suf(gs[["50"]])) == 1 && max(suf(gs[["50"]])) == 100,
      "anchor spans the min and max eligible stat")
d <- diff(sort(suf(gs[["50"]])))
check(max(d) - min(d) <= 1, "anchor rank gaps are uniform (evenly spaced)")

cat("select_nested_guides: determinism incl. tie-breaking by id\n")
stat_ties <- mk_stat(c(rep(5, 30), rep(9, 30)))            # heavy ties
a <- select_nested_guides(stat_ties, K = 1, sizes = c(15, 30))$guide_sets
b <- select_nested_guides(stat_ties, K = 1, sizes = c(15, 30))$guide_sets
check(identical(a, b), "same inputs -> identical selection")

cat("select_nested_guides: K floor filters the pool\n")
selK <- select_nested_guides(mk_stat(1:100), K = 51, sizes = c(25, 50))
check(selK$pool_size == 50,                         "pool = guides with stat >= 51 (50 of 100)")
check(all(suf(selK$guide_sets[["50"]]) >= 51),      "every selected guide clears the K floor")

cat("select_nested_guides: guards\n")
expect_error(select_nested_guides(mk_stat(1:30),  K = 1, sizes = c(20, 40)),
             "pool (30) < anchor (40) errors")
expect_error(select_nested_guides(mk_stat(1:100), K = 1, sizes = c(30, 40)),
             "non-divisible sizes (40 %% 30) error")


# ===========================================================================
cat("compute_guide_perturbation_stats: counts cells with value >= threshold\n")
# Only uses nrow()/rownames()/`[i, ]`, so a plain matrix stands in for the ODM.
m <- rbind(
  gA = c(0, 5, 6, 4, 9),     # >=5 : cells 2,3,5      -> 3
  gB = c(1, 2, 3, 4, 4),     # >=5 : none             -> 0
  gC = c(5, 5, 5, 5, 5),     # >=5 : all 5            -> 5
  gD = c(4, 5, 4, 5, 4)      # boundary: 4 excluded, 5 counted -> 2
)
st <- compute_guide_perturbation_stats(m, threshold = 5)
check(identical(names(st), c("gA", "gB", "gC", "gD")), "stat named by guide id, in row order")
check(st[["gA"]] == 3 && st[["gB"]] == 0 && st[["gC"]] == 5 && st[["gD"]] == 2,
      "per-row >= T counts correct (T boundary: ==T counts, ==T-1 does not)")


# ===========================================================================
cat("extract_guide_counts: correct guides x cells, order, names, values\n")
cell_names <- paste0("CELL_", 1:5)
ex <- extract_guide_counts(m, guides = c("gC", "gA"), cell_names = cell_names)  # subset + reorder
check(inherits(ex, "CsparseMatrix"),               "returns a sparse (dgC) matrix")
check(identical(rownames(ex), c("gC", "gA")),      "rownames = requested guides, in order")
check(identical(colnames(ex), cell_names),         "colnames = cell_names")
check(sum(abs(as.matrix(ex) - m[c("gC", "gA"), ])) == 0, "values match the source rows exactly")


# ===========================================================================
cat("validate_guide_dataset: passes valid, catches an empty row\n")
good <- extract_guide_counts(m, guides = c("gA", "gC", "gD"), cell_names = cell_names)
check(isTRUE(validate_guide_dataset(good, cell_names, "good")), "valid dataset passes")
mz   <- rbind(gA = c(0, 5, 6, 4, 9), gZERO = c(0, 0, 0, 0, 0))          # gZERO: all-zero row
bad  <- extract_guide_counts(mz, guides = c("gA", "gZERO"), cell_names = cell_names)
expect_error(validate_guide_dataset(bad, cell_names, "bad"), "empty (all-zero) row is rejected")

check(isTRUE(validate_guide_dataset(bad, cell_names, "bad", allow_empty_rows = TRUE)),
      "empty row is ACCEPTED when allow_empty_rows=TRUE (full-matrix build)")

cat("\nALL UNIT TESTS PASSED\n")
