#!/usr/bin/env Rscript
# Tests for resample_counts() in resample-utils.R.
#
#   Rscript tests/test-resample.R      # from grna-simulator/; on Betty via run-in-image.sh
#
# Uses a small synthetic "real" matrix, so it needs no data and runs in seconds.

suppressPackageStartupMessages(library(Matrix))

.self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here  <- if (length(.self) == 1) dirname(normalizePath(sub("^--file=", "", .self))) else getwd()
source(file.path(normalizePath(file.path(here, "..")), "resample-utils.R"))

check <- function(cond, label) {
  if (!isTRUE(cond)) stop("FAIL: ", label, call. = FALSE)
  cat("  ok:", label, "\n")
}

# A "real" matrix: 60 guides x 400 cells, 0-25 nonzeros per cell, counts 1-500,
# including some cells with no nonzeros at all.
set.seed(7)
G0 <- 60L; N0 <- 400L
k0 <- sample(0:25, N0, replace = TRUE)
real <- sparseMatrix(
  i = unlist(lapply(k0, function(k) sample.int(G0, k))),
  j = rep.int(seq_len(N0), k0),
  x = as.numeric(sample.int(500L, sum(k0), replace = TRUE)),
  dims = c(G0, N0))
cell_counts <- function(m, j) sort(m@x[(m@p[j] + 1L):m@p[j + 1L]][seq_len(m@p[j + 1L] - m@p[j])])


cat("\n[1] every simulated cell is its source cell, with guides redrawn\n")
N <- 250L; G <- 40L
sim <- resample_counts(real, N, G, seed = 1)
m <- sim$counts
check(identical(dim(m), c(G, N)), "dimensions are n_guides x n_cells")
check(!anyDuplicated(sim$source_cells), "no real cell is used twice")
check(all(diff(m@p) == diff(real@p)[sim$source_cells]),
      "each cell's number of nonzeros equals its source cell's")
check(all(vapply(seq_len(N), function(j)
        identical(cell_counts(m, j), cell_counts(real, sim$source_cells[j])), logical(1))),
      "each cell's counts are exactly its source cell's counts")
check(isTRUE(all.equal(unname(colSums(m)), unname(colSums(real)[sim$source_cells]))),
      "each cell's UMI total equals its source cell's")
check(sum(real[, sim$source_cells] >= 5) == sum(m >= 5),
      "thresholded MOI equals the source cells' (at any threshold, e.g. 5)")
check(identical(rownames(m)[1:2], c("grna_1", "grna_2")) &&
        identical(colnames(m)[1:2], c("CELL_1", "CELL_2")),
      "pipeline dimnames")


cat("\n[2] n_guides only changes where counts land\n")
for (Gi in c(max(k0), 100L, 1000L)) {
  mi <- resample_counts(real, N, Gi, seed = 1)$counts
  check(all(diff(mi@p) == diff(real@p)[sim$source_cells]),
        sprintf("G = %d: nonzeros per cell unchanged", Gi))
}


cat("\n[3] refusals\n")
check(inherits(try(resample_counts(real, N, max(k0) - 1L, seed = 1), silent = TRUE), "try-error"),
      "n_guides below the busiest cell's nonzeros stops (cells are never dropped)")
check(inherits(try(resample_counts(real, N0 + 1L, G, seed = 1), silent = TRUE), "try-error"),
      "n_cells above the number of real cells stops")


cat("\n[4] seeds\n")
check(identical(resample_counts(real, N, G, seed = 1), sim), "same seed, identical result")
check(!identical(resample_counts(real, N, G, seed = 2)$counts, m), "different seed, different result")


cat("\n[5] at the real n_guides and n_cells, it is the real matrix with guides shuffled\n")
full <- resample_counts(real, N0, G0, seed = 3)$counts
check(isTRUE(all.equal(sort(full@x), sort(real@x))), "same multiset of counts")
check(length(full@x) == length(real@x), "same number of nonzeros")


cat("\n[6] weighted placement\n")
# A real matrix whose guides are very uneven: guide g's share grows with g, and
# guides 1-5 have no nonzeros at all.
set.seed(11)
G1 <- 80L; N1 <- 3000L
gw <- c(rep(0, 5), seq(1, 10, length.out = G1 - 5L))
k1 <- sample(1:12, N1, replace = TRUE)
uneven <- sparseMatrix(
  i = unlist(lapply(k1, function(k) sample.int(G1, k, prob = gw))),
  j = rep.int(seq_len(N1), k1),
  x = as.numeric(sample.int(50L, sum(k1), replace = TRUE)),
  dims = c(G1, N1))
wt <- resample_counts(uneven, 2500L, 40L, seed = 4, placement = "weighted")
mw <- wt$counts
check(all(diff(mw@p) == diff(uneven@p)[wt$source_cells]),
      "each cell's number of nonzeros equals its source cell's")
check(all(vapply(seq_len(2500L), function(j)
        identical(cell_counts(mw, j), cell_counts(uneven, wt$source_cells[j])), logical(1))),
      "each cell's counts are exactly its source cell's counts")
check(length(wt$source_guides) == 40L && !anyDuplicated(wt$source_guides),
      "each simulated guide takes its weight from a different real guide")
check(!any(wt$source_guides %in% 1:5), "real guides with no nonzeros are never used")
check(identical(wt$guide_weights, tabulate(uneven@i + 1L, G1)[wt$source_guides]),
      "weights are the real guides' nonzero counts")
load <- Matrix::rowSums(mw > 0)
check(cor(load, wt$guide_weights) > 0.9,
      sprintf("per-guide load tracks the weights (cor %.3f)", cor(load, wt$guide_weights)))
uni <- Matrix::rowSums(resample_counts(uneven, 2500L, 40L, seed = 4)$counts > 0)
check(sd(load) / mean(load) > 2 * sd(uni) / mean(uni),
      "weighted loads are far more uneven than uniform ones")
check(is.null(resample_counts(uneven, 100L, 40L, seed = 4)$source_guides),
      "uniform placement records no source guides")
check(identical(resample_counts(real, N, G, seed = 1, placement = "uniform"), sim),
      "placement = \"uniform\" is the default")
many <- resample_counts(uneven, 500L, 200L, seed = 5, placement = "weighted")
check(length(many$source_guides) == 200L && all(many$source_guides > 5),
      "more guides than the real screen has: weights drawn with replacement, never from empty guides")

cat("\n[7] placement = \"none\" keeps the real guides\n")
nn <- resample_counts(real, N, G0, seed = 1, placement = "none")
check(identical(unname(as.matrix(nn$counts)), unname(as.matrix(real[, nn$source_cells]))),
      "the result is exactly the real matrix restricted to the drawn cells")
check(identical(nn$source_cells, sim$source_cells) &&
        identical(nn$source_cells, resample_counts(real, N, G0, seed = 1, placement = "weighted")$source_cells),
      "same seed, same cells as the uniform and weighted placements")
check(is.null(nn$source_guides), "no source guides or weights recorded")
check(inherits(try(resample_counts(real, N, G0 - 1L, seed = 1, placement = "none"), silent = TRUE),
               "try-error"),
      "n_guides other than the real number of guides stops")

cat("\nAll tests passed.\n")
