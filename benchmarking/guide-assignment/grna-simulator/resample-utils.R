#!/usr/bin/env Rscript
# Resampling generator for the computational-scaling datasets. PROTOTYPE.
#
# A simulated dataset of n_guides x n_cells is built from a real gRNA count
# matrix by:
#   1. picking n_cells real cells, uniformly without replacement;
#   2. taking each picked cell's nonzero counts as they are, dropping which
#      guides they were on;
#   3. placing each cell's counts on distinct guides out of n_guides, chosen
#      without replacement, either uniformly or (placement = "weighted") with
#      probability proportional to a per-guide weight. Each simulated guide's
#      weight is the number of nonzeros of a real guide drawn at random, so
#      per-guide loads are as uneven as the real screen's instead of all equal.
#
# placement = "none" skips step 3: each cell keeps its real guides, so the result
# is the real matrix restricted to the drawn cells (n_guides must be the real
# number of guides). With the same seed it draws the same cells as the other
# placements, so it is the baseline for measuring what redrawing guides changes.
#
# Placement is without replacement even when weighted: two entries of a cell
# can never land on the same guide, so nonzeros per cell and MOI stay exactly
# the real cell's. The cost is that a guide can take at most one entry per
# cell, which slightly caps the heaviest guides in cells with many nonzeros.
#
# So every simulated cell is a real cell with its guide labels redrawn: its
# MOI, nonzeros, signal and ambient counts are exactly the real cell's, at any
# n_cells and n_guides. Nothing is calibrated and no threshold is used; the
# threshold below is only for the summary statistics.
#
# Not modelled: guide identity (each simulated guide's column mixes counts from
# many real guides), and collisions -- in a real screen with few guides some
# ambient molecules would land on the same guide, giving fewer ambient nonzeros.
# The second is small for n_guides >= ~500.
#
# Packages are called with `::`, so sourcing this file attaches nothing.

# Real datasets that can be resampled, and the count threshold used to call an
# entry perturbed in the summary statistics (as in measure-real-targets.R).
RESAMPLE_SOURCES <- list(
  gasperini      = list(threshold = 5L),
  `replogle-rd7` = list(threshold = 10L)
)

#' Path of a real dataset's full gRNA count matrix (guides x cells): fishashplus's
#' .rds (fast to read; data-preprocessing/add_fishashplus_rds.R checked it equals
#' the .mtx) or cleanser's .mtx (what cleanser runs on).
#' Needs .get_config_path() (from ~/.Rprofile).
real_counts_path <- function(source, format = c("rds", "mtx")) {
  format <- match.arg(format)
  root <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment", "input_data", source)
  if (format == "rds") file.path(root, "fishashplus", "grna_matrix.rds")
  else file.path(root, "cleanser", "grna_matrix.mtx")
}

#' Read a real dataset's full count matrix as a dgCMatrix (columns are cells),
#' with explicit zeros dropped so each column's entries are its nonzeros.
load_real_counts <- function(source, format = c("rds", "mtx")) {
  if (!source %in% names(RESAMPLE_SOURCES))
    stop("unknown source: ", source, "; one of ", paste(names(RESAMPLE_SOURCES), collapse = ", "))
  fp <- real_counts_path(source, match.arg(format))
  cat("reading", fp, "\n")
  m <- if (grepl("\\.rds$", fp)) readRDS(fp) else Matrix::readMM(fp)
  Matrix::drop0(methods::as(methods::as(m, "CsparseMatrix"), "generalMatrix"))
}

#' Simulate one n_guides x n_cells dataset by resampling cells of `real`.
#'
#' @param real      dgCMatrix, guides x cells, from load_real_counts().
#' @param n_cells   number of cells; at most ncol(real).
#' @param n_guides  number of guides; must be at least the number of nonzeros
#'                  of every picked cell, or it stops (it never drops cells,
#'                  which would lower MOI).
#' @param seed      seed for the cell draw and the guide placement.
#' @param placement "uniform", "weighted" or "none" (see the top of this file).
#' @return list(counts = dgCMatrix n_guides x n_cells, source_cells = the
#'         column of `real` each simulated cell came from, source_guides = the
#'         real guide each simulated guide took its weight from (NULL unless
#'         weighted), guide_weights = those weights (NULL unless weighted)).
resample_counts <- function(real, n_cells, n_guides, seed,
                            placement = c("uniform", "weighted", "none")) {
  stopifnot(methods::is(real, "dgCMatrix"))
  placement <- match.arg(placement)
  n_cells <- as.integer(n_cells); n_guides <- as.integer(n_guides)
  if (n_cells > ncol(real))
    stop("n_cells = ", n_cells, " exceeds the ", ncol(real), " real cells")
  if (placement == "none" && n_guides != nrow(real))
    stop("placement = \"none\" keeps the real guides, so n_guides must be ",
         nrow(real), ", not ", n_guides)

  set.seed(seed)
  # 1. cells. Without replacement for now, so no real cell is used twice and at
  # n_guides = nrow(real) the result is a real cell subset with its guide
  # labels shuffled. Drawing with replacement (a bootstrap) would allow
  # n_cells > ncol(real); revisit once the datasets have been looked at.
  cells <- sample.int(ncol(real), n_cells)
  p <- real@p
  k <- p[cells + 1L] - p[cells]              # nonzeros of each picked cell
  if (max(k) > n_guides)
    stop(sprintf(paste0("n_guides = %d is below the %d nonzeros of the busiest ",
                        "picked cell; raise n_guides (cells are never dropped)"),
                 n_guides, max(k)))
  # 2. their counts, cell by cell, in column order
  x <- real@x[sequence(k, from = p[cells] + 1L)]
  # 3. distinct guides for each cell
  source_guides <- NULL; w <- NULL
  if (placement == "none") {
    i <- real@i[sequence(k, from = p[cells] + 1L)] + 1L     # each entry's real guide
  } else if (placement == "uniform") {
    i <- unlist(lapply(k, function(kk) sample.int(n_guides, kk)), use.names = FALSE)
  } else {
    # Weights: the nonzeros of real guides drawn at random -- without
    # replacement while there are enough real guides, with replacement beyond.
    # Real guides with no nonzeros are left out: a zero weight could never be
    # chosen, and too many of them would leave a cell too few guides to use.
    guide_nnz <- tabulate(real@i + 1L, nbins = nrow(real))
    pool <- which(guide_nnz > 0)
    source_guides <- pool[sample.int(length(pool), n_guides,
                                     replace = n_guides > length(pool))]
    w <- guide_nnz[source_guides]
    i <- unlist(lapply(k, function(kk) sample.int(n_guides, kk, prob = w)), use.names = FALSE)
  }
  j <- rep.int(seq_len(n_cells), k)

  counts <- Matrix::sparseMatrix(
    i = i, j = j, x = x, dims = c(n_guides, n_cells),
    dimnames = list(paste0("grna_", seq_len(n_guides)),   # pipeline convention
                    paste0("CELL_", seq_len(n_cells))))
  # Distinct guides within a cell mean no entry was summed into another.
  stopifnot(length(counts@x) == sum(k))
  list(counts = counts, source_cells = cells, source_guides = source_guides,
       guide_weights = w)
}

#' Dataset id. Keeps the source name, because run_cleanser.py picks its model
#' by matching "gasperini" / "replogle" in the id, and the placement, so uniform
#' and weighted datasets of the same size never overwrite each other.
resample_dataset_id <- function(source, placement, n_guides, n_cells)
  sprintf("rs_%s_%s_g%d_c%d", source, placement, n_guides, n_cells)
