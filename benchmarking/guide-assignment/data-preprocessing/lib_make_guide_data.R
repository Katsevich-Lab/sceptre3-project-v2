# lib_make_guide_data.R
#
# PARTLY DORMANT -- this file does double duty, so do not delete it with the rest
# of the paused subset benchmark:
#   LIVE  extract_guide_counts(), validate_guide_dataset(), write_h5ad_methods()
#         are used by build_full_h5ad.R, which builds the CURRENT full-matrix
#         inputs. validate_guide_dataset(allow_empty_rows=) exists for that path.
#   DORMANT  select_nested_guides(), compute_guide_perturbation_stats() and
#         build_guide_assignment_datasets() belong to the paused nested
#         guide-subset benchmark (100/200/400/800 guides by perturbed sample
#         size). Unused at present; kept because that work was tested and may
#         resume.
#
# Build the guide-assignment benchmark input datasets -- crispat + pertpy
# (grna_matrix.h5ad) and cleanser (grna_matrix.mtx) only. SCEPTRE IS DEFERRED:
# no response matrix / covariates / target df / formula, so this needs nothing
# from the association lib -- just the raw grna count matrix.
#
# One parametrized builder, driven by build_gasperini.R / build_replogle.R. For
# each source dataset we write a NESTED sequence of guide-count sizes
# (100/200/400/800) over the SAME full cell universe (all cells), so a
# size-to-size difference in the pipeline's memory/time is purely the guide
# count, not a change in composition.
#
# Guide selection axis = "perturbed sample size": per guide, the number of cells
# with count >= T (default 5). Cells with 1-4 UMIs are mostly ambient, so >= T
# proxies the genuinely-perturbed population. We keep guides with at least K such
# cells (default 10), then pick guides EVENLY SPACED IN RANK across that pool so
# the sizes span the range of perturbed sample sizes.
#
# Deps: ondisc (grna ODM), Matrix, and R_to_h5ad() from ./convert_odm_to_h5ad.R.


# ---- pass 1: per-guide perturbation stat (streaming, low memory) -----------
#' Number of cells with count >= threshold, for every guide, read one row at a
#' time so we never hold the whole grna matrix in memory. Returns a numeric
#' vector named by guide id, in ODM row order.
compute_guide_perturbation_stats <- function(grna_odm, threshold) {
  n_guides  <- nrow(grna_odm)
  guide_ids <- rownames(grna_odm)
  stat <- numeric(n_guides)
  for (i in seq_len(n_guides)) {
    stat[i] <- sum(grna_odm[i, ] >= threshold)   # row read, then discarded
    if (i %% 500 == 0) cat("  statted", i, "/", n_guides, "guides\n")
  }
  names(stat) <- guide_ids
  stat
}


# ---- guide selection: eligible pool -> nested even-coverage sizes ----------
#' Sort the eligible pool (stat >= K) by the stat, take the ANCHOR = max(sizes)
#' guides evenly spaced in rank, then define every smaller size as an exact
#' DECIMATION of that anchor -- guaranteeing size_100 subset of size_200 subset
#' of ... subset of size_max.
#'
#' Convention: the anchor is ALWAYS the final (largest) dataset in `sizes` --
#' there is NO phantom larger anchor. If you later extend `sizes` to a bigger
#' top size, the smaller sets are re-derived from the new, larger anchor and can
#' shift slightly; extend the sequence deliberately.
#'
#' @return list(guide_sets = named list size->guide_ids, pool_size, anchor_size)
select_nested_guides <- function(stat, K, sizes) {
  sizes  <- sort(unique(as.integer(sizes)))
  anchor <- max(sizes)                     # convention: anchor = final dataset
  eligible  <- stat[stat >= K]
  pool_size <- length(eligible)
  cat("Eligible pool: ", pool_size, " guides with stat >= ", K,
      " (of ", length(stat), " total).\n", sep = "")

  if (pool_size < anchor) {
    stop("Eligible pool (", pool_size, ") is smaller than the anchor / largest size (",
         anchor, "). Lower K, lower the top size, or raise T's tolerance.")
  }
  bad <- sizes[anchor %% sizes != 0]       # smaller sizes must divide the anchor
  if (length(bad) > 0) {
    stop("Sizes ", paste(bad, collapse = ", "), " do not evenly divide the anchor (",
         anchor, "); pick sizes that are integer fractions of max(sizes).")
  }

  # Rank the pool by stat (ascending), break ties by guide id for determinism,
  # then take the anchor grid of evenly-spaced ranks.
  ranked_ids   <- names(eligible)[order(eligible, names(eligible))]
  anchor_ranks <- round(seq(1, pool_size, length.out = anchor))
  stopifnot(!any(duplicated(anchor_ranks)))          # even spacing must be 1:1
  anchor_ids   <- ranked_ids[anchor_ranks]

  # Each size is a strided subset of the anchor -> exact nesting.
  guide_sets <- lapply(sizes, function(m) anchor_ids[seq(1, anchor, by = anchor / m)])
  names(guide_sets) <- as.character(sizes)

  for (j in seq_along(guide_sets)[-1]) {
    stopifnot(all(guide_sets[[j - 1]] %in% guide_sets[[j]]))     # strictly nested
  }
  stopifnot(vapply(seq_along(sizes), function(j) length(guide_sets[[j]]) == sizes[j], logical(1)))

  for (m in sizes) {
    s <- stat[guide_sets[[as.character(m)]]]
    cat("  size ", m, ": stat range [", min(s), ", ", max(s), "], median ",
        stats::median(s), "\n", sep = "")
  }
  list(guide_sets = guide_sets, pool_size = pool_size, anchor_size = anchor)
}


# ---- pass 2: extract selected guides x all cells (sparse) ------------------
#' Row-wise ODM subset -> guides x cells dgCMatrix, read one row at a time.
extract_guide_counts <- function(grna_odm, guides, cell_names) {
  ilist <- jlist <- xlist <- vector("list", length(guides))
  for (i in seq_along(guides)) {
    row <- grna_odm[guides[i], ]
    nz  <- which(row > 0)
    if (length(nz) > 0) {
      ilist[[i]] <- rep.int(i, length(nz)); jlist[[i]] <- nz; xlist[[i]] <- row[nz]
    }
  }
  m <- Matrix::sparseMatrix(
    i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
    dims = c(length(guides), length(cell_names))
  )
  rownames(m) <- guides
  colnames(m) <- cell_names
  m
}


# ---- validation (asserts + FYI report) -------------------------------------
#' Empty COLUMNS are expected (a cell need not express any selected guide) so
#' they are reported, not errored. Empty ROWS cannot occur (every selected guide
#' has >= K perturbed cells) but are asserted anyway. No data is modified.
#' @param allow_empty_rows  FALSE (default) for SELECTED guide sets, where the K
#'   floor makes an all-zero guide impossible, so one indicates a real bug. TRUE
#'   for the FULL matrix, which by definition keeps every guide including any
#'   that were never detected (Gasperini has 3 such guides of 13,077). Those are
#'   kept deliberately: the cleanser .mtx contains all rows, so dropping them
#'   from the h5ad would give crispat/pertpy a different guide set than cleanser
#'   and silently misalign positional guide indices between the two encodings.
validate_guide_dataset <- function(grna_counts, cell_names, dataset_id,
                                   allow_empty_rows = FALSE) {
  n_guides <- nrow(grna_counts)
  n_cells  <- ncol(grna_counts)
  stopifnot(!is.null(rownames(grna_counts)), !is.null(colnames(grna_counts)))
  stopifnot(!any(duplicated(rownames(grna_counts))))
  stopifnot(identical(colnames(grna_counts), cell_names))
  stopifnot(length(grna_counts@x) == 0 || min(grna_counts@x) >= 0)

  empty_rows <- sum(Matrix::rowSums(grna_counts > 0) == 0)
  if (!allow_empty_rows) stopifnot(empty_rows == 0)
  empty_cols <- sum(Matrix::colSums(grna_counts) == 0)
  cat(sprintf("  [validate] %s: %d guides x %d cells | empty rows=%d | empty cols=%d (%.1f%%)\n",
              dataset_id, n_guides, n_cells, empty_rows, empty_cols,
              100 * empty_cols / n_cells))
  invisible(TRUE)
}


# ---- per-method writers ----------------------------------------------------
#' crispat/pertpy h5ad (cells x guides). Written once, copied to the sibling
#' method dir -- the two consume an identical grna_matrix.h5ad.
write_h5ad_methods <- function(grna_counts, out_dir, methods) {
  h5ad_methods <- intersect(c("crispat", "pertpy"), methods)
  if (length(h5ad_methods) == 0) return(invisible())

  first_dir <- file.path(out_dir, h5ad_methods[1])
  dir.create(first_dir, recursive = TRUE, showWarnings = FALSE)
  first_fp <- file.path(first_dir, "grna_matrix.h5ad")
  cat("  writing", h5ad_methods[1], "h5ad ->", first_fp, "\n")
  R_to_h5ad(grna_counts, first_fp)                 # transposes to cells x guides

  for (m in h5ad_methods[-1]) {
    m_dir <- file.path(out_dir, m)
    dir.create(m_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(first_fp, file.path(m_dir, "grna_matrix.h5ad"), overwrite = TRUE)
    cat("  copied h5ad ->", file.path(m_dir, "grna_matrix.h5ad"), "\n")
  }
}

#' cleanser Matrix Market file (guides x cells).
write_cleanser_method <- function(grna_counts, out_dir) {
  cdir <- file.path(out_dir, "cleanser")
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  fp <- file.path(cdir, "grna_matrix.mtx")
  cat("  writing cleanser mtx ->", fp, "\n")
  Matrix::writeMM(grna_counts, fp)
}


# ---- top-level builder -----------------------------------------------------
#' Build the nested sequence of guide-assignment inputs for ONE source dataset.
#'
#' @param source_name  ODM source under <input_root>/<source_name>/sceptre-pipeline/
#' @param out_prefix   output dir prefix (must contain "gasperini"/"replogle" --
#'                     the pipeline greps the dataset_id for it)
#' @param threshold    T: count >= T counts a cell as perturbed (default 5)
#' @param K            min perturbed cells for a guide to be eligible (default 10)
#' @param sizes        guide counts to materialize; max(sizes) is the anchor
#' @param methods      which method inputs to write
#' @param input_root   root holding the source ODM (<input_root>/<source_name>/sceptre-pipeline/)
#' @param output_root  root the size datasets are written under (default = input_root)
build_guide_assignment_datasets <- function(
    source_name, out_prefix,
    threshold = 5, K = 10, sizes = c(100, 200, 400, 800),
    methods = c("crispat", "pertpy", "cleanser"),
    input_root = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                           "guide_assignment/input_data"),
    output_root = input_root) {

  stopifnot(all(methods %in% c("crispat", "pertpy", "cleanser")))
  cat("\n=== Building guide-assignment data for", source_name, "===\n")
  cat("T =", threshold, " K =", K, " sizes =", paste(sizes, collapse = ","),
      " anchor =", max(sizes), "(final dataset)  methods:", paste(methods, collapse = ", "), "\n")

  # --- source grna ODM -----------------------------------------------------
  grna_odm <- ondisc::initialize_odm_from_backing_file(
    file.path(input_root, source_name, "sceptre-pipeline", "grna.odm")
  )
  n_cells    <- ncol(grna_odm)
  cell_names <- paste0("CELL_", seq_len(n_cells))     # ODM barcodes not exposed
  cat("Source:", nrow(grna_odm), "guides x", n_cells, "cells.\n")

  # --- select nested guide sets --------------------------------------------
  cat("\n--- Pass 1: statting all guides ---\n")
  stat <- compute_guide_perturbation_stats(grna_odm, threshold)
  sel  <- select_nested_guides(stat, K, sizes)

  # --- extract the largest set once (others are row subsets) ---------------
  cat("\n--- Pass 2: extracting largest set (", max(sizes), " guides x all cells) ---\n", sep = "")
  big_counts <- extract_guide_counts(grna_odm, sel$guide_sets[[as.character(max(sizes))]], cell_names)

  # --- write each size ------------------------------------------------------
  for (m in sizes) {
    counts_m   <- big_counts[sel$guide_sets[[as.character(m)]], , drop = FALSE]   # nested subset
    dataset_id <- paste0(out_prefix, "_T", threshold, "_K", K, "_", m)
    out_dir    <- file.path(output_root, dataset_id)
    cat("\n### ", dataset_id, ": ", nrow(counts_m), " guides x ", ncol(counts_m), " cells ###\n", sep = "")

    validate_guide_dataset(counts_m, cell_names, dataset_id)
    write_h5ad_methods(counts_m, out_dir, methods)
    if ("cleanser" %in% methods) write_cleanser_method(counts_m, out_dir)
  }

  cat("\n=== Done:", source_name, "-> ", output_root, "===\n")
  invisible(sel)
}
