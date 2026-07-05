# lib/io.R
# I/O and output-contract layer for association benchmarking dataset prep.
#
# Contains, in order:
#   - fingerprint_cells()            : shared order-sensitive cell-universe hash
#   - read_source_data()             : load ODMs + covariates + target df + gene stats
#   - load_assignment()              : load assignment matrix, verify universe, zero excluded cells
#   - compute_genes_passing_qc()     : gene QC (verbatim old formula; NT included)
#   - odm_to_sparse_matrix()         : row-wise ODM subset -> sparse matrix  [ported verbatim]
#   - prepare_frperturb_covariates() : FR-Perturb covariate prep            [ported verbatim]
#   - setup_anndata_environment()    : conda env for h5ad writing           [ported verbatim]
#   - write_sceptre_output()         : SCEPTRE format writer                [ported verbatim]
#   - write_mixscale_output()        : Mixscale format writer               [ported verbatim]
#   - write_frperturb_output()       : FR-Perturb (h5ad) format writer      [ported verbatim]
#
# The three writers + odm_to_sparse_matrix + prepare_frperturb_covariates are
# copied UNCHANGED from the legacy utils_data_prep.R: they define the on-disk
# format the downstream pipeline consumes, and keeping them byte-identical is
# what makes the golden-master diff meaningful.


# ---- cell-universe fingerprint --------------------------------------------
# Order-sensitive hash of a per-cell quantity, used to detect a reordered /
# re-exported cell universe for Replogle (which has no cell barcodes). MUST
# match the definition in the assign-for-assoc scripts. as.integer() collapses
# the double/integer split so an ODM-derived total and the stored integer
# grna_n_umis hash identically -- only values-in-order matter.
stopifnot(requireNamespace("digest", quietly = TRUE))
fingerprint_cells <- function(x) digest::digest(as.integer(x))


# ---- read raw source data -------------------------------------------------
#' Load the ODM-backed sceptre object, ODMs, covariates, target df, gene stats.
#'
#' @param dataset "gasperini" or "replogle-rd7"
#' @return list(scep, response_odm, grna_odm, cell_covariates, grna_target_df,
#'              gene_summary_stats, low_moi)
read_source_data <- function(dataset) {
  path_to_data <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "guide_assignment/input_data", dataset, "sceptre-pipeline"
  )

  scep <- sceptre::read_ondisc_backed_sceptre_object(
    sceptre_object_fp    = file.path(path_to_data, "sceptre_object.rds"),
    response_odm_file_fp = file.path(path_to_data, "response.odm"),
    grna_odm_file_fp     = file.path(path_to_data, "grna.odm")
  )

  response_odm <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
  grna_odm     <- ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
  gene_summary_stats <- readr::read_csv(
    file.path(path_to_data, "gene_summary_stats.csv"), show_col_types = FALSE
  )

  list(
    scep               = scep,
    response_odm       = response_odm,
    grna_odm           = grna_odm,
    cell_covariates    = scep@covariate_data_frame,
    grna_target_df     = scep@grna_target_data_frame,
    gene_summary_stats = gene_summary_stats,
    low_moi            = scep@low_moi
  )
}


# ---- load assignment + apply global cell exclusion ------------------------
#' Load a saved guide-assignment matrix, verify it against the current cell
#' universe, and zero out the columns of any globally-excluded (junk) cells.
#'
#' This is the SINGLE place the exclusion is applied. Zeroing (not dropping)
#' preserves column positions, so cell_idx keeps meaning "column i of the
#' original ODM" -- required by odm_to_sparse_matrix() and cells_to_remove.csv.
#' After this, no selector needs to know about exclusion: an excluded cell
#' expresses no guide, so it can never be selected.
#'
#' @param dataset            "gasperini" or "replogle-rd7"
#' @param assignment_method  subdir under outputs/<dataset>/ (e.g. "maximum",
#'                           "thresholding-5")
#' @param cell_covariates    the covariate frame from read_source_data(); used
#'                           to verify the cell universe (fingerprint / barcodes)
#' @param excluded_required  if TRUE, an excluded_cells.rds MUST be present
#' @return list(assn, assn_raw, excluded_idx, usable_cells, metadata)
#'   assn      : assignment matrix with excluded columns zeroed (use everywhere)
#'   assn_raw  : the matrix as loaded, no zeroing (for gene QC, to match legacy)
load_assignment <- function(dataset, assignment_method, cell_covariates,
                            excluded_required = FALSE) {
  assign_dir <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "guide_assignment/outputs", dataset, assignment_method
  )

  assn <- readRDS(file.path(assign_dir, "grna_assignment_matrix.rds"))
  meta <- readRDS(file.path(assign_dir, "assignment_metadata.rds"))

  # --- verify the assignment matches the current cell universe -------------
  stopifnot(meta$n_cells == ncol(assn))
  stopifnot(nrow(cell_covariates) == ncol(assn))

  # Gasperini: cells have barcodes -> check colnames against the covariate rows.
  if (!is.null(colnames(assn))) {
    if (!identical(colnames(assn), rownames(cell_covariates))) {
      stop("load_assignment: assignment colnames do not match cell_covariates rownames.")
    }
  }
  # Replogle: no barcodes -> the metadata carries an order-sensitive fingerprint.
  if (!is.null(meta$cell_fingerprint)) {
    if (!identical(meta$cell_fingerprint, fingerprint_cells(cell_covariates$grna_n_umis))) {
      stop("load_assignment: cell-universe fingerprint mismatch vs assignment_metadata. ",
           "Do NOT use these indices.")
    }
  }

  # --- load + verify the global exclusion list -----------------------------
  excl_fp <- file.path(assign_dir, "excluded_cells.rds")
  if (file.exists(excl_fp)) {
    excl <- readRDS(excl_fp)
    stopifnot(excl$n_cells == ncol(assn))
    if (!is.null(excl$cell_fingerprint)) {
      if (!identical(excl$cell_fingerprint, fingerprint_cells(cell_covariates$grna_n_umis))) {
        stop("load_assignment: excluded_cells.rds fingerprint mismatch vs cell universe.")
      }
    }
    excluded_idx <- as.integer(excl$excluded_idx)
    stopifnot(all(excluded_idx >= 1L), all(excluded_idx <= ncol(assn)))
  } else {
    if (excluded_required) {
      stop("load_assignment: excluded_cells.rds required but not found in ", assign_dir)
    }
    excluded_idx <- integer(0)
  }

  # --- zero the excluded columns (single mutation point) -------------------
  # sceptre's get_grna_assignments() returns an ngRMatrix: a PATTERN matrix
  # (no stored values) in row-compressed (CSR) storage. Subassigning FALSE into
  # a pattern matrix is fragile, and CSR is the wrong direction for column ops,
  # so coerce ONCE to a logical column-compressed matrix (lgCMatrix). colSums,
  # subsetting, and the downstream `+ 0` are identical across these storage
  # types, so no result changes. Keep assn_raw as the as-loaded ngRMatrix so
  # the gene-QC colSums matches the legacy path exactly.
  assn_raw <- assn
  assn <- as(as(assn, "lMatrix"), "CsparseMatrix")   # ngRMatrix -> lgCMatrix
  if (length(excluded_idx) > 0) {
    assn[, excluded_idx] <- FALSE
  }
  usable_cells <- setdiff(seq_len(ncol(assn)), excluded_idx)

  cat("load_assignment: ", ncol(assn), " cells; ", length(excluded_idx),
      " excluded; ", length(usable_cells), " usable.\n", sep = "")

  list(
    assn         = assn,        # excluded cells zeroed -- use for all selection
    assn_raw     = assn_raw,    # as-loaded -- use only for legacy-matching gene QC
    excluded_idx = excluded_idx,
    usable_cells = usable_cells,
    metadata     = meta
  )
}


# ---- gene QC --------------------------------------------------------------
#' Genes passing QC. Lifted verbatim from the legacy drivers (behaviour
#' unchanged): a threshold on the expected number of perturbed cells, for a
#' typical target, in which the gene is nonzero. NT IS included in the guide /
#' target counts (only nt_off_target and unknown are dropped), exactly as before.
#' moi_observed is taken over the as-loaded matrix (all cells).
#'
#' Because gene_n_nonzero is the only per-gene term, this is algebraically a
#' cutoff on gene_n_nonzero; the printed count lets you calibrate gene_qc_thresh.
#'
#' @param assn_mat  the as-loaded assignment matrix (pass load_assignment()$assn_raw)
compute_genes_passing_qc <- function(gene_summary_stats, assn_mat, grna_target_df,
                                     gene_qc_thresh,
                                     ignore_targets = c("nt_off_target", "unknown")) {
  moi_observed <- mean(Matrix::colSums(assn_mat))

  keep <- !grna_target_df$grna_target %in% ignore_targets
  total_num_guides <- length(unique(grna_target_df$grna_id[keep]))
  avg_num_guides_per_target <- grna_target_df[keep, , drop = FALSE] |>
    dplyr::group_by(grna_target) |>
    dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
    dplyr::pull(n) |>
    mean()

  qc_stat <- gene_summary_stats$gene_n_nonzero * moi_observed /
    total_num_guides * avg_num_guides_per_target
  genes <- gene_summary_stats$gene[qc_stat >= gene_qc_thresh]

  cat(length(genes), " genes (out of ", nrow(gene_summary_stats),
      ") pass gene QC (thresh=", gene_qc_thresh, ", moi=", round(moi_observed, 3), ").\n", sep = "")
  genes
}


# ===========================================================================
# Below: ported UNCHANGED from legacy utils_data_prep.R (output contract).
# ===========================================================================

# ---- odm_to_sparse_matrix -------------------------------------------------
#' Convert OnDisk Matrix subset to sparse matrix (genes x cells).
odm_to_sparse_matrix <- function(odm, genes, cell_idx, set_rownames = TRUE) {
  ilist <- jlist <- xlist <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    curr_umis <- odm[genes[i], ][cell_idx]
    is_positive <- curr_umis > 0
    num_entries <- sum(is_positive)
    if (num_entries > 0) {
      ilist[[i]] <- rep(i, num_entries)
      jlist[[i]] <- which(is_positive)
      xlist[[i]] <- curr_umis[is_positive]
    }
  }

  result <- Matrix::sparseMatrix(
    i = unlist(ilist),
    j = unlist(jlist),
    x = unlist(xlist),
    dims = c(length(genes), length(cell_idx))
  )

  if (set_rownames) {
    rownames(result) <- genes
  }

  return(result)
}


# ---- prepare_frperturb_covariates -----------------------------------------
#' Prepare covariates for FR-Perturb format (log1p specified cols + perturbation).
prepare_frperturb_covariates <- function(cell_covariates, grna_targets, covariates_to_log1p) {
  stopifnot(all(covariates_to_log1p %in% names(cell_covariates)))
  stopifnot(nrow(cell_covariates) == length(grna_targets))

  cell_covs_frpert <- cell_covariates
  for (col in covariates_to_log1p) {
    new_col_name <- paste0(col, "_log1p")
    cell_covs_frpert[[new_col_name]] <- log1p(cell_covariates[[col]])
  }

  cell_covs_frpert <- cell_covs_frpert |>
    dplyr::select(-dplyr::all_of(covariates_to_log1p)) |>
    dplyr::mutate(perturbation = grna_targets)

  stopifnot(!any(is.na(cell_covs_frpert)))

  return(cell_covs_frpert)
}


# ---- setup_anndata_environment --------------------------------------------
#' Set up conda environment for AnnData / FR-Perturb h5ad writing.
setup_anndata_environment <- function(env_name = "r-anndata") {
  library(reticulate)
  curr_envs <- conda_list()$name
  if (!env_name %in% curr_envs) {
    cat("Creating conda environment:", env_name, "\n")
    conda_create(env_name, packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"), envname = env_name, pip = TRUE)
  }
  use_condaenv(env_name, required = TRUE)
  py_config()
  return(TRUE)
}


# ---- write_sceptre_output -------------------------------------------------
write_sceptre_output <- function(response_matrix, grna_matrix, cell_covariates,
                                 grna_target_df, discovery_pairs, formula_object,
                                 output_path) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  saveRDS(response_matrix, file.path(output_path, "response_matrix.rds"))
  saveRDS(grna_matrix, file.path(output_path, "grna_matrix.rds"))
  write.csv(cell_covariates, file.path(output_path, "cell_covariates.csv"), row.names = FALSE)
  write.csv(grna_target_df, file.path(output_path, "grna_target_data_frame.csv"), row.names = FALSE)
  saveRDS(formula_object, file.path(output_path, "formula_object.rds"))
  if (!is.null(discovery_pairs)) {
    saveRDS(discovery_pairs, file.path(output_path, "discovery_pairs.rds"))
  }

  cat("   SCEPTRE format written to:", output_path, "\n")
  if (!is.null(discovery_pairs)) {
    cat("   Discovery pairs:", nrow(discovery_pairs), "\n")
  }
}


# ---- write_mixscale_output ------------------------------------------------
write_mixscale_output <- function(response_matrix, cell_info, output_path, all_targets = NULL) {
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  cell_names <- paste0("cell_idx_", cell_info$cell_idx)

  response_mat_mixscale <- response_matrix |> `colnames<-`(cell_names)
  assign_vec <- cell_info$grna_target |> setNames(cell_names)

  saveRDS(response_mat_mixscale, file.path(output_path, "response_matrix.rds"))
  saveRDS(assign_vec, file.path(output_path, "assignments.rds"))
  if (!is.null(all_targets)) {
    saveRDS(all_targets, file.path(output_path, "mixscale_nt_targets.rds"))
  }

  cat("   Mixscale format written to:", output_path, "\n")
  if (!is.null(all_targets)) {
    cat("   Mixscale targets:", length(all_targets), "\n")
  }
}


# ---- write_frperturb_output -----------------------------------------------
write_frperturb_output <- function(response_matrix, cell_names, cell_covariates_frpert,
                                   output_path, perturbation_in_covariates = TRUE) {
  if (perturbation_in_covariates && !"perturbation" %in% colnames(cell_covariates_frpert)) {
    stop('"perturbation" is expected but not found as a colname of `cell_covariates_frpert`.')
  }
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  setup_anndata_environment()

  library(SingleCellExperiment)
  library(zellkonverter)

  stopifnot(length(cell_names) == ncol(response_matrix))
  stopifnot(nrow(cell_covariates_frpert) == ncol(response_matrix))

  response_subset_frpert <- response_matrix |> `colnames<-`(cell_names)

  sce <- SingleCellExperiment(
    assays  = list(counts = response_subset_frpert),
    colData = cell_covariates_frpert
  )

  writeH5AD(sce, file = file.path(output_path, "response_matrix.h5ad"), X_name = "counts")

  cat("   FR-Perturb format written to:", output_path, "\n")
}
