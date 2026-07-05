# lib/assemble.R
# The ONE dataset-construction core. Every builder (neg / comp / pos, both
# regimes) runs through build_dataset(): it takes a DatasetSpec (what to
# include) + a dataset config (how) + the source bundle, and does steps 4-11:
# extract -> build grna indicator -> enforce -> filter -> metadata ->
# discovery pairs -> validate -> write. filter_all_objects() and
# validate_dataset() run UNCONDITIONALLY, so the legacy drift/bugs are
# structurally impossible.


#' Build one benchmark dataset end to end and write all requested formats.
#'
#' @param spec         DatasetSpec from a select_*() function
#' @param dataset      dataset config (dataset_gasperini / dataset_replogle)
#' @param src          bundle from read_source_data() (needs response_odm, cell_covariates)
#' @param assn         zeroed assignment matrix from load_assignment()$assn
#' @param excluded_idx integer excluded-cell indices (for the leak assertion)
#' @param dataset_name output subdirectory name
#' @param output_root  base dir for this dataset TYPE (e.g. .../neg-control/input_data)
#' @param random_seed  seed for enforce_single_guide (no-op under "maximum")
#' @param concat_string FR-Perturb multi-target delimiter (high-MOI only)
#' @param methods_to_skip character subset of c("sceptre","mixscale","frperturb")
build_dataset <- function(spec, dataset, src, assn, excluded_idx,
                          dataset_name, output_root, random_seed,
                          concat_string = "@", methods_to_skip = character(0)) {

  response_odm    <- src$response_odm
  cell_covariates <- src$cell_covariates
  moi             <- dataset$moi

  # --- 4. extract response (no drop yet) -----------------------------------
  cat("\n=== Extracting response matrix ===\n")
  response_subset <- odm_to_sparse_matrix(response_odm, spec$genes, spec$candidate_cells)
  cat("Extracted:", nrow(response_subset), "genes x", ncol(response_subset), "cells (unfiltered)\n")

  # --- 5-6. grna indicator on the same candidate cells ---------------------
  grna_indicator <- assn[spec$guides, spec$candidate_cells, drop = FALSE] + 0   # -> dgCMatrix

  # --- 7. enforce single guide (low-MOI; no-op under "maximum") ------------
  if (moi$enforce_single_guide) {
    grna_indicator <- enforce_single_guide_per_cell(grna_indicator, random_seed)
  }

  # --- 8. consistency filter (always) --------------------------------------
  cat("\n=== Filtering all objects ===\n")
  filtered <- filter_all_objects(
    response_matrix = response_subset,
    grna_indicator  = grna_indicator,
    grna_target_df  = spec$grna_target_df,
    cell_idx        = spec$candidate_cells,
    low_moi         = moi$low_moi
  )
  response_subset <- filtered$response_matrix
  grna_indicator  <- filtered$grna_indicator
  grna_target_df  <- filtered$grna_target_df
  cell_idx        <- filtered$cell_idx
  genes           <- filtered$genes
  cat("After filtering:", nrow(response_subset), "genes x", ncol(response_subset),
      "cells x", nrow(grna_indicator), "guides\n")

  # --- 9. cell metadata + covariates ---------------------------------------
  cell_info <- moi$build_cell_info(grna_indicator, grna_target_df, cell_covariates,
                                   cell_idx, concat_string)
  cell_covariates_subset <- subset_cell_covariates(cell_covariates, cell_idx, dataset$batch_col)

  # --- 10. discovery pairs (surviving discovery targets x final genes) -----
  if (!is.null(spec$discovery_targets)) {
    surviving_targets <- intersect(spec$discovery_targets, grna_target_df$grna_target)
    discovery_pairs <- expand.grid(
      grna_target = surviving_targets,
      response_id = genes,
      stringsAsFactors = FALSE
    )
  } else {
    discovery_pairs <- NULL
  }

  # --- 11a. validate -------------------------------------------------------
  validate_dataset(
    response_matrix = response_subset, grna_indicator = grna_indicator,
    grna_target_df = grna_target_df, cell_info = cell_info,
    cell_covariates_subset = cell_covariates_subset, cell_idx = cell_idx,
    discovery_pairs = discovery_pairs, dataset = dataset, excluded_idx = excluded_idx
  )

  # --- 11b. write ----------------------------------------------------------
  write_dataset_outputs(
    write_fp = file.path(output_root, dataset_name),
    response_subset = response_subset, grna_indicator = grna_indicator,
    grna_target_df = grna_target_df, discovery_pairs = discovery_pairs,
    cell_info = cell_info, cell_covariates_subset = cell_covariates_subset,
    dataset = dataset, spec = spec, concat_string = concat_string,
    methods_to_skip = methods_to_skip
  )

  cat("\nDataset creation complete! Output:", file.path(output_root, dataset_name), "\n")
  invisible(list(
    path          = file.path(output_root, dataset_name),
    cell_idx      = cell_idx,          # final kept cells (into the original ODM)
    n_cells_total = ncol(assn)         # full universe size (for cells_to_remove)
  ))
}


#' LIGHTWEIGHT sceptre-pipeline input generator. Writes ONLY the two files the
#' full-ODM pipeline consumes, in the dirs their consumers expect:
#'   sceptre/discovery_pairs.rds          (read by set_analysis_parameters.R --
#'                                         same file the direct sceptre run uses)
#'   sceptre-pipeline/cells_to_remove.csv (read by run-sceptre-pipeline.sh)
#' It does NOT form any response/grna matrix, call odm_to_sparse_matrix, or touch
#' the ODM. Feed it a spec from select_comp() called with a CHARACTER gene vector:
#' the selection reads only the in-memory assignment matrix, so nothing heavy runs.
#'
#' cells_to_remove is the complement of the (pre-filter) candidate cells; the
#' pipeline's own permissive run_qc handles any degenerate cells, so we skip the
#' all-zero-response filter (which would require an ODM read).
write_pipeline_inputs <- function(spec, n_cells_total, dataset_name, output_root) {
  stopifnot(!is.null(spec$discovery_targets))
  ds_dir <- file.path(output_root, dataset_name)

  scep_dir <- file.path(ds_dir, "sceptre")
  dir.create(scep_dir, showWarnings = FALSE, recursive = TRUE)
  discovery_pairs <- expand.grid(grna_target = spec$discovery_targets,
                                 response_id = spec$genes, stringsAsFactors = FALSE)
  saveRDS(discovery_pairs, file.path(scep_dir, "discovery_pairs.rds"))

  write_cells_to_remove(ds_dir, spec$candidate_cells, n_cells_total)   # -> sceptre-pipeline/

  cat("pipeline inputs (no matrix):", nrow(discovery_pairs), "pairs (npairs);",
      length(spec$candidate_cells), "kept cells ->", ds_dir, "\n")
  invisible(ds_dir)
}


#' Write the sceptre-pipeline cells_to_remove.csv = the complement of the kept
#' cells within the FULL cell universe. Globally-excluded (junk) cells are
#' automatically included (they're never among the kept cells). Enforces the
#' partition invariant: kept cells and removed cells tile the universe exactly.
write_cells_to_remove <- function(dataset_dir, kept_cell_idx, n_cells_total) {
  stopifnot(!any(duplicated(kept_cell_idx)))
  stopifnot(max(kept_cell_idx) <= n_cells_total)
  pipe_dir <- file.path(dataset_dir, "sceptre-pipeline")
  dir.create(pipe_dir, showWarnings = FALSE, recursive = TRUE)
  cells_to_remove <- setdiff(seq_len(n_cells_total), kept_cell_idx)
  stopifnot(length(cells_to_remove) + length(kept_cell_idx) == n_cells_total)  # partition
  writeLines(as.character(cells_to_remove), file.path(pipe_dir, "cells_to_remove.csv"))
  cat("   cells_to_remove.csv:", length(cells_to_remove), "removed /", n_cells_total, "total\n")
}


#' Write the shared metadata files + each requested format (sceptre / mixscale /
#' frperturb). Writer selection: sceptre + frperturb always; mixscale only if the
#' MOI regime uses it. methods_to_skip removes any of them.
write_dataset_outputs <- function(write_fp, response_subset, grna_indicator,
                                  grna_target_df, discovery_pairs, cell_info,
                                  cell_covariates_subset, dataset, spec,
                                  concat_string, methods_to_skip) {
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)

  # shared, top-level
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)

  writers <- c("sceptre", "frperturb")
  if (dataset$moi$writes_mixscale) writers <- c(writers, "mixscale")
  writers <- setdiff(writers, methods_to_skip)

  # SCEPTRE
  if ("sceptre" %in% writers) {
    write_sceptre_output(
      response_matrix = response_subset, grna_matrix = grna_indicator,
      cell_covariates = cell_covariates_subset, grna_target_df = grna_target_df,
      discovery_pairs = discovery_pairs, formula_object = dataset_formula(dataset),
      output_path = file.path(write_fp, "sceptre")
    )
  }

  # Mixscale (low-MOI only; needs cell_info$grna_target). The mixscale wrapper
  # derives its target loop from unique(assignments), so we do NOT write a
  # separate mixscale_nt_targets.rds (it was redundant with assignments.rds).
  if ("mixscale" %in% writers) {
    write_mixscale_output(
      response_matrix = response_subset, cell_info = cell_info,
      output_path = file.path(write_fp, "mixscale")
    )
  }

  # FR-Perturb
  if ("frperturb" %in% writers) {
    covariates_to_log1p <- setdiff(names(cell_covariates_subset), dataset$batch_col)
    cell_covariates_frpert <- prepare_frperturb_covariates(
      cell_covariates = cell_covariates_subset,
      grna_targets = dataset$moi$perturbation_of(cell_info),
      covariates_to_log1p = covariates_to_log1p
    )
    write_frperturb_output(
      response_matrix = response_subset,
      cell_names = dataset$frperturb_cell_names(cell_info),
      cell_covariates_frpert = cell_covariates_frpert,
      output_path = file.path(write_fp, "frperturb")
    )
  } else {
    cat("Skipping frperturb.\n")
  }
}
