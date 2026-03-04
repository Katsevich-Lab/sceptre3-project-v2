

make_computational_replogle <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    num_genes,
    num_targets,
    max_num_cells,
    nt_name = "non-targeting",
    random_seed = 243535
) {
  # Negative control data preparation for Replogle (LOW MOI):
  # 1. Select genes and targets randomly
  # 2. Get all targeting and non-targeting guides, get cells, downsample
  # 3. Get response matrix, remove all-0 cells, update objects
  # 4. make grna indicator matrix, enforce low moi
  # 5. get metadata, make discovery pairs
  # 6. Write outputs in SCEPTRE, Mixscale, and FR-Perturb formats
  
  # Section 0: Setup
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")
  
  # Section 1: Gene and target selection
  genes <- select_genes_random(response_odm, num_genes, random_seed)
  
  targets <- select_targets_random(
    grna_target_df, num_targets,
    exclude = c(nt_name, "nt_off_target", "unknown"),
    random_seed = random_seed
  )
  
  # 2. get all guides, and get the cells for these guides with possible downsampling
  guides = grna_target_df$grna_id[grna_target_df$grna_target %in% c(targets, nt_name)]
  
  candidate_cells = which(Matrix::colSums(scep_assn_mat[guides,]) > 0)
  cat("Initially", length(candidate_cells), "cells found.\n")
  if(length(candidate_cells) >= max_num_cells) {
    cat("Downsampling to", max_num_cells, "cells.\n")
    candidate_cells <- sample(candidate_cells, max_num_cells)
  } else {
    # cat("Fewer than", max_num_cells, "cells; randomly sampling more.\n")
    # all_other_cells = setdiff(1:ncol(scep_assn_mat), candidate_cells)
    # candidate_cells <- c(candidate_cells, sample(all_other_cells, max_num_cells - length(candidate_cells)))
    cat("Fewer than", max_num_cells, "cells; all kept.\n")
  }
  
  
  # 3. get response matrix, remove cells with all-0 expression, update cells
  response_subset = odm_to_sparse_matrix(response_odm, genes, candidate_cells)
  non_zero_cells = Matrix::colSums(response_subset) > 0
  response_subset = response_subset[, non_zero_cells]
  cell_idx = candidate_cells[non_zero_cells]
  cat("Response matrix created.", length(cell_idx), "cells remain after removing all-0 expression cells.\n")
  
  # 4. make grna_matrix and enforce low MOI
  grna_indicator = scep_assn_mat[guides, cell_idx] |>
    enforce_single_guide_per_cell(random_seed = random_seed)
  cat("grna indicator matrix made and low moi enforced.\n")
  
  # 5. prepare metadata and discovery pairs
  metadata = prepare_cell_metadata_low_moi(grna_indicator, grna_target_df, cell_covariates, cell_idx)
  cat("Cell metadata created.\n")
  
  cell_info <- metadata$cell_info
  cell_covariates_subset = metadata$cell_covariates_subset |>
    dplyr::select(-batch)
  
  
  grna_target_df_subset = grna_target_df[grna_target_df$grna_id %in% rownames(grna_indicator), ]
  
  
  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))
  
  # Create discovery pairs (Cartesian product)
  discovery_pairs <- expand.grid(
    grna_target = targets,
    response_id = genes,
    stringsAsFactors = FALSE
  )
  
  
  # Section 6: Write outputs
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/computational/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)
  
  # Write shared metadata files
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")
  
  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"),
            row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")
  
  # Formula object
  formula_string <- "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1)"
  # if(include_batch) {
  #   formula_string <- paste(formula_string, "+ batch")
  # } else {
  #   # it was still included when they were saved to disc,
  #   # but not included in any of the model specific datasets
  #   cell_covariates_subset$batch <- NULL
  # }
  
  # Write SCEPTRE format
  write_sceptre_output(
    response_matrix = response_subset,
    grna_matrix = grna_indicator,
    cell_covariates = cell_covariates_subset,
    grna_target_df = grna_target_df_subset,
    discovery_pairs = discovery_pairs,
    formula_object = formula_string,
    output_path = file.path(write_fp, "sceptre")
  )
  
  # Write Mixscale format
  write_mixscale_output(
    response_matrix = response_subset,
    cell_info = cell_info,
    output_path = file.path(write_fp, "mixscale")
  )
  
  
  # Prepare FR-Perturb covariates (Replogle: log1p all covariates)
  # regardless of whether batch is used, do not log it
  covariates_to_log1p <- setdiff(names(cell_covariates_subset), "batch")
  cell_covariates_frpert_with_perturbation <- prepare_frperturb_covariates(
    cell_covariates = cell_covariates_subset,
    grna_targets = cell_info$grna_target,
    covariates_to_log1p = covariates_to_log1p
  )
  
  write_frperturb_output(
    response_matrix = response_subset,
    cell_names =  paste0("cell_idx_", cell_info$cell_idx),
    cell_covariates_frpert = cell_covariates_frpert_with_perturbation,
    output_path = file.path(write_fp, "frperturb")
  )
  
  cat("\n Computational benchmarking dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")
}