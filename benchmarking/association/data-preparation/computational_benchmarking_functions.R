
# select_genes_more_expressed <- function(response_odm, num_genes, B = 5000, random_seed = NULL) {
#   if (!is.null(random_seed)) set.seed(random_seed)
#   
#   cat("Sampling genes: estimating gene total umi distribution...")
#   # Pilot sample of row sums (CDF reference)
#   locs <- sample(nrow(response_odm), B)
#   tots <- numeric(B)
#   for (i in seq_len(B)) {
#     tots[i] <- sum(response_odm[locs[i], ])
#   }
#   cat("done.\nSampling genes...")
#   
#   genes <- character(0)
#   all_genes <- sample(rownames(response_odm))  # shuffle
#   
#   for (curr_gene in all_genes) {
#     curr_total <- sum(response_odm[curr_gene, ])
#     prob_keep <- mean(tots <= curr_total)              # estimated percentile
#     if (runif(1) <= prob_keep) {
#       genes <- c(genes, curr_gene)
#       if (length(genes) >= num_genes) break
#     }
#   }
#   cat("done.", num_genes, "found.\n")
#   
#   if (length(genes) < num_genes) {
#     stop("in select_genes_more_expressed() failed to sample enough genes.")
#   }
#   genes
# }



sample_genes_with_expression <- function(response_odm, num_genes, cell_idx, random_seed = NULL) {
  if (!is.null(random_seed)) set.seed(random_seed)
  
  genes <- character(0)
  all_genes <- sample(rownames(response_odm))  # shuffle
  
  for (curr_gene in all_genes) {
    has_umis <- sum(response_odm[curr_gene, ][cell_idx]) > 0
    if(has_umis) {
      genes <- c(genes, curr_gene)
      if (length(genes) >= num_genes) break
    }
  }
  genes
}


prepare_cell_metadata_high_moi <- function(grna_indicator_matrix, grna_target_df, cell_covariates, cell_idx,
                                           fr_perturb_concat_string) {
  # Create cell_info with guide assignments
  expressed_guide_concat = expressed_target_concat = character(ncol(grna_indicator_matrix))
  
  num_assigned_grnas = Matrix::colSums(grna_indicator_matrix)
  if(any(num_assigned_grnas == 0)) {
    stop("In `prepare_cell_metadata_high_moi()`, some cells have no assigned guides, which should not happen here.")
  }
  for(i in seq_len(ncol(grna_indicator_matrix))) {
    curr_guides = rownames(grna_indicator_matrix)[grna_indicator_matrix[,i] == 1]
    expressed_guide_concat[i] = paste0(
      curr_guides,
      collapse=fr_perturb_concat_string
    )
    expressed_target_concat[i] = paste0(
      grna_target_df$grna_target[grna_target_df$grna_id %in% curr_guides],
      collapse=fr_perturb_concat_string
    )
    
  }
    
  
  cell_info <- data.frame(
    cell_idx = cell_idx,
    expressed_guide_concat,
    expressed_target_concat,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(cell_name = rownames(cell_covariates)[cell_idx])
  
  # Subset and prepare covariates (NO BATCH for Replogle but YES batch for gasperini)
  cell_covariates_subset <- cell_covariates[cell_idx, ] |>
    dplyr::transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full = response_n_umis,
      grna_n_nonzero_full = grna_n_nonzero,
      grna_n_umis_full = grna_n_umis,
      prep_batch = prep_batch
    )
  
  return(list(
    cell_info = cell_info,
    cell_covariates_subset = cell_covariates_subset
  ))
}


make_computational_replogle <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    genes,
    num_targets,
    max_num_cells,
    nt_name = "non-targeting",
    random_seed = 243535,
    methods_to_skip=""
) {
  
  # 1. sample targets
  # 2. get all cells for those targets
  # 3. sample genes, and only keep genes that have some umi counts
  # 4. get response matrix
  # 5. drop cells with zero expression [this won't cause any genes to now have zero because these cells were all 0]
  # 6. continue
  
  
  # Section 0: Setup
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")
  
  # 1. get targets
  targets <- select_targets_random(
    grna_target_df, num_targets,
    exclude = c(nt_name, "nt_off_target", "unknown"),
    random_seed = random_seed
  )
  
  # 2. get cells for these targets
  guides = grna_target_df$grna_id[grna_target_df$grna_target %in% c(targets, nt_name)]
  
  candidate_cells = which(Matrix::colSums(scep_assn_mat[guides,]) > 0)
  cat("Initially", length(candidate_cells), "cells found.\n")
  
  if(is.finite(max_num_cells)) {
    if(length(candidate_cells) >= max_num_cells) {
      cat("Downsampling to", max_num_cells, "cells.\n")
      candidate_cells <- sample(candidate_cells, max_num_cells)
    } else {
      # cat("Fewer than", max_num_cells, "cells; randomly sampling more.\n")
      # all_other_cells = setdiff(1:ncol(scep_assn_mat), candidate_cells)
      # candidate_cells <- c(candidate_cells, sample(all_other_cells, max_num_cells - length(candidate_cells)))
      # cat("Fewer than", max_num_cells, "cells; all kept.\n")
      stop("Not enough cells for these parameters.")
    }
  } else {
    cat("max_num_cells = Inf so all cells for these targets are kept.\n")
  }
  
  # 3. use the provided genes, or sample genes and only keep genes that have some umi counts
  if(is.numeric(genes) && length(genes) == 1) {
    num_genes = genes
    genes <- sample_genes_with_expression(response_odm, num_genes, candidate_cells)
  } else if(!is.character(genes)) {
    stop("`genes` must be either a character vector or an integer.\n")
  }
  cat("Using", length(genes), "genes.\n")
  
  # 4. response matrix and drop empty cells
  response_subset = odm_to_sparse_matrix(response_odm, genes, candidate_cells)
  # remove any cells that have all 0 expression
  non_zero_cells = Matrix::colSums(response_subset) > 0
  response_subset = response_subset[, non_zero_cells]
  cell_idx = candidate_cells[non_zero_cells]
  cat("Response matrix created.", length(cell_idx), "cells remain after removing all-0 expression cells.\n")
  
  # 5. make grna_matrix and enforce low MOI
  grna_indicator = (scep_assn_mat[guides, cell_idx] + 0) |> # +0 converts to integer from logical
    enforce_single_guide_per_cell(random_seed = random_seed)
  cat("grna indicator matrix made and low moi enforced.\n")
  
  # 6. prepare metadata and discovery pairs
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
    response_id = rownames(response_subset),
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
  
  if(! "mixscale" %in% methods_to_skip) {
    # Write Mixscale format
    write_mixscale_output(
      response_matrix = response_subset,
      cell_info = cell_info,
      output_path = file.path(write_fp, "mixscale")
    )
  } else {
    cat("Skipping mixscale.\n")
  }
  
  if(! "frperturb" %in% methods_to_skip) {
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
  } else {
    cat("Skipping FR-Perturb.\n")
  }
  

  
  cat("\n Computational benchmarking dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")
  
  
  
  
}




# dataset_name = "gasp-comp-test"
# response_odm
# grna_odm
# cell_covariates
# scep_assn_mat
# grna_target_df
# genes = 10
# num_targets = 12
# max_num_cells = 234
# nt_name = "non-targeting"
# random_seed = 243535



# make_computational_gasperini(
#   dataset_name = "gasp-comp-test",
#   response_odm=response_odm,
#   grna_odm=grna_odm,
#   cell_covariates=cell_covariates,
#   scep_assn_mat=scep_assn_mat,
#   grna_target_df=grna_target_df,
#   genes=10,
#   num_targets=50,
#   max_num_cells=523,
#   nt_name = "non-targeting",
#   force_nt_inclusion = FALSE,
#   fr_perturb_concat_string = "@",
#   random_seed = 243535
# )


make_computational_gasperini <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    genes,
    num_targets,
    max_num_cells,
    nt_name = "non-targeting",
    force_nt_inclusion = FALSE,
    fr_perturb_concat_string = "@",
    random_seed = 243535
) {
  

  if(any(grepl(fr_perturb_concat_string, grna_target_df$grna_target))) {
    stop("`fr_perturb_concat_string` detected in `grna_target_df$grna_target`.")
  }

  
  # Section 0: Setup
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")
  
  # 1. get targets
  targets <- select_targets_random(
    grna_target_df, num_targets,
    exclude = nt_name,
    random_seed = random_seed
  )
  
  # 2. get cells for these targets
  # for high moi, we do not force NT inclusion
  if(force_nt_inclusion) {
    guides = grna_target_df$grna_id[grna_target_df$grna_target %in% c(targets, nt_name)]
  } else {
    guides = grna_target_df$grna_id[grna_target_df$grna_target %in% targets]
  }
  
  
  candidate_cells = which(Matrix::colSums(scep_assn_mat[guides,]) > 0)
  cat("Initially", length(candidate_cells), "cells found.\n")
  if(length(candidate_cells) >= max_num_cells) {
    cat("Downsampling to", max_num_cells, "cells.\n")
    candidate_cells <- sample(candidate_cells, max_num_cells)
  } else {
    cat("Fewer than", max_num_cells, "cells; all kept.\n")
  }
  
  # 3. use the provided genes, or sample genes and only keep genes that have some umi counts
  if(is.numeric(genes) && length(genes) == 1) {
    num_genes = genes
    genes <- sample_genes_with_expression(response_odm, num_genes, candidate_cells)
  } else if(!is.character(genes)) {
    stop("`genes` must be either a character vector or an integer.\n")
  }
  cat("Using", length(genes), "genes.\n")
  
  # 4. response matrix and drop empty cells
  response_subset = odm_to_sparse_matrix(response_odm, genes, candidate_cells)
  # remove any cells that have all 0 expression
  non_zero_cells = Matrix::colSums(response_subset) > 0
  response_subset = response_subset[, non_zero_cells]
  cell_idx = candidate_cells[non_zero_cells]
  cat("Response matrix created.", length(cell_idx), "cells remain after removing all-0 expression cells.\n")
  
  # 5. make grna_matrix and enforce low MOI
  grna_indicator = scep_assn_mat[guides, cell_idx] + 0 # +0 converts to integer from logical
  cat("high MOI grna indicator matrix made.\n")
  
  # 6. prepare metadata and discovery pairs
  metadata = prepare_cell_metadata_high_moi(grna_indicator, grna_target_df, cell_covariates, cell_idx,
                                            fr_perturb_concat_string=fr_perturb_concat_string)
  cat("Cell metadata created.\n")
  
  cell_info <- metadata$cell_info
  cell_covariates_subset = metadata$cell_covariates_subset 
  
  
  grna_target_df_subset = grna_target_df[grna_target_df$grna_id %in% rownames(grna_indicator), ]
  
  
  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))
  
  # Create discovery pairs (Cartesian product)
  discovery_pairs <- expand.grid(
    grna_target = targets,
    response_id = rownames(response_subset),
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
  formula_string <- "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1) + prep_batch"

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
  

  # Prepare FR-Perturb covariates, with all but batch log1p'd
  covariates_to_log1p <- setdiff(names(cell_covariates_subset), "prep_batch")
  cell_covariates_frpert_with_perturbation <- prepare_frperturb_covariates(
    cell_covariates = cell_covariates_subset,
    grna_targets = cell_info$expressed_target_concat,
    covariates_to_log1p = covariates_to_log1p
  )
  
  # with gasperini, we have real cell names
  write_frperturb_output(
    response_matrix = response_subset,
    cell_names =  cell_info$cell_name,
    cell_covariates_frpert = cell_covariates_frpert_with_perturbation,
    output_path = file.path(write_fp, "frperturb")
  )
  
  cat("\n Computational benchmarking dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")

}
