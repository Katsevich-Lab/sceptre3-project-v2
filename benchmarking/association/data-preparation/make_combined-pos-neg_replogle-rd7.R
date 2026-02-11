library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/replogle-rd7/sceptre-pipeline"
)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/replogle-rd7/sceptre-pipeline"
)

# Helper function to convert ODM to sparse matrix
odm_to_sparse_matrix <- function(odm, genes, cell_idx) {
  ilist <- jlist <- xlist <- vector("list", length(genes))
  for(i in seq_along(genes)) {
    curr_umis <- odm[genes[i], ][cell_idx]
    is_positive <- curr_umis > 0
    num_entries <- sum(is_positive)
    if(num_entries > 0) {
      ilist[[i]] <- rep(i, num_entries)
      jlist[[i]] <- which(is_positive)
      xlist[[i]] <- curr_umis[is_positive]
    }
  }
  Matrix::sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
                       dims = c(length(genes), length(cell_idx))) |>
    `rownames<-`(genes)
}


# for stepping thru the function
# num_genes = 10
# 
# num_pos_genes = num_genes
# min_cells_per_target = 50
# num_nt_guides = 20
# num_extra_genes = 10
# pos_and_neg_pairs_only = TRUE
# dataset_name = paste0("replogle-rd7_discrimination_", num_genes, "genes")
# response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
# grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
# cell_covariates = scep@covariate_data_frame
# scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
# grna_target_df = scep@grna_target_data_frame
# 
# pos_and_neg_pairs_only = TRUE
# random_seed = 12345

make_discrimination_replogle <- function(
  num_pos_genes,
  min_cells_per_target,
  num_nt_guides,
  num_extra_genes,
  pos_and_neg_pairs_only = TRUE,
  dataset_name,
  response_odm,
  grna_odm,
  cell_covariates,
  scep_assn_mat,
  grna_target_df,
  num_mixscale_nts,
  random_seed = 12345
) {
  # Discrimination dataset for Replogle (LOW MOI):
  # Combines positive control (cells with targeting guides) and
  # negative control (cells with NT guides only) in a single dataset
  # for discrimination analysis

  # Section 0: Setup and random sampling -------------------
  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")

  # Identify targeting guides and NT guides
  non_targeting_target_labels <- c("non-targeting", "nt_off_target", "unknown")
  stopifnot(non_targeting_target_labels %in% grna_target_df$grna_target) # all 3 should be there
  targeting_guides <- grna_target_df$grna_id[!grna_target_df$grna_target %in% non_targeting_target_labels]
  nt_guides <- grna_target_df$grna_id[grna_target_df$grna_target == "non-targeting"]

  cat("Found", length(targeting_guides), "targeting guides\n")
  cat("Found", length(nt_guides), "NT guides\n")


  # Section 1: Select positive control genes (random) ---------------
  # Get genes that have targeting guides (intersection)
  on_target_genes <- setdiff(grna_target_df$grna_target, non_targeting_target_labels)
  all_response_genes <- rownames(response_odm)
  candidate_genes <- intersect(on_target_genes, all_response_genes)

  cat("Found", length(candidate_genes), "genes with targeting guides\n")

  # Shuffle candidate genes for randomness
  candidate_genes_shuffled <- sample(candidate_genes)

  # Iterate through genes, keeping those with sufficient coverage
  selected_pos_genes <- character()
  targeting_guides_to_keep <- character()
  # gene = "ENSG00000131037" # has 5 guides
  for(gene in candidate_genes_shuffled) {
    # Stop if we've reached the target number
    if(length(selected_pos_genes) >= num_pos_genes) {
      break
    }

    # Get all guides targeting this gene
    gene_guides <- grna_target_df$grna_id[grna_target_df$grna_target == gene]

    # Check if every guide has sufficient cells
    has_sufficient_coverage <- any(
      Matrix::rowSums(scep_assn_mat[gene_guides,,drop=FALSE]) >= min_cells_per_target
    )

    if(has_sufficient_coverage) {
      selected_pos_genes <- c(selected_pos_genes, gene)
      targeting_guides_to_keep <- union(targeting_guides_to_keep, gene_guides)
    }
  }
  # ordering them the same way they appear as columns
  # selected_pos_genes <- all_response_genes[all_response_genes %in% selected_pos_genes]
  cat("Selected", length(selected_pos_genes), "positive control genes with sufficient coverage\n")
  cat("Keeping", length(targeting_guides_to_keep), "targeting guides\n")

  if(length(selected_pos_genes) == 0) {
    stop("No genes met the minimum cell coverage requirement")
  }

  # Section 2: Collect positive control cells ------------------
  # any cell that has an on-target guide is ok, even if it has others too.
  candidate_pos_cell_idx <- which(Matrix::colSums(scep_assn_mat[targeting_guides_to_keep, ]) >= 1)

  cat("Found", length(candidate_pos_cell_idx), "candidate positive control cells\n")

  # Section 3: Select and collect negative control cells --------
  # Sample NT guides
  if(num_nt_guides > length(nt_guides)) {
    warning("Requested more NT guides than available, using all")
    selected_nt_guides <- nt_guides
  } else {
    selected_nt_guides <- sample(nt_guides, num_nt_guides)
  }

  cat("Selected", length(selected_nt_guides), "NT guides\n")

  # Filter cells: NO targeting guides + at least one selected NT guide
  cells_no_targeting <- which(Matrix::colSums(scep_assn_mat[targeting_guides, ]) == 0)
  cells_with_nt <- which(Matrix::colSums(scep_assn_mat[selected_nt_guides, ]) > 0)
  candidate_neg_cell_idx <- intersect(cells_no_targeting, cells_with_nt)

  cat("Found", length(candidate_neg_cell_idx), "candidate negative control cells\n")

  # Section 4: Combine cells and filter -----------------------
  # Combine cell indices
  candidate_all_cell_idx <- c(candidate_pos_cell_idx, candidate_neg_cell_idx)
  candidate_cell_status <- rep(c("on_target", "negative_control"), c(length(candidate_pos_cell_idx), length(candidate_neg_cell_idx)))
  # Extract response matrix for selected genes
  response_subset <- odm_to_sparse_matrix(response_odm, selected_pos_genes, candidate_all_cell_idx)

  # Remove cells with zero UMI
  cells_with_umis <- Matrix::colSums(response_subset) > 0
  response_subset <- response_subset[, cells_with_umis]
  all_cell_idx <- candidate_all_cell_idx[cells_with_umis]
  cell_status <- candidate_cell_status[cells_with_umis]
  # remove genes and corresponding guides with no expression too
  genes_with_umis <- Matrix::rowSums(response_subset) > 0
  response_subset <- response_subset[genes_with_umis, ]
  selected_pos_genes <- rownames(response_subset)  
  targeting_guides_to_keep <- grna_target_df$grna_id[grna_target_df$grna_target %in% selected_pos_genes]
  
  cat("Removed", sum(!cells_with_umis), "cells with zero UMI\n")
  cat("Removed", sum(!genes_with_umis), "genes with zero UMI\n")

  # Add extra genes for additional negative control pairs
  if(num_extra_genes > 0) {
    all_extra_genes <- setdiff(all_response_genes, selected_pos_genes)
    if(num_extra_genes > length(all_extra_genes)) {
      warning("Requested more extra genes than available, using all")
      selected_extra_genes <- all_extra_genes
    } else {
      selected_extra_genes <- sample(all_extra_genes, num_extra_genes)
    }

    # cat("Selected", length(selected_extra_genes), "extra genes for negative control pairs\n")

    # Extract response matrix for extra genes
    response_extra <- odm_to_sparse_matrix(response_odm, selected_extra_genes, all_cell_idx)
    # even though most n_trt == 0 still, we'll require at least one cell
    response_extra <- response_extra[Matrix::rowSums(response_extra) > 0, , drop=FALSE]
    selected_extra_genes <- rownames(response_extra)
    cat("Selected", length(selected_extra_genes), "extra genes for negative control pairs\n")
    # Combine with main response matrix
    response_subset <- rbind(response_subset, response_extra)

    cat("Total genes in dataset:", nrow(response_subset),
        "(", length(selected_pos_genes), "pos +", length(selected_extra_genes), "extra)\n")
  }

  all_guides_to_keep <- c(targeting_guides_to_keep, selected_nt_guides)
  # now ordering them as they appear
  # all_guides_to_keep <- rownames(scep_assn_mat)[rownames(scep_assn_mat) %in% all_guides_to_keep]
  dummy_grna_indicator <- scep_assn_mat[all_guides_to_keep, all_cell_idx] + 0L

  # x <- dummy_grna_indicator |> as.matrix()
  # image(
  #   t(x[nrow(x):1, ]),            # flip rows so x[1,1] ends up top-left
  #   col = c("white", "black"),
  #   axes = FALSE,
  #   asp = 1
  # )
  # enforcing low moi
  for(cell_idx in seq_len(ncol(dummy_grna_indicator))) {
    expressed_guides <- which(dummy_grna_indicator[, cell_idx] == 1)
    if(length(expressed_guides) > 1) {
      keep_guide <- sample(expressed_guides, 1)
      dummy_grna_indicator[setdiff(expressed_guides, keep_guide), cell_idx] <- 0
    }
  }

  # Verify enforcement
  stopifnot(all(Matrix::colSums(dummy_grna_indicator) == 1))
  stopifnot(all(Matrix::rowSums(dummy_grna_indicator) >= 1))
  
  # Section 5: Create unified gRNA target mapping --------------
  # Create pseudo-targets for NT guides
  nt_pseudo_targets <- data.frame(
    grna_id = selected_nt_guides,
    grna_target = paste0("non-targeting", seq_along(selected_nt_guides)),
    stringsAsFactors = FALSE
  )

  # Get real target mappings
  real_target_mappings <- grna_target_df[grna_target_df$grna_id %in% targeting_guides_to_keep, ] |>
    dplyr::select(grna_id, grna_target)
    

  # Combine
  unified_grna_target_df <- rbind(real_target_mappings, nt_pseudo_targets) |>
    as.data.frame()  # it's a data.table

  cat("Unified target mapping:", nrow(unified_grna_target_df), "total guides\n")
  cat("  Real targets:", nrow(real_target_mappings), "\n")
  cat("  Pseudo-targets:", nrow(nt_pseudo_targets), "\n")

  # Verify dimensions
  stopifnot(ncol(dummy_grna_indicator) == ncol(response_subset))
  stopifnot(nrow(dummy_grna_indicator) == nrow(unified_grna_target_df))

  # Section 6: Create metadata and discovery pairs -------------
  # Cell info with cell type labels
  cell_info <- data.frame(
    cell_idx = all_cell_idx,
    cell_type = cell_status,
    cell_name = rownames(cell_covariates)[all_cell_idx],
    grna_id = rownames(dummy_grna_indicator)[apply(dummy_grna_indicator, 2, function(cell) which(cell == 1))],
    stringsAsFactors = FALSE
  ) |>
    left_join(unified_grna_target_df, by = "grna_id")
  
  stopifnot(setequal(
    selected_pos_genes,
    intersect(rownames(response_subset), unified_grna_target_df$grna_target)
  ))
  
            
  # Create discovery pairs based on mode
  if(pos_and_neg_pairs_only) {
    cat("Creating discovery pairs: on-target + negative control only\n")

    # Only on-target pairs + negative control pairs
    on_target_pairs <- data.frame(
      grna_target = selected_pos_genes,
      response_id = selected_pos_genes,
      pair_type = "on_target",
      stringsAsFactors = FALSE
    )

    neg_control_pairs <- expand.grid(
      grna_target = nt_pseudo_targets$grna_target,
      response_id = rownames(response_subset),
      stringsAsFactors = FALSE
    ) |>
      mutate(
        pair_type = ifelse(
          response_id %in% selected_pos_genes,
          "negative_control_on",
          "negative_control_off"
        )
      )
    discovery_pairs <- rbind(on_target_pairs, neg_control_pairs)

  } else {
    stop("not tested yet")
    cat("Creating discovery pairs: full Cartesian product\n")

    # Full Cartesian product
    discovery_pairs <- expand.grid(
      grna_target = unified_grna_target_df$grna_target,
      response_id = rownames(response_subset),
      stringsAsFactors = FALSE
    )

    # Add pair type labels
    discovery_pairs$pair_type <- ifelse(
      discovery_pairs$grna_target %in% nt_pseudo_targets$grna_target,
      "negative_control",
      ifelse(
        discovery_pairs$grna_target == discovery_pairs$response_id,
        "on_target",
        "off_target"
      )
    )
  }

  cat("Discovery pairs created:", nrow(discovery_pairs), "total\n")
  cat("  ", table(discovery_pairs$pair_type), "\n")

  # Section 7: Prepare covariates (NO BATCH for Replogle) ------
  cell_covariates_subset <- cell_covariates[all_cell_idx, ] |>
    transmute(
      response_n_nonzero_full = response_n_nonzero,
      response_n_umis_full = response_n_umis,
      grna_n_nonzero_full = grna_n_nonzero,
      grna_n_umis_full = grna_n_umis
    )

  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))

  # Section 8: Write outputs ------------------------------------
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/combined-pos-neg/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)

  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")

  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")

  ## Section 8a: SCEPTRE format ---------------------------------
  write_sceptre_fp <- file.path(write_fp, "sceptre")
  dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)

  # Formula string (NO BATCH for Replogle)
  formula_string <- "~ 1 + log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + log(grna_n_umis_full + 1)"

  cell_covariates_sceptre <- cell_covariates_subset

  saveRDS(response_subset, file.path(write_sceptre_fp, "response_matrix.rds"))
  saveRDS(dummy_grna_indicator, file.path(write_sceptre_fp, "grna_matrix.rds"))
  write.csv(cell_covariates_sceptre, file.path(write_sceptre_fp, "cell_covariates.csv"), row.names = FALSE)
  write.csv(unified_grna_target_df, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names = FALSE)
  saveRDS(formula_string, file.path(write_sceptre_fp, "formula_object.rds"))
  saveRDS(discovery_pairs, file.path(write_sceptre_fp, "discovery_pairs.rds"))
  cat("   SCEPTRE format written.\n")

  ## Section 8b: Mixscale format --------------------------------
  write_mixscale_fp <- file.path(write_fp, "mixscale")
  dir.create(write_mixscale_fp, showWarnings = FALSE, recursive = TRUE)

  # Prepare response matrix with cell names
  cell_names <- paste0("cell_idx_", cell_info$cell_idx)
  response_mat_mixscale <- response_subset |> `colnames<-`(cell_names)

  # Create assignments vector: named vector mapping cell names to targets
  assign_vec <- cell_info$grna_target |> setNames(cell_names)

  # For Mixscale: use pseudo-targets for looping (neg control portion)
  if(num_mixscale_nts > nrow(nt_pseudo_targets)) {
    warning("Requested more mixscale NTs than available, using all")
    mixscale_nt_targets <- nt_pseudo_targets$grna_target
  } else {
    mixscale_nt_targets <- sample(nt_pseudo_targets$grna_target, num_mixscale_nts)
  }

  saveRDS(response_mat_mixscale, file.path(write_mixscale_fp, "response_matrix.rds"))
  saveRDS(assign_vec, file.path(write_mixscale_fp, "assignments.rds"))
  saveRDS(mixscale_nt_targets, file.path(write_mixscale_fp, "mixscale_nt_targets.rds"))
  cat("   Mixscale format written.\n")

  ## Section 8c: FR-Perturb format ------------------------------
  write_frperturb_fp <- file.path(write_fp, "frperturb")
  dir.create(write_frperturb_fp, showWarnings = FALSE, recursive = TRUE)

  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)

  # Setup conda environment
  env_name <- "r-anndata"
  curr_envs <- conda_list()$name
  if(!env_name %in% curr_envs) {
    conda_create("r-anndata", packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"),
               envname = "r-anndata", pip = TRUE)
  }
  use_condaenv("r-anndata", required = TRUE)
  py_config()

  # Prepare data for FR-Perturb
  response_mat_frpert <- response_subset |> `colnames<-`(cell_names)

  # Build covariate dataframe with log-transformed covariates
  cell_covs_frpert <- cell_covariates_subset |>
    mutate_all(log1p) %>%
    `names<-`(paste0(names(.), "_log1p")) |>
    mutate(perturbation = cell_info$grna_target)

  stopifnot(!any(is.na(cell_covs_frpert)))

  sce <- SingleCellExperiment(
    assays  = list(counts = response_mat_frpert),
    colData = cell_covs_frpert
  )

  writeH5AD(
    sce,
    file = file.path(write_frperturb_fp, "response_matrix.h5ad"),
    X_name = "counts"
  )
  cat("   FR-Perturb format written.\n")

  cat("\n`combined-pos-neg` dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")
}

# Run the function with test configuration
# num_genes = 10
# make_discrimination_replogle(
#   num_pos_genes = num_genes,
#   min_cells_per_target = 50,
#   num_extra_genes = 10,
#   num_nt_guides = 20,
#   pos_and_neg_pairs_only = TRUE,
#   dataset_name = paste0("replogle-rd7_combined-pos-neg_", num_genes, "genes"),
#   response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
#   grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
#   cell_covariates = scep@covariate_data_frame,
#   scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
#   grna_target_df = scep@grna_target_data_frame
# )
# 
# num_genes = 1000
make_discrimination_replogle(
  num_pos_genes = 100,
  min_cells_per_target = 10,
  num_nt_guides = 113,
  num_extra_genes = 1500,
  pos_and_neg_pairs_only = TRUE,
  dataset_name = "replogle-rd7_combined-pos-neg_medium",
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = readRDS(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame,
  num_mixscale_nts = 3
)
