library(tidyverse)

rm(list=ls())
source("~/.Rprofile")

dataset_name <- "gasperini_pos_control"

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/gasperini/sceptre-pipeline"  
)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)

path_to_assigns <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs/gasperini/sceptre-pipeline"  
)

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
                       dims = c(length(genes), length(cell_idx)))
}

make_pos_control_gasperini <- function(dataset_name, response_odm, grna_odm, cell_covariates, scep_assn_mat, grna_target_df) {
  # idea:
  # 1. determine which targets are on-targets
  # 2. determine which cells express these, or NT cells
  # 3. take the subset of response.odm for these on-target genes and the cells expressing them, as well as NT cells.
  # 4. write what each method needs
  
  # 1. getting on-targets
  on_targets <- intersect(grna_target_df$grna_target, rownames(response_odm))
  cat("There are", length(on_targets), "on-target genes.\n")
  
  # 2a. i need to map these on_target genes to their respective guides
  guides_for_on_targets <- grna_target_df$grna_id[grna_target_df$grna_target %in% on_targets]
  cat("There are", length(guides_for_on_targets), "guides for on-target genes.\n")
  
  # 2b. determine which cells express these
  cells_expressing_on_target_guides <- vector("list", length(guides_for_on_targets)) |>
    setNames(guides_for_on_targets)
  for(guide in guides_for_on_targets) {
    cells_expressing_on_target_guides[[guide]] <- which(scep_assn_mat[guide,])
  }
  cat("There are", unlist(cells_expressing_on_target_guides) |> unique() |> length(), "cells targeting these guides.\n")
  
  
  # 2c. add in NT cells
  cells_expressing_nt_guides <- list()
  nt_guides <- grna_target_df$grna_id[grepl("non-targeting", grna_target_df$grna_target)]
  cat("There are", length(nt_guides), "NT guides.\n")
  
  for(nt_guide in nt_guides) {
    cells_expressing_nt_guides[[nt_guide]] <- which(scep_assn_mat[nt_guide,])
  }
  cat("There are", unlist(cells_expressing_nt_guides) |> unique() |> length(), "cells targeting NT guides.\n")
  
  # 2d. collect into long-form cell_info data.frame
  cells_expressing_on_target_or_nt_guides <- c(cells_expressing_on_target_guides, cells_expressing_nt_guides)
  rm(cells_expressing_on_target_guides, cells_expressing_nt_guides)
  cell_info <- lapply(
    names(cells_expressing_on_target_or_nt_guides),
    function(guide) data.frame(cell_id = cells_expressing_on_target_or_nt_guides[[guide]], grna_id = guide)
  ) |>
    do.call(what = rbind) |>
    left_join(grna_target_df |> select(grna_id, grna_target), by = "grna_id") |>
    mutate(cell_name = rownames(cell_covariates)[cell_id])
  
  
  # 3. subset the response odm, save as sparse matrix
  all_cell_idx <- unique(cell_info$cell_id) |> sort()  # only 20k cells trimmed off by this
  
  stopifnot(setequal(on_targets, setdiff(unique(cell_info$grna_target), "non-targeting")))
  
  
  response_subset <- odm_to_sparse_matrix(odm = response_odm, genes = on_targets, cell_idx = all_cell_idx) |>
    `rownames<-`(on_targets)
  cat("response matrix subset made with", nrow(response_subset), "genes and", ncol(response_subset), "cells.\n")
  
  # 4. now i need to save this for sceptre and FR-Perturb
  cell_covariates_subset <- cell_covariates[all_cell_idx, ]
  
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/pos-control/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")
  
  write.csv(cell_covariates_subset, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")
  
  ## 4a. sceptre -----------------------------------------------------
  write_sceptre_fp <- file.path(write_fp, "sceptre")
  dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)
  
  grna_target_df_kept <- grna_target_df |>
    dplyr::filter(grna_id %in% names(cells_expressing_on_target_or_nt_guides))
  
  # Create binary assignment indicator matrix (not UMI counts)
  # grna_indicator_matrix[grna, cell] = 1 if cell is assigned to grna, 0 otherwise
  locs <- lapply(seq_along(cells_expressing_on_target_or_nt_guides), function(i) {
    grna <- names(cells_expressing_on_target_or_nt_guides)[i]
    assigned_cells <- cells_expressing_on_target_or_nt_guides[[grna]]
    list(
      i = rep(i, length(assigned_cells)),
      j = match(assigned_cells, all_cell_idx),
      x = rep(1, length(assigned_cells))
    )
  })
  grna_indicator_matrix <- Matrix::sparseMatrix(
    i = lapply(locs, `[[`, "i") |> unlist(),
    j = lapply(locs, `[[`, "j") |> unlist(),
    x = lapply(locs, `[[`, "x") |> unlist(),
    dims = c(length(cells_expressing_on_target_or_nt_guides), length(all_cell_idx))
  )|>
    `rownames<-`(names(cells_expressing_on_target_or_nt_guides))
  
  # Compute sceptre-specific covariates from the subsetted gRNA UMI data
  # Extract actual UMI counts from grna_odm for this subset of gRNAs and cells
  grna_matrix_subset <- odm_to_sparse_matrix(grna_odm, names(cells_expressing_on_target_or_nt_guides), all_cell_idx)
  
  
  # Compute covariates from the UMI subset (not from binary grna_indicator_matrix)
  grna_n_nonzero_subset <- Matrix::colSums(grna_matrix_subset != 0)
  grna_n_umis_subset <- Matrix::colSums(grna_matrix_subset)
  
  # Create sceptre-specific cell covariates with subset-based gRNA metrics
  cell_covariates_sceptre <- cbind(
    cell_covariates_subset,
    data.frame(
      grna_n_nonzero_subset = grna_n_nonzero_subset,
      grna_n_umis_subset = grna_n_umis_subset
    )
  )
  
  saveRDS(response_subset, file.path(write_sceptre_fp, "response_matrix.rds"))
  saveRDS(grna_indicator_matrix, file.path(write_sceptre_fp, "grna_matrix.rds"))
  write.csv(cell_covariates_sceptre, file.path(write_sceptre_fp, "cell_covariates.csv"), row.names = FALSE)
  write.csv(grna_target_df_kept, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names = FALSE)
  cat("   sceptre written.\n")
  
  ## 4b. FR-Perturb -----------------------------------------------------
  write_frperturb_fp <- file.path(write_fp, "frperturb")
  dir.create(write_frperturb_fp, showWarnings = FALSE, recursive = TRUE)
  
  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)
  env_name <- "r-anndata"
  curr_envs <- conda_list()$name
  if(!env_name %in% curr_envs) {
    conda_create("r-anndata", packages = c("python=3.12"))
    py_install(c("numpy", "scipy", "h5py", "anndata"),
               envname = "r-anndata", pip = TRUE)
  }
  use_condaenv("r-anndata", required = TRUE)
  py_config()
  
  # actually writing
  # see here for input specifications: https://github.com/douglasyao/FR-Perturb 
  # using Option 2 from here
  # this is high MOI so we need to concat to make the perturbation column
  cell_names <- rownames(cell_covariates_subset)
  response_subset_frpert <- response_subset |> `colnames<-`(cell_names)
  
  # these get added to .obs of the anndata object
  cell_covs_frpert <- dplyr::select(cell_covariates_sceptre,
                                    grna_n_nonzero_subset, grna_n_umis_subset, response_p_mito_full = response_p_mito) 
  # getting perturbation indicator
  
  stopifnot(!any(grepl(":", on_targets)))  # ensure NO targets contain ":", making it safe as delimiter
  # for each cell, we need to get the perturbations it got
  perturb_df <- cell_info |> group_by(cell_name) |> summarise(perturbation = paste0(grna_target, collapse = ":"))
  cell_covs_frpert <- left_join(
    cell_covs_frpert %>% mutate(cell_name = rownames(.)),
    perturb_df,
    by = "cell_name"
  )
  
  sce <- SingleCellExperiment(
    assays  = list(counts = response_subset_frpert),
    colData = cell_covs_frpert   # -> AnnData .obs
  )
  
  writeH5AD(
    sce,
    file = file.path(write_frperturb_fp, "response_matrix.h5ad"),
    X_name = "counts" # the name of the actual data in my `sce`
  )
  cat("   FR-perturb written.\n")
}

# response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm"))
# grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm"))
# cell_covariates = scep@covariate_data_frame
# scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds"))
# grna_target_df = scep@grna_target_data_frame

make_pos_control_gasperini(
  dataset_name,
  response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "response.odm")),
  grna_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_data, "grna.odm")),
  cell_covariates = scep@covariate_data_frame,
  scep_assn_mat = read_rds(file.path(path_to_assigns, "grna_assignment_matrix.rds")),
  grna_target_df = scep@grna_target_data_frame
)