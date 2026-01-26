rm(list=ls())
source("~/.Rprofile")

library(Matrix)
library(tidyverse)

make_sparse_row <- function(x) {
  j <- which(x != 0)
  sparseMatrix(
    i = rep.int(1L, length(j)),  # row index (all 1s)
    j = j,                       # column indices
    x = x[j],                    # nonzero values
    dims = c(1L, length(x))
  )
}

source_data <- "replogle-rd7"

input_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data", source_data, "sceptre-pipeline/"
)
output_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs", source_data, "sceptre-pipeline/"
)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(input_fp, "sceptre_object.rds"),
  response_odm_file_fp = paste0(input_fp, "response.odm"),
  grna_odm_file_fp = paste0(input_fp, "grna.odm")
)


# inputs
num_pc_cells = 200000
num_nt_cells = 10000
response_odm = ondisc::initialize_odm_from_backing_file(paste0(input_fp, "response.odm"))
grna_odm = ondisc::initialize_odm_from_backing_file(paste0(input_fp, "grna.odm"))
grna_assign_mat = readRDS(file.path(output_fp, "grna_assignment_matrix.rds"))
grna_target_df = scep@grna_target_data_frame
cell_covariates = scep@covariate_data_frame
NT_name = "non-targeting"
dataset_name = paste0(source_data, "_", "large")
min_num_cells_per_target = 20


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
  sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = unlist(xlist),
              dims = c(length(genes), length(cell_idx)))
}

# NOTE
# this has low MOI pretty hard-coded in, since we only take the 
# cells expressing exactly one gRNA
make_pos_control_replogle_rd7 <- function(num_pc_cells, num_nt_cells, min_num_cells_per_target,
                              dataset_name, response_odm,
                              grna_odm, grna_assign_mat, grna_target_df,
                              cell_covariates, NT_name) {
  
  ## 0. determining which cells express exactly one gRNA. 
  ##    This is the population we subsample from.
  
  num_grnas_expressed_per_cell <- Matrix::colSums(grna_assign_mat) # potentially big operation
  expresses_exactly_one_grna <- which(num_grnas_expressed_per_cell == 1)
  grna_assign_mat_to_use <- grna_assign_mat[,expresses_exactly_one_grna]
  
  ## 1. getting  pos control cells ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # key idea: for each response gene name, in the order 
  # they appear in the data, get the cells expressing any gRNA
  # targeting that gene. Keep doing this until we have at least
  # `num_pc_cells` many cells.
  all_response_gene_names <- rownames(response_odm) 
  target_to_grna_ids <- split(grna_target_df$grna_id, grna_target_df$grna_target)
  num_cells_per_grna <- c()
  cells_per_grna <- list()
  genes_kept <- c()
  for(gene in all_response_gene_names) {

    # a. get the corresponding guide rnas
    curr_grna_ids <- target_to_grna_ids[[gene]]
    if(is.null(curr_grna_ids)) {
      next  # skip this gene if it's not a target
    } 
    
    # b. get the cells that express each of these gRNAs
    
    # now we see how many cells for this grna_target
    # we only want targets with above our threshold number of cells
    # so we need to see how big it is first
    
    curr_cells_per_grna <- list()
    curr_num_cells_per_grna <- c()
    for(grna in curr_grna_ids) {
      curr_cells_per_grna[[grna]] <- expresses_exactly_one_grna[which(grna_assign_mat_to_use[grna, ])]
      curr_num_cells_per_grna <- c(curr_num_cells_per_grna, length(curr_cells_per_grna[[grna]]))
    }
    if(sum(curr_num_cells_per_grna) < min_num_cells_per_target) {
      next
    }
    # now we know we are keeping this cell
    # so we update our tracking
    
    genes_kept <- c(genes_kept, gene)
    cells_per_grna <- c(cells_per_grna, curr_cells_per_grna)
    num_cells_per_grna <- c(num_cells_per_grna, curr_num_cells_per_grna)

    # break once we've saved enough cells
    if(sum(num_cells_per_grna) >= num_pc_cells) {
      break
    }
  }
  
  
  ## 2. adding in the NT cells
  nt_grna_ids <- target_to_grna_ids[[NT_name]]
  num_nt_so_far <- 0
  for(grna in nt_grna_ids) {
    cells_per_grna[[grna]] <- expresses_exactly_one_grna[which(grna_assign_mat_to_use[grna, ])]
    num_cells_per_grna <- c(num_cells_per_grna, length(cells_per_grna[[grna]]))
    num_nt_so_far <- num_nt_so_far + length(cells_per_grna[[grna]])
    if(num_nt_so_far >= num_nt_cells) {
      break
    }
  }
  
  
  
  ## 3. making cell info dataframe 
  all_cell_idx <- unlist(cells_per_grna) |> setNames(NULL)
  cell_info <- data.frame(stringsAsFactors = FALSE,
                          cell_idx = all_cell_idx,
                          grna_id = rep(names(cells_per_grna), num_cells_per_grna)
  ) |>
    dplyr::left_join(grna_target_df, by = "grna_id")

  cell_covariates <- cell_covariates[all_cell_idx, ] |>
    setNames(paste0(names(cell_covariates), "_full"))
  
  ## 4. getting the final data and writing
  cat("Starting to read `response_odm`...")
  t0 <- Sys.time()
  response_mat <- odm_to_sparse_matrix(response_odm, genes_kept, all_cell_idx) |>
    `rownames<-`(genes_kept)
  t1 <- Sys.time()
  cat(" read in", round(t1-t0, 2), units(t1-t0), "\n")
  
  
  # remove any cells that didn't have any non-zero gene UMIs
  total_cell_expressions <- Matrix::colSums(response_mat)
  cells_without_expression <- total_cell_expressions == 0
  if(any(cells_without_expression)) {
    num_removed <- sum(cells_without_expression)
    # cat("  Removing", num_removed, "cells with zero expression across all genes\n")

    # Keep only cells with expression
    cells_to_keep <- !cells_without_expression

    # Update response_mat: remove columns
    response_mat <- response_mat[, cells_to_keep, drop = FALSE]

    # Update cell_info and cell_covariates: remove rows
    cell_info <- cell_info[cells_to_keep, ]
    cell_covariates <- cell_covariates[cells_to_keep, ]

    # Update all_cell_idx: remove elements
    removed_cell_indices <- all_cell_idx[cells_without_expression]
    all_cell_idx <- all_cell_idx[cells_to_keep]

    # Update cells_per_grna: remove filtered cell indices from each gRNA's list
    cells_per_grna <- lapply(cells_per_grna, function(cells) {
      setdiff(cells, removed_cell_indices)
    })

    # Recalculate num_cells_per_grna after filtering
    num_cells_per_grna <- sapply(cells_per_grna, length)

    # Remove gRNAs with no cells after filtering
    empty_grnas <- num_cells_per_grna == 0
    if(any(empty_grnas)) {
      num_empty <- sum(empty_grnas)
      cat("  Removing", num_empty, "gRNAs with no cells after expression filtering\n")
      cells_per_grna <- cells_per_grna[!empty_grnas]
      num_cells_per_grna <- num_cells_per_grna[!empty_grnas]
    }
  }

  cat("Getting ready to write", dataset_name, "\b...\n")
  
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/pos-control/input_data", dataset_name
  )
  dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)
  write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)
  cat("   `cell_info.csv` written.\n")
  
  write.csv(cell_covariates, file.path(write_fp, "cell_covariates.csv"), row.names = FALSE)
  cat("   `cell_covariates.csv` written.\n")
  
  ## 4a. sceptre -----------------------------------------------------
  write_sceptre_fp <- file.path(write_fp, "sceptre")
  dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)

  grna_target_df_kept <- grna_target_df |>
    dplyr::filter(grna_id %in% names(cells_per_grna))

  # Create binary assignment indicator matrix (not UMI counts)
  # grna_indicator_matrix[grna, cell] = 1 if cell is assigned to grna, 0 otherwise
  locs <- lapply(seq_along(cells_per_grna), function(i) {
    grna <- names(cells_per_grna)[i]
    assigned_cells <- cells_per_grna[[grna]]
    list(
      i = rep(i, length(assigned_cells)),
      j = match(assigned_cells, all_cell_idx),
      x = rep(1, length(assigned_cells))
    )
  })
  grna_indicator_matrix <- sparseMatrix(
    i = lapply(locs, `[[`, "i") |> unlist(),
    j = lapply(locs, `[[`, "j") |> unlist(),
    x = lapply(locs, `[[`, "x") |> unlist(),
    dims = c(length(cells_per_grna), length(all_cell_idx))
  )|>
    `rownames<-`(names(cells_per_grna))

  # Compute sceptre-specific covariates from the subsetted gRNA UMI data
  # Extract actual UMI counts from grna_odm for this subset of gRNAs and cells
  grna_matrix_subset <- odm_to_sparse_matrix(grna_odm, names(cells_per_grna), all_cell_idx)
  

  # Compute covariates from the UMI subset (not from binary grna_indicator_matrix)
  grna_n_nonzero_subset <- Matrix::colSums(grna_matrix_subset != 0)
  grna_n_umis_subset <- Matrix::colSums(grna_matrix_subset)

  # Create sceptre-specific cell covariates with subset-based gRNA metrics
  cell_covariates_sceptre <- cbind(
    cell_covariates,
    data.frame(
      grna_n_nonzero_subset = grna_n_nonzero_subset,
      grna_n_umis_subset = grna_n_umis_subset
    )
  )

  saveRDS(response_mat, file.path(write_sceptre_fp, "response_matrix.rds"))
  saveRDS(grna_indicator_matrix, file.path(write_sceptre_fp, "grna_matrix.rds"))
  write.csv(cell_covariates_sceptre, file.path(write_sceptre_fp, "cell_covariates.csv"), row.names = FALSE)
  write.csv(grna_target_df_kept, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names = FALSE)
  cat("   sceptre written.\n")
  
  ## 4b. Mixscale -----------------------------------------------------
  write_mixscale_fp <- file.path(write_fp, "mixscale")
  dir.create(write_mixscale_fp, showWarnings = FALSE, recursive = TRUE)
  cell_names <- paste0("CELL_", cell_info$cell_idx)
  
  saveRDS(response_mat |> `colnames<-`(cell_names),
          file.path(write_mixscale_fp, "response_matrix.rds"))
  saveRDS(setNames(cell_info$grna_target, cell_names),
          file.path(write_mixscale_fp, "assignments.rds"))
  cat("   mixscale written.\n")
  
  ## 4c. FR-Perturb -----------------------------------------------------
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
  # this is low MOI so the perturbation column is simple
  response_mat_frpert <- response_mat |> `colnames<-`(cell_names)
  
  # these get added to .obs of the anndata object
  cell_covs_frpert <- dplyr::select(cell_covariates_sceptre,
                                    grna_n_nonzero_subset, grna_n_umis_subset, response_p_mito_full) |>
    mutate(perturbation = cell_info$grna_target)
  
  sce <- SingleCellExperiment(
    assays  = list(counts = response_mat_frpert),
    colData = cell_covs_frpert   # -> AnnData .obs
  )
  
  writeH5AD(
    sce,
    file = file.path(write_frperturb_fp, "response_matrix.h5ad"),
    X_name = "counts" # the name of the actual data in my `sce`
  )
  cat("   FR-perturb written.\n")

  cat(dataset_name, "finished.\n\n")
}


make_pos_control_replogle_rd7(
  num_pc_cells = 10000, num_nt_cells = 2000, min_num_cells_per_target = 100,
  dataset_name = paste0(source_data, "_", "small"),
  response_odm = response_odm, grna_odm = grna_odm,
  grna_assign_mat = grna_assign_mat, grna_target_df = grna_target_df,
  cell_covariates = cell_covariates, NT_name = NT_name
)

make_pos_control_replogle_rd7(
  num_pc_cells = 50000, num_nt_cells = 5000, min_num_cells_per_target = 300,
  dataset_name = paste0(source_data, "_", "medium"),
  response_odm = response_odm, grna_odm = grna_odm,
  grna_assign_mat = grna_assign_mat, grna_target_df = grna_target_df,
  cell_covariates = cell_covariates, NT_name = NT_name
)

# make_pos_control_replogle_rd7(
#   num_pc_cells = 50000, num_nt_cells = 5000, min_num_cells_per_target = 300,
#   dataset_name = paste0(source_data, "_", "medium_v2"),
#   response_odm = response_odm, grna_odm = grna_odm,
#   grna_assign_mat = grna_assign_mat, grna_target_df = grna_target_df,
#   cell_covariates = cell_covariates, NT_name = NT_name
# )



# make_pos_control_replogle_rd7(
#   num_pc_cells = 200000, num_nt_cells = 10000, min_num_cells_per_target = 20,
#   dataset_name = paste0(source_data, "_", "large"),
#   response_odm = response_odm, grna_odm = grna_odm,
#   grna_assign_mat = grna_assign_mat, grna_target_df = grna_target_df,
#   cell_covariates = cell_covariates, NT_name = NT_name
# )





## some testing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fp1 <- "~/data/projects/sceptre3/benchmarking/association/pos-control/input_data/replogle-rd7_medium/"
# fp2 <- "~/data/projects/sceptre3/benchmarking/association/pos-control/input_data/replogle-rd7_medium_v2/"
# fps <- c(fp1, fp2)
# 
# cc <- lapply(fps, function(fp) read.csv(file.path(fp, "cell_info.csv")))
# stopifnot(identical(cc[[1]], cc[[2]]))
# 
# cc <- lapply(fps, function(fp) read.csv(file.path(fp, "cell_covariates.csv")))
# stopifnot(identical(cc[[1]], cc[[2]]))
# 
# cc <- lapply(fps, function(fp) read.csv(file.path(fp, "cell_covariates.csv")))
# stopifnot(identical(cc[[1]], cc[[2]]))
# 
# cc <- lapply(fps, function(fp) readRDS(file.path(fp, "/sceptre/grna_matrix.rds")))
# stopifnot(identical(cc[[1]], cc[[2]]))
# 
# cc <- lapply(fps, function(fp) readRDS(file.path(fp, "/sceptre/response_matrix.rds")))
# stopifnot(identical(cc[[1]], cc[[2]]))
# 
# cc <- lapply(fps, function(fp) read.csv(file.path(fp, "/sceptre/grna_target_data_frame.csv")))
# stopifnot(identical(cc[[1]], cc[[2]]))
