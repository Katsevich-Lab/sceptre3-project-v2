
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)


rm(list=ls())
source("~/.Rprofile")

# Source utility functions
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "../../association/data-preparation/utils_data_prep.R"))

output_base_dir = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")

make_grna_subset <- function(source_data, dataset_name, output_base_dir, num_guides, num_cells, batch_col = NULL) {
  
  # 0. load objects ~~~~~~~~~~~~~~~~~~~~~~~
  path_to_data <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "guide_assignment/input_data", source_data,
    "sceptre-pipeline"
  )
  cat("Using", source_data, "\n")
  grna_odm <- ondisc::initialize_odm_from_backing_file(
    file.path(path_to_data, "grna.odm")
  )
  cat("grna.odm :", nrow(grna_odm), "guides ;", ncol(grna_odm), "cells.\n")
  
  scep <- sceptre::read_ondisc_backed_sceptre_object(
    sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
    response_odm_file_fp = file.path(path_to_data, "response.odm"),
    grna_odm_file_fp = file.path(path_to_data, "grna.odm")
  )
  grna_summary_stats <- read_csv(file.path(path_to_data, "grna_summary_stats.csv"), show_col_types = FALSE) |>
    mutate(grna=rownames(grna_odm))
  
  # 1. get guides in a linear seq of decerasing sparsity
  tot_num_guides = nrow(grna_summary_stats)
  guide_idx = seq(50, tot_num_guides-50, length = num_guides) |> round()
  stopifnot(!any(duplicated(guide_idx)))
  
  guides = (grna_summary_stats |>
              arrange(grna_nonzero) |>
              pull(grna))[guide_idx]
  
  
  # 2. get cells and remove all-0 columns
  candidate_cells <- sample(ncol(grna_odm), num_cells)
  
  grna_subset <- odm_to_sparse_matrix(grna_odm, guides, candidate_cells)
  cells_with_umis <- Matrix::colSums(grna_subset) > 0
  
  cell_idx <- candidate_cells[cells_with_umis]
  grna_subset <- grna_subset[,cells_with_umis]
  
  cat(" ~~ Keeping", nrow(grna_subset), "guides and", ncol(grna_subset), "cells ~~ \n")
  # 3. prepare objects
  
  if(is.null(colnames(grna_odm))) {
    cell_names <- paste0("cell_idx_", cell_idx)
  } else {
    cell_names <- colnames(grna_odm)[cell_idx]
  }
  colnames(grna_subset) <- cell_names
  
  # 4. write
  output_dir <- file.path(output_base_dir, dataset_name)
  
  setup_anndata_environment()
  
  # 4a,b. crispat and pertpy ~~~~~~~~~~~~~~~~~~~~~~
  
  sce <- SingleCellExperiment(
    assays  = list(counts = grna_subset)
  )
  
  crispat_dir <- file.path(output_dir, "crispat")
  dir.create(crispat_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing crispat to", crispat_dir))
  
  writeH5AD(
    sce,
    file = file.path(crispat_dir, "grna_matrix.h5ad")
  )
  
  pertpy_dir <- file.path(output_dir, "pertpy")
  dir.create(pertpy_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing pertpy to", pertpy_dir))
  
  writeH5AD(
    sce,
    file = file.path(pertpy_dir, "grna_matrix.h5ad")
  )
  
  # 4c: cleanser  ~~~~~~~~~~~~~~~~~~~~~~
  
  cleanser_dir <- file.path(output_dir, "cleanser")
  dir.create(cleanser_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing cleanser to", cleanser_dir))
  Matrix::writeMM(grna_subset, file.path(cleanser_dir, "grna_matrix.mtx"))
  
  # 4d. sceptre  ~~~~~~~~~~~~~~~~~~~~~~
  
  sceptre_dir <- file.path(output_dir, "sceptre")
  dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing sceptre gRNA matrix to", sceptre_dir))
  
  covariates_subset <- scep@covariate_data_frame[cell_idx,] |>
    dplyr::select(response_n_umis_full = response_n_umis, response_n_nonzero_full = response_n_nonzero) |>
    as.data.frame()
  if(!is.null(batch_col)) {
    covariates_subset[,batch_col] = as.data.frame(scep@covariate_data_frame)[cell_idx, batch_col] # hacky just to get this done
  }
  
  dummy_response_mat <- Matrix::sparseMatrix(
    i = integer(0),  # No non-zero entries
    j = integer(0),
    x = numeric(0),
    dims = c(1L, length(cell_names)),
    dimnames = list(
      "FAKE_GENE_1",  # Single fake gene
      cell_names  # Use same cell names as grna_matrix
    )
  )
  
  assign_fmla <- "~ log(response_n_umis_full + 1) + log(response_n_nonzero_full + 1) + log(grna_n_umis + 1) + log(grna_n_nonzero + 1)"
  if(!is.null(batch_col)) {
    assign_fmla <- paste0(assign_fmla, " + ", batch_col)
  }
  
  grna_target_df <- scep@grna_target_data_frame |>
    dplyr::filter(grna_id %in% rownames(grna_subset)) |>
    as.data.frame()
  if(! "non-targeting" %in% grna_target_df$grna_target) {
    grna_target_df[1,"grna_target"] <- "non-targeting"
  }
  
  write_sceptre_output(
    response_matrix = dummy_response_mat,
    grna_matrix = grna_subset,
    cell_covariates = covariates_subset,
    grna_target_df,
    discovery_pairs = NULL,
    formula_object = assign_fmla,
    output_path = sceptre_dir
  )
  
  
}


dataset_params = data.frame(
  source_data = rep(c( "replogle-rd7", "gasperini"), each = 3),
  num_cells=  100000,
  num_guides = rep(c(100,200,400), times = 2)
) |>
  mutate(
    dataset_name = paste0(source_data, "_assign_nguides=", num_guides, "_ncells=", num_cells / 1000, "k")
  )

for(i in 1:nrow(dataset_params)) {
  source_data = dataset_params$source_data[i]
  if(source_data == "gasperini") {
    batch_col = "prep_batch"
  } else {
    batch_col = NULL
  }
  make_grna_subset(
    source_data = source_data,
    dataset_name = dataset_params$dataset_name[i],
    output_base_dir = output_base_dir,
                   num_guides = dataset_params$num_guides[i],
    num_cells = dataset_params$num_cells[i],
    batch_col = batch_col
    )
  cat("Done with", i," \n")
}


for(i in 1:nrow(dataset_params)) {
  cat("scp -r", dataset_params$dataset_name[i], " hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/input_data \n")
}
cat("\n")
