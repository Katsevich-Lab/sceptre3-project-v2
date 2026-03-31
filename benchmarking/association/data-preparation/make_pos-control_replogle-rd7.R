rm(list=ls())
source("~/.Rprofile")

library(tidyverse)

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "utils_data_prep.R"))
source(file.path(script_dir, "neg_control_functions.R"))
source(file.path(script_dir, "computational_benchmarking_functions.R"))

make_pos_control_replogle_rd7 <- function(
    dataset_name,
    response_odm,
    grna_odm,
    cell_covariates,
    scep_assn_mat,
    grna_target_df,
    genes_passing_qc,
    num_targets,
    max_num_cells,
    only_consider_targets_with_this_many_cells,
    nt_name = "non-targeting",
    random_seed = 34534) {

  set.seed(random_seed)
  cat("Random seed set to", random_seed, "\n")
  
  # 1. get all on-targets that pass QC and have enough cells, and sample from there
  all_on_targets <- intersect(
    grna_target_df$grna_target,
    genes_passing_qc
  )
  all_on_target_guides = grna_target_df$grna_id[grna_target_df$grna_target %in% all_on_targets]
  
  # because we will likely be downsampling to get enough on-targets, 
  # we want to only take ones with enough cells to survive this downsampling
  on_targets_with_enough_cells <- data.frame(
    grna_id = all_on_target_guides,
    cells_per_guide = Matrix::rowSums(scep_assn_mat[all_on_target_guides,]),
    stringsAsFactors = FALSE
  ) |>
    left_join(dplyr::select(grna_target_df, -grna_group), by = "grna_id") |>
    group_by(grna_target) |>
    summarize(cells_per_target = sum(cells_per_guide)) |>
    dplyr::filter(cells_per_target >= only_consider_targets_with_this_many_cells)
  if(nrow(on_targets_with_enough_cells) < num_targets) {
    stop("Not enough on-targets have enough cells per target.")
  }
  
  targets <- sample(on_targets_with_enough_cells$grna_target, num_targets)
  genes = targets # the same for on-target
  cat("There are", length(all_on_targets), "on-targets.\n")
  cat("   ", nrow(on_targets_with_enough_cells), "have at least", only_consider_targets_with_this_many_cells, "cells.\n")
  cat("   ", length(targets), "on-targets sampled.\n")
  
  # 2. get all guides for these targets, and all cells expressing these guides
  guides = grna_target_df$grna_id[grna_target_df$grna_target %in% c(targets, nt_name)]
  candidate_cells = which(Matrix::colSums(scep_assn_mat[guides,]) > 0)
  
  cat("Initially", length(candidate_cells), "cells found.\n")
  if(length(candidate_cells) >= max_num_cells) {
    cat("Downsampling to", max_num_cells, "cells.\n")
    candidate_cells <- sample(candidate_cells, max_num_cells)
  } else {
    stop("Not enough cells for these parameters.")
  }
  
  # 3. get the response matrix subset, and remove any cells with all-0 response
  response_subset = odm_to_sparse_matrix(response_odm, genes, candidate_cells)
  # remove any cells that have all 0 expression
  non_zero_cells = Matrix::colSums(response_subset) > 0
  response_subset = response_subset[, non_zero_cells]
  cell_idx = candidate_cells[non_zero_cells]
  cat("Response matrix created.", length(cell_idx), "cells remain after removing all-0 expression cells.\n")
  
  # 4. make dummy grna matrix and enforce low MOI
  grna_indicator = (scep_assn_mat[guides, cell_idx] + 0) |> # +0 converts to integer from logical
    enforce_single_guide_per_cell(random_seed = random_seed)
  cat("grna indicator matrix made and low moi enforced.\n")
  
  # 5. prepare metadata and the final necessary pieces 
  metadata = prepare_cell_metadata_low_moi(grna_indicator, grna_target_df, cell_covariates, cell_idx)
  
  cell_info <- metadata$cell_info
  cell_covariates_subset = metadata$cell_covariates_subset |>
    dplyr::select(-batch)
  
  grna_target_df_subset = grna_target_df[grna_target_df$grna_id %in% rownames(grna_indicator), ]
  
  stopifnot(nrow(cell_covariates_subset) == ncol(response_subset))
  
  # what's the least sample size for an on-target?
  num_cells_per_target <- cell_info |>
    group_by(grna_target) |>
    summarize(n=n()) |>
    pull(n) |>
    sort()
  cat("Smallest num. cells per target:", paste0(num_cells_per_target[1:5], collapse=", "), "\n")
  plt <- data.frame(num_cells_per_target = num_cells_per_target) |>
    ggplot(aes(x = num_cells_per_target)) +
    geom_histogram(bins=60) +
    scale_y_log10() +
    theme_bw() +
    labs(x = "num. cells for each on-target",
         title = "Sample size for on-targets")
  print(plt)
  
  
  # 6: Write outputs
  write_fp <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "association/pos-control/input_data", dataset_name
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

  # Write SCEPTRE format
  write_sceptre_output(
    response_matrix = response_subset,
    grna_matrix = grna_indicator,
    cell_covariates = cell_covariates_subset,
    grna_target_df = grna_target_df_subset,
    discovery_pairs = NULL, # not needed for pos control
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

  cat("\n Positive control dataset creation complete!\n")
  cat("Output directory:", write_fp, "\n")
}



source_data <- "replogle-rd7"
num_targets = 1500
max_num_cells = 100000
gene_qc_thresh = 7
# almost 500k cells initially are selected, so we will be downsampling to ~20%
only_consider_targets_with_this_many_cells = 150 
random_seed = 876567

dataset_name = paste0(source_data, "_pos_ctrl_ntargets=", num_targets, "_ncells=", max_num_cells / 1000, "k_gene_thresh=",gene_qc_thresh)

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

response_odm = ondisc::initialize_odm_from_backing_file(paste0(input_fp, "response.odm"))
grna_odm = ondisc::initialize_odm_from_backing_file(paste0(input_fp, "grna.odm"))
scep_assn_mat = readRDS(file.path(output_fp, "grna_assignment_matrix.rds"))
grna_target_df = scep@grna_target_data_frame
cell_covariates = scep@covariate_data_frame

moi_observed <- Matrix::colSums(scep_assn_mat) |> mean()

total_num_guides <- grna_target_df |>
  filter(! grna_target %in% c("nt_off_target", "unknown")) |>
  pull(grna_id) |>
  unique() |>
  length()

avg_num_guides_per_target <- grna_target_df |>
  filter(! grna_target %in% c("nt_off_target", "unknown")) |>
  group_by(grna_target) |>
  summarize(n=n()) |>
  pull(n) |>
  mean()

gene_summary_stats = read_csv(file.path(input_fp, "gene_summary_stats.csv"), show_col_types = FALSE)

genes_passing_qc <- gene_summary_stats |>
  filter(gene_n_nonzero * moi_observed / total_num_guides * avg_num_guides_per_target >= gene_qc_thresh) |>
  pull(gene)

cat("There are", length(genes_passing_qc), "genes passing QC.\n")

make_pos_control_replogle_rd7(
  dataset_name = dataset_name,
  response_odm = response_odm,
  grna_odm = grna_odm,
  cell_covariates = cell_covariates,
  scep_assn_mat = scep_assn_mat,
  grna_target_df = grna_target_df,
  genes_passing_qc = genes_passing_qc,
  num_targets = num_targets,
  max_num_cells = max_num_cells,
  only_consider_targets_with_this_many_cells = only_consider_targets_with_this_many_cells,
  nt_name = "non-targeting",
  random_seed = random_seed
)
