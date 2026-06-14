
# rm(list=ls())
# source("~/.Rprofile")
# library(tidyverse)
# library(SingleCellExperiment)
# library(zellkonverter)
# 
# # for object writing utils
# script_dir <- dirname(sys.frame(1)$ofile)
# source(file.path(script_dir, "../../association/data-preparation/utils_data_prep.R"))
# 

grna_simulator_iid_sum_process <- function(num_guides, lib_sizes_scaled, params, seed, return_as_matrices = TRUE) {
  
  set.seed(seed)
  pert_sparse_pieces <- y_sparse_pieces <- list()
  num_cells = length(lib_sizes_scaled)
  
  for (i in seq_len(num_guides)) {
    
    # 1. generate non-perturbed counts
    is_endo <- which(rbinom(num_cells, 1, params$prob_endo) == 1)
    y <- rpois(num_cells, lib_sizes_scaled * params$mu_exo) 
    
    y[is_endo] <- y[is_endo] +
      rnbinom(
        n = length(is_endo),
        size = params$theta_endo,
        mu = lib_sizes_scaled[is_endo] * params$mu_endo
      )
    
    # 2. select perturbed cells and add true counts
    is_pert <- which(rbinom(num_cells, 1, params$prob_pert) == 1)
    
    y[is_pert] <- y[is_pert] + 
      rnbinom(
        n = length(is_pert),
        size = params$theta_pert,
        mu = lib_sizes_scaled[is_pert] * params$mu_pert
      )
    
    # 3. collect what we need for sparse matrices
    y_pos <- which(y > 0) 
    
    pert_sparse_pieces[[i]] <- list(
      row_id = rep.int(i, length(is_pert)),
      col_id = is_pert
    )
    
    y_sparse_pieces[[i]] <- list(
      row_id = rep.int(i, length(y_pos)),
      col_id = y_pos,
      vals = y[y_pos]
    )
  }
  
  out <- list(
    lib_sizes_scaled = lib_sizes_scaled,
    num_guides = num_guides, 
    params = params,
    seed = seed
  )
  
  if (!return_as_matrices) {
    
    out$pert_sparse_pieces <- pert_sparse_pieces
    out$y_sparse_pieces <- y_sparse_pieces
    
  } else {
    
    pert_i <- unlist(lapply(pert_sparse_pieces, `[[`, "row_id"), use.names = FALSE)
    pert_j <- unlist(lapply(pert_sparse_pieces, `[[`, "col_id"), use.names = FALSE)
    
    y_i <- unlist(lapply(y_sparse_pieces, `[[`, "row_id"), use.names = FALSE)
    y_j <- unlist(lapply(y_sparse_pieces, `[[`, "col_id"), use.names = FALSE)
    y_x <- unlist(lapply(y_sparse_pieces, `[[`, "vals"), use.names = FALSE)
    
    true_perts <- Matrix::sparseMatrix(
      i = pert_i,
      j = pert_j,
      x = rep.int(1L, length(pert_i)),
      dims = c(num_guides, num_cells)
    )
    
    grna_matrix_sim <- Matrix::sparseMatrix(
      i = y_i,
      j = y_j,
      x = y_x,
      dims = c(num_guides, num_cells)
    )
    
    out$true_perts <- true_perts
    out$grna_matrix_sim <- grna_matrix_sim
  }
  
  return(out)
}

make_grna_simulation_sum_process <- function(
    dataset_name, params_list, num_cells, num_guides_per_param, scep, seed, output_base_dir = NULL
) {
  
  set.seed(seed)
  if(is.null(names(params_list))) {
    stop("`params_list` must have names. They are used for guide names.")
  }
  if(is.null(output_base_dir)) {
    output_base_dir =  file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")
  }
  real_cell_locs = which(scep@covariate_data_frame$grna_n_nonzero > 0) |>
    sample(num_cells, replace = FALSE) |>
    sort()
  cell_names <- paste0("cell_idx_", real_cell_locs)
  
  lib_sizes_scaled = scep@covariate_data_frame$grna_n_nonzero[real_cell_locs] %>%
    `/`(mean(.))
  
  sim_results <- vector("list", length(params_list))
  
  for(param_idx in seq_along(params_list)) {
    sims_curr <- grna_simulator_iid_sum_process(
      num_guides = num_guides_per_param, 
      lib_sizes_scaled = lib_sizes_scaled,
      params = params_list[[param_idx]],
      seed = seed + param_idx
    )
    
    guides_names <- paste0("sim_", names(params_list)[[param_idx]], "_", 1:num_guides_per_param)

    sim_results[[param_idx]] <- list(
      grna_matrix = sims_curr$grna_matrix_sim |>
        `dimnames<-`(list(guides_names, cell_names)),
      true_perts = sims_curr$true_perts |>
        `dimnames<-`(list(guides_names, cell_names))
    )
  }

  grna_sims <- lapply(sim_results, function(l) l$grna_matrix) |>
    do.call(what = rbind) |>
    as("RsparseMatrix")
  true_pert_matrix <- lapply(sim_results, function(l) l$true_perts) |>
    do.call(what = rbind) |>
    as("RsparseMatrix")

  output_dir <- file.path(output_base_dir, dataset_name)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  setup_anndata_environment()
  
  
  # saving objects ~~~~~~~~~~~~~~~~
  
  cov_dat <- data.frame(
    real_cell_locs = real_cell_locs,
    lib_sizes_scaled =  lib_sizes_scaled
  )
  sim_params <- list(
    params_list = params_list, num_guides_per_param = num_guides_per_param,
    seed = seed, num_cells = num_cells
  )
  
  saveRDS(true_pert_matrix, file.path(output_dir, "true_pert_matrix.rds"))
  saveRDS(cov_dat, file.path(output_dir, "cell_scaling_and_locs.rds"))
  saveRDS(sim_params, file.path(output_dir, "sim_params.rds"))
  
  # 4a,b. crispat and pertpy ~~~~~~~~~~~~~~~~~~~~~~
  
  sce <- SingleCellExperiment(
    assays  = list(counts = grna_sims)
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
  Matrix::writeMM(grna_sims, file.path(cleanser_dir, "grna_matrix.mtx"))
  
  # 4d. sceptre  ~~~~~~~~~~~~~~~~~~~~~~
  
  sceptre_dir <- file.path(output_dir, "sceptre")
  dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
  print(paste("  Writing sceptre gRNA matrix to", sceptre_dir))
  
  covariates_subset <- scep@covariate_data_frame[real_cell_locs,] |>
    dplyr::select(
      response_n_umis_full = response_n_umis,
      response_n_nonzero_full = response_n_nonzero,
      grna_n_umis_full = grna_n_umis,
      grna_n_nonzero_full = grna_n_nonzero
    ) |>
    as.data.frame()
  
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
  
  grna_target_df <- data.frame(
    grna_id = rownames(grna_sims),
    grna_target = paste0("grna_target_", 1:nrow(grna_sims))
  )
  grna_target_df$grna_target[1] = "non-targeting"
  
  write_sceptre_output(
    response_matrix = dummy_response_mat,
    grna_matrix = grna_sims,
    cell_covariates = covariates_subset,
    grna_target_df,
    discovery_pairs = NULL,
    formula_object = assign_fmla,
    output_path = sceptre_dir
  )
  
  cat("scp -r", dataset_name, "hpc3:/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
}




# path_to_data <- file.path(
#   .get_config_path("LOCAL_BENCHMARKING_DIR"),
#   "guide_assignment/input_data",
#   "replogle-rd7",
#   "sceptre-pipeline"
# )
# scep_repl <- sceptre::read_ondisc_backed_sceptre_object(
#   sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
#   response_odm_file_fp = file.path(path_to_data, "response.odm"),
#   grna_odm_file_fp = file.path(path_to_data, "grna.odm")
# )
# 
# 
# params_np_highvar <- list(
#   prob_endo = 0.002,
#   mu_endo = 10,
#   theta_endo = 1,
#   mu_exo = 0.01
# )
# 
# params_np_lowvar <- list(
#   prob_endo = 0.001,
#   mu_endo = 5,
#   theta_endo = 10,
#   mu_exo = 0.01
# )
# 
# params_p_3 <- list(
#   prob_pert = 0.01,
#   mu_pert = 1000,
#   theta_pert = 5
# )
# 
# params_p_2 <- list(
#   prob_pert = 0.01,
#   mu_pert = 750,
#   theta_pert = 2.5
# )
# 
# params_p_1 <- list(
#   prob_pert = 0.01,
#   mu_pert = 500,
#   theta_pert = 1
# )
# 
# params_list <- list(
#   NPlowvar_P1 = c(params_np_lowvar, params_p_1),
#   NPlowvar_P2 = c(params_np_lowvar, params_p_2),
#   NPlowvar_P3 = c(params_np_lowvar, params_p_3),
#   NPhighvar_P1 = c(params_np_highvar, params_p_1),
#   NPhighvar_P2 = c(params_np_highvar, params_p_2),
#   NPhighvar_P3 = c(params_np_highvar, params_p_3)
# )
# 
# dataset_name = "new-sims-test"
# make_grna_simulation_sum_process(
#   dataset_name = dataset_name,
#   params_list = params_list,
#   num_cells = 50000,
#   num_guides_per_param = 100,
#   scep = scep_repl,
#   seed = 654
# )


plot_umi_histogram_real_vs_sim <- function(
    umis_real,
    umis_sim,
    is_pert,
    bins = 80,
    title = NULL,
    x_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                 2000, 5000, 10000, 20000, 50000, 100000)
    # x_breaks = c(0,1,2,3,4,5,6,7,8,9,10 * 2^(0:14))
) {
  stopifnot(all(umis_real >= 0, na.rm = TRUE))
  stopifnot(all(umis_sim >= 0, na.rm = TRUE))
  stopifnot(length(umis_sim) == length(is_pert))
  
  umis_real <- umis_real[is.finite(umis_real)]
  
  keep_sim <- is.finite(umis_sim) & !is.na(is_pert)
  umis_sim <- umis_sim[keep_sim]
  is_pert <- is_pert[keep_sim]
  
  x_real <- log1p(umis_real)
  x_sim <- log1p(umis_sim)
  x_all <- c(x_real, x_sim)
  
  # Shared x-axis bins across real and simulated data.
  h_all <- hist(x_all, breaks = bins, plot = FALSE)
  breaks <- h_all$breaks
  
  make_bin_df <- function(x, group, panel) {
    bin_id <- cut(
      x,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    )
    
    df <- data.frame(
      bin_id = bin_id,
      group = group
    )
    
    df <- df[!is.na(df$bin_id), , drop = FALSE]
    
    tab <- as.data.frame(table(df$bin_id, df$group))
    names(tab) <- c("bin_id", "group", "count")
    
    tab$bin_id <- as.integer(as.character(tab$bin_id))
    tab$count <- as.integer(tab$count)
    tab <- tab[tab$count > 0, , drop = FALSE]
    
    tab$xmin <- breaks[tab$bin_id]
    tab$xmax <- breaks[tab$bin_id + 1L]
    tab$panel <- panel
    
    tab
  }
  
  df_real <- make_bin_df(
    x = x_real,
    group = "Real",
    panel = "Real"
  )
  
  df_sim <- make_bin_df(
    x = x_sim,
    group = ifelse(is_pert, "Perturbed", "Non-perturbed"),
    panel = "Simulated"
  )
  
  # Stack simulated bars manually so total bar height is the true bin count.
  # Real panel has only one group.
  df <- rbind(df_real, df_sim)
  
  df$group <- factor(
    df$group,
    levels = c("Real", "Non-perturbed", "Perturbed")
  )
  
  df <- df[order(df$panel, df$bin_id, df$group), , drop = FALSE]
  
  df$ymin <- 0
  df$ymax <- 0
  
  split_ids <- split(seq_len(nrow(df)), paste(df$panel, df$bin_id, sep = "___"))
  
  for (idx in split_ids) {
    counts <- df$count[idx]
    cum_counts <- cumsum(counts)
    
    df$ymin[idx] <- c(0, head(cum_counts, -1))
    df$ymax[idx] <- cum_counts
  }
  
  # Log y-scale cannot display 0, so bars start visually at 1.
  # The top of each stacked bar is still the true bin count.
  df$ymin_plot <- pmax(df$ymin, 1)
  df$ymax_plot <- df$ymax
  
  df <- df[df$ymax_plot > df$ymin_plot, , drop = FALSE]
  
  # Keep only breaks inside the observed range.
  max_count <- max(c(umis_real, umis_sim), na.rm = TRUE)
  x_breaks <- x_breaks[x_breaks <= max_count]
  if (!0 %in% x_breaks) {
    x_breaks <- c(0, x_breaks)
  }
  
  facet_scales <- if (length(umis_real) == length(umis_sim)) {
    "fixed"
  } else {
    "free_y"
  }
  
  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin_plot,
        ymax = ymax_plot,
        fill = group
      )
    ) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = facet_scales) +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      limits = range(breaks),
      name = "UMI count"
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      name = "Number of cells"
    ) +
    ggplot2::labs(
      title = title,
      fill = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 7,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )
}


# # load the grna matrix and odm, make side-by-side 
# path_to_new_data = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data", dataset_name)
# grna_mat = file.path(path_to_new_data,  "sceptre/grna_matrix.rds") |>
#   read_rds()
# true_perts = file.path(path_to_new_data,  "true_pert_matrix.rds") |>
#   read_rds()
# 
# grna_odm = ondisc::initialize_odm_from_backing_file(
#   file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
#             "guide_assignment/input_data/replogle-rd7/sceptre-pipeline/grna.odm")
# )
# 
# i=1
# j=1
# plot_umi_histogram_real_vs_sim(
#   umis_real = grna_odm[j,],
#   umis_sim = grna_mat[i,],
#   is_pert = true_perts[i,],
#   title = paste0("real grna ", j, "; sim grna ", i)
# )


