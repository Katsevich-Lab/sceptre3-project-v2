library(MASS)
library(Matrix)
source("~/.Rprofile")
source(file.path(
  .get_config_path("LOCAL_SCEPTRE3_REPO_DIR"),
  "benchmarking/guide-assignment/data-preprocessing/convert_odm_to_h5ad.R"
))

simulate_grna_matrix <- function(
    grna_matrix_rows,                    # Matrix: rows = targets/gRNAs, cols = cells
    formula,                             # One-sided formula for the covariates upon which gRNA UMI counts are regressed
    cell_covariates,                     # Data frame with covariates (rows = cells)
    num_grnas_per_perturbation = 1,          # int, the amount of gRNAs simulated per each real gRNA & perturbation effect size pair
    perturbation_effect_sizes = data.frame(non_pert_effect = -4, pert_effect = 6, theta = 4.5),
    prob_perturbed = 0.02                 # Probability a cell receives the perturbation
) {
  
  num_cells <- nrow(cell_covariates)
  num_targets <- nrow(grna_matrix_rows)
  
  # Validate inputs
  stopifnot(
    is.matrix(grna_matrix_rows),
    is.data.frame(cell_covariates),
    ncol(grna_matrix_rows) == num_cells,
    is.numeric(num_grnas_per_perturbation),
    num_grnas_per_perturbation > 0,
    is.numeric(prob_perturbed),
    length(prob_perturbed) >= 1,
    all(prob_perturbed >= 0 & prob_perturbed <= 1),
    inherits(formula, "formula"),
    length(formula) == 2,  # One-sided formula
    is.data.frame(perturbation_effect_sizes),
    all(c("non_pert_effect", "pert_effect", "theta") %in% names(perturbation_effect_sizes)),
    nrow(perturbation_effect_sizes) >= 1
  )
  
  # Construct full formula internally (response name doesn't matter to user)
  full_formula <- as.formula(paste("..y.. ~", as.character(formula)[2]))
  
  # Storage for simulated gRNA rows, perturbation indicators, and model info
  simulated_rows <- list()
  perturbation_rows <- list()
  row_names <- character()
  model_info <- list()
  
  # For each perturbation probability (outermost loop)
  for (prob_idx in seq_along(prob_perturbed)) {
    prob_pert <- prob_perturbed[prob_idx]
    
    # For each target (input gRNA row)
    for (target_idx in seq_len(num_targets)) {
      
      # Get count vector for this target
      count_vector <- grna_matrix_rows[target_idx, ]
      
      # Fit negative binomial model
      fit_data <- cbind(..y.. = count_vector, cell_covariates)
      mod <- glm.nb(full_formula, data = fit_data)
      
      # Get target name (or use index if no rownames)
      sim_grna_name <- if (!is.null(rownames(grna_matrix_rows))) {
        paste0("sim_grna=", rownames(grna_matrix_rows)[target_idx])
      } else {
        paste0("sim_grna_", target_idx)
      }
      
      # sim_grna_name <- paste0("sim_grna_", target_idx)
      
      # Store model information
      model_info[[sim_grna_name]] <- list(
        coefficients = summary(mod)$coef,
        theta_from_model = mod$theta
      )
      
      # Get baseline log mean for each cell, before adding in perturbation effects
      log_mu <- predict(mod, newdata = cell_covariates)
      
      # For each perturbation effect configuration
      for (effect_idx in seq_len(nrow(perturbation_effect_sizes))) {
        non_pert_effect <- perturbation_effect_sizes$non_pert_effect[effect_idx]
        pert_effect <- perturbation_effect_sizes$pert_effect[effect_idx]
        theta_new <- perturbation_effect_sizes$theta[effect_idx]

        # Simulate num_grnas_per_perturbation gRNAs for this target/effect combination
        for (grna_idx in seq_len(num_grnas_per_perturbation)) {

          # Generate perturbation indicators (which cells are perturbed)
          is_perturbed <- rbinom(num_cells, size = 1, prob = prob_pert)

          # Simulate counts: add effect to perturbed and non-perturbed cells
          log_mu_pert <- log_mu + pert_effect * is_perturbed + non_pert_effect * (1 - is_perturbed)
          simulated_counts <- rnegbin(
            n = num_cells,
            mu = exp(log_mu_pert),
            theta = theta_new
          )

          simulated_rows[[length(simulated_rows) + 1]] <- row_to_sparse_matrix(simulated_counts)
          perturbation_rows[[length(perturbation_rows) + 1]] <- row_to_sparse_matrix(is_perturbed)

          # Create row name with all three effect parameters
          row_name <- paste0(sim_grna_name, "_probpert_", prob_pert,
                            "_nonpert_", round(non_pert_effect, 4),
                            "_pert_", round(pert_effect, 4),
                            "_theta_", round(theta_new, 4),
                            "_grna_", grna_idx)
          row_names <- c(row_names, row_name)
        }
      }
      cat("Progress: prob_pert", prob_idx, "/", length(prob_perturbed),
          "- gRNA target", target_idx, "/", num_targets, "done.\n")
    }
  } # End prob_perturbed loop
  
  # Combine the sparse row vectors into sparse matrices (rows = gRNAs, cols = cells)
  simulated_matrix <- do.call(rbind, simulated_rows)
  perturbation_indicators <- do.call(rbind, perturbation_rows)
  
  # Set row and column names for both matrices
  rownames(simulated_matrix) <- row_names
  rownames(perturbation_indicators) <- row_names
  
  # Preserve cell names if they exist
  if (!is.null(colnames(grna_matrix_rows))) {
    colnames(simulated_matrix) <- colnames(grna_matrix_rows)
    colnames(perturbation_indicators) <- colnames(grna_matrix_rows)
  }
  
  # Return list with simulated matrix, ground truth indicators, and model information
  return(list(
    simulated_matrix = simulated_matrix,
    perturbation_indicators = perturbation_indicators,
    model_info = model_info,
    sim_params = list(formula = as.character(formula) |> paste(collapse=" "),
                      perturbation_effect_sizes = perturbation_effect_sizes,
                      prob_perturbed = prob_perturbed)
  ))
}

plot_umi_hist_by_pert <- function(grna_matrix_rows, sim_results, grna_idx, sim_idx) {
  
  grna_real <- grna_matrix_rows[grna_idx,]
  grna_sim <- sim_results$simulated_matrix[sim_idx,]
  # make sure this real grna was used for this sim grna
  is_pert <- sim_results$perturbation_indicators[sim_idx, ]
  sim_name <- rownames(sim_results$simulated_matrix)[sim_idx]
  
  grepl(pattern = rownames(grna_matrix_rows)[grna_idx], x = sim_name) |>
    stopifnot()
  stopifnot(length(grna_real) == length(grna_sim) && length(grna_real) == length(is_pert))
  
  logp1_trans <- function(base = 10) {
    trans <- function(x) log(x + 1, base)
    inv   <- function(x) base^x - 1  # correct inverse
    
    # Major breaks: at base^k - 1
    brks <- function(limits) {
      rng <- limits + 1
      scales::log_breaks(base = base)(rng) - 0
    }
    
    scales::trans_new(
      name = paste0("logp1-", base),
      transform = trans,
      inverse   = inv,
      breaks    = brks,
      domain = c(0, Inf)
    )
  }
  
  # parse sim name
  sim_name_pieces <- strsplit(sim_name, "_")[[1]]
  keep <- sim_name_pieces[which(sim_name_pieces == "probpert"):(which(sim_name_pieces == "grna")-1)]
  sim_name_plot <- paste0(keep[c(1,3,5,7)], "=", keep[c(2,4,6,8)], collapse=" ") 
  
  plt1 <- data.frame(
    umi = c(grna_real, grna_sim),
    group = rep(c("real", "sim"), each = length(grna_real))
  ) |>
    ggplot(aes(x = umi, fill = group)) +
    geom_histogram(bins=30) +
    facet_wrap(~group) +
    theme_bw()  +
    scale_y_continuous(trans = logp1_trans(10)) +
    scale_x_continuous(trans = logp1_trans(10)) +
    ggtitle(paste0("Real vs Sim: ", sim_name_plot))
  
  plt2 <- data.frame(
    umi = grna_sim,
    pert_status = ifelse(is_pert == 1, "pert", "not pert")
  ) |>
    ggplot(aes(x = umi, fill = pert_status)) +
    geom_histogram(bins=30) +
    facet_wrap(~pert_status) +
    theme_bw()  +
    scale_y_continuous(trans = logp1_trans(10)) +
    scale_x_continuous(trans = logp1_trans(10)) +
    ggtitle("Simulated UMI counts by pert status")
  
  cowplot::plot_grid(plt1, plt2, nrow = 2)
}

