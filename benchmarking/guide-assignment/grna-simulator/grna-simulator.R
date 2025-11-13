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
    perturbation_effect_sizes = log(3), # effect size of an actual perturbation
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
    prob_perturbed >= 0 && prob_perturbed <= 1,
    inherits(formula, "formula"),
    length(formula) == 2  # One-sided formula
  )
  
  # Construct full formula internally (response name doesn't matter to user)
  full_formula <- as.formula(paste("..y.. ~", as.character(formula)[2]))
  
  # Storage for simulated gRNA rows, perturbation indicators, and model info
  simulated_rows <- list()
  perturbation_rows <- list()
  row_names <- character()
  model_info <- list()
  
  # For each target (input gRNA row)
  for (target_idx in seq_len(num_targets)) {
    
    # Get count vector for this target
    count_vector <- grna_matrix_rows[target_idx, ]
    
    # Fit negative binomial model
    fit_data <- cbind(..y.. = count_vector, cell_covariates)
    mod <- glm.nb(full_formula, data = fit_data)
    
    # Get target name (or use index if no rownames)
    sim_grna_name <- if (!is.null(rownames(grna_matrix_rows))) {
      paste0("sim_", rownames(grna_matrix_rows)[target_idx])
    } else {
      paste0("sim_grna_", target_idx)
    }
    
    # sim_grna_name <- paste0("sim_grna_", target_idx)
    
    # Store model information
    model_info[[sim_grna_name]] <- list(
      coefficients = summary(mod)$coef,
      theta = mod$theta
    )
    
    # Get baseline log mean for each cell, before adding in perturbation effects
    log_mu <- predict(mod, newdata = cell_covariates)
    
    # For each perturbation effect size
    for (effect_idx in seq_along(perturbation_effect_sizes)) {
      effect_size <- perturbation_effect_sizes[effect_idx]
      
      # Simulate num_grnas_per_perturbation gRNAs for this target/effect combination
      for (grna_idx in seq_len(num_grnas_per_perturbation)) {
        
        # Generate perturbation indicators (which cells are perturbed)
        is_perturbed <- rbinom(num_cells, size = 1, prob = prob_perturbed)
        
        # Simulate counts: add effect_size to log_mu for perturbed cells
        simulated_counts <- rnegbin(
          n = num_cells,
          mu = exp(log_mu + effect_size * is_perturbed),
          theta = mod$theta
        )
        
        # Store simulated counts and perturbation indicators
        # simulated_rows[[length(simulated_rows) + 1]] <- simulated_counts
        # perturbation_rows[[length(perturbation_rows) + 1]] <- is_perturbed
        
        simulated_rows[[length(simulated_rows) + 1]] <- row_to_sparse_matrix(simulated_counts)
        perturbation_rows[[length(perturbation_rows) + 1]] <- row_to_sparse_matrix(is_perturbed)
        
        # Create row name
        row_name <- sprintf(
          "%s_effect_%.3f_grna_%d",
          sim_grna_name,
          effect_size,
          grna_idx
        )
        row_names <- c(row_names, row_name)
      }
    }
    cat("Progress: gRNA target", target_idx, "/", num_targets, "done.\n")
  }

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
    model_info = model_info
  ))
}