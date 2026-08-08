make_data <- function() {


  # Parameters
  num_cells <- 500
  num_grnas <- 100
  num_NTs <- 10
  num_targets <- 10  # should divide num_grnas
  num_responses <- 100
  num_batches <- 4  # must divide num_cells

  set.seed(42)  # For reproducibility

  # Create gRNA and cell names
  grna_names <- c(paste0("gRNA_", seq_len(num_grnas)),
                  paste0("NT_", seq_len(num_NTs)))
  cell_names <- sprintf("CELL_%04d", seq_len(num_cells))

  # 1. Simulate gRNA matrix (gRNAs x cells)
  grna_matrix <- matrix(0, nrow = num_grnas + num_NTs, ncol = num_cells)

  grnas_expressed <- sample(num_grnas + num_NTs, num_cells, replace = TRUE)

  for (i in seq_len(num_cells)) {
    grna_matrix[grnas_expressed[i], i] <- sample(10:15, 1)  # Primary expression

    background_counts <- sample(1:8, num_grnas + num_NTs - 1, TRUE)
    background_mask <- rbinom(num_grnas + num_NTs - 1, 1, 0.5)

    grna_matrix[-grnas_expressed[i], i] <- background_counts * background_mask
  }

  colnames(grna_matrix) <- cell_names
  rownames(grna_matrix) <- grna_names

  # 2. Simulate response matrix (responses x cells)
  response_names <- paste0("response_", seq_len(num_responses))
  response_matrix <- matrix(rpois(num_responses * num_cells, lambda = 3),
                            nrow = num_responses, ncol = num_cells)
  rownames(response_matrix) <- response_names
  colnames(response_matrix) <- cell_names

  # 3. Create gRNA-to-target mapping
  grnas_per_target <- num_grnas %/% num_targets
  targets <- c(
    rep(paste0("target_", 1:num_targets), each = grnas_per_target),
    rep("non-targeting", num_NTs)
  )


  grna_target_df <- data.frame(
    grna_id = grna_names,
    grna_target = targets,
    stringsAsFactors = FALSE
  )

  # 4. create covariate_data_frame
  covariate_data_frame <- data.frame(
    batch = paste0("batch_", rep(1:num_batches, each = num_cells / num_batches))
  )

  # 4. return
  list(grna_matrix=grna_matrix,
       response_matrix=response_matrix,
       grna_target_data_frame = grna_target_df,
       covariate_data_frame = covariate_data_frame)
}
