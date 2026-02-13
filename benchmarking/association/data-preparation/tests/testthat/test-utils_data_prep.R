# test-utils_data_prep.R
# Unit tests for shared utility functions in utils_data_prep.R

library(testthat)
library(Matrix)

# Source the utility functions
source("../../utils_data_prep.R")


# Test enforce_single_guide_per_cell() ------------------------------------

test_that("enforce_single_guide_per_cell removes multi-guide cells", {
  # Create test matrix: 3 guides × 4 cells
  # Cell 1: guide1, Cell 2: guide2, Cell 3: guide1+guide2 (multi), Cell 4: guide2+guide3 (multi)
  mat <- Matrix::Matrix(c(
    1, 0, 1, 0,  # guide1
    0, 1, 1, 1,  # guide2
    0, 0, 0, 1,   # guide3
    0, 0, 1, 1,
    0, 0, 1, 1
  ), nrow = 5, byrow = TRUE, sparse = TRUE)
  rownames(mat) <- c("guide1", "guide2", "guide3", "guide4", "guide5")

  result <- enforce_single_guide_per_cell(mat, random_seed = 123)
  set.seed(123)
  col3_is_1 <- which(mat[,3] == 1)
  col4_is_1 <- which(mat[,4] == 1)
  result_manual <- mat
  result_manual[setdiff(col3_is_1, sample(col3_is_1, 1)), 3] <- 0
  result_manual[setdiff(col4_is_1, sample(col4_is_1, 1)), 4] <- 0
  
    # All cells should have exactly 1 guide
  expect_equal(as.vector(Matrix::colSums(result)), c(1, 1, 1, 1))
  expect_equal(result, result_manual)

  # Result should be sparse
  expect_s4_class(result, "sparseMatrix")
})


test_that("enforce_single_guide_per_cell handles already-low-MOI data", {
  # Matrix where each cell already has exactly 1 guide
  mat <- Matrix::Matrix(c(
    1, 0, 0,  # guide1
    0, 1, 0,  # guide2
    0, 0, 1   # guide3
  ), nrow = 3, byrow = TRUE, sparse = TRUE)

  result <- enforce_single_guide_per_cell(mat, random_seed = 123)

  # Should be unchanged
  expect_equal(result, mat)
})


# Test prepare_frperturb_covariates() ------------------------------------

test_that("prepare_frperturb_covariates transforms specified columns", {
  covariates <- data.frame(
    response_n_umis = c(1000, 2000, 1500),
    grna_n_umis = c(50, 100, 75),
    batch = c("A", "B", "A")
  )

  grna_targets <- c("target1", "target2", "target3")
  covariates_to_log1p <- c("response_n_umis", "grna_n_umis")

  result <- prepare_frperturb_covariates(
    cell_covariates = covariates,
    grna_targets = grna_targets,
    covariates_to_log1p = covariates_to_log1p
  )

  # Should have log1p-transformed columns
  expect_true("response_n_umis_log1p" %in% names(result))
  expect_true("grna_n_umis_log1p" %in% names(result))

  # Should NOT have original untransformed columns
  expect_false("response_n_umis" %in% names(result))
  expect_false("grna_n_umis" %in% names(result))

  # Batch should remain untransformed
  expect_true("batch" %in% names(result))
  expect_equal(result$batch, c("A", "B", "A"))

  # Should have perturbation column
  expect_true("perturbation" %in% names(result))
  expect_equal(result$perturbation, grna_targets)

  # Check transformation is correct
  expect_equal(result$response_n_umis_log1p, log1p(c(1000, 2000, 1500)))
})

test_that("prepare_frperturb_covariates errors on mismatched lengths", {
  covariates <- data.frame(
    response_n_umis = c(1000, 2000, 1500)
  )

  grna_targets <- c("target1", "target2")  # Length mismatch

  expect_error(
    prepare_frperturb_covariates(
      cell_covariates = covariates,
      grna_targets = grna_targets,
      covariates_to_log1p = "response_n_umis"
    )
  )
})

test_that("prepare_frperturb_covariates errors on invalid column names", {
  covariates <- data.frame(
    response_n_umis = c(1000, 2000, 1500)
  )

  grna_targets <- c("target1", "target2", "target3")

  expect_error(
    prepare_frperturb_covariates(
      cell_covariates = covariates,
      grna_targets = grna_targets,
      covariates_to_log1p = "nonexistent_column"
    )
  )
})


# Test odm_to_sparse_matrix() ---------------------------------------------
# Can test with regular matrices since [ operator works the same way

test_that("odm_to_sparse_matrix extracts correct subset", {
  # Create test matrix: 5 genes × 10 cells
  test_mat <- matrix(c(
    0, 5, 0, 10, 0, 0, 3, 0, 0, 0,   # gene1
    0, 0, 0,  0, 2, 0, 0, 0, 0, 4,   # gene2
    1, 0, 0,  0, 0, 0, 0, 0, 0, 0,   # gene3
    0, 0, 3,  0, 0, 2, 0, 0, 0, 0,   # gene4
    0, 0, 0,  0, 0, 0, 0, 1, 0, 0    # gene5
  ), nrow = 5, byrow = TRUE)
  rownames(test_mat) <- paste0("GENE_", 1:5)

  # Extract genes 1, 3, 5 and cells 2, 4, 6, 8
  genes <- c("GENE_1", "GENE_3", "GENE_5")
  cell_idx <- c(2, 4, 6, 8)

  result <- odm_to_sparse_matrix(test_mat, genes, cell_idx, set_rownames = TRUE)

  # Check dimensions
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 4)

  # Check it's sparse
  expect_s4_class(result, "sparseMatrix")

  # Check rownames
  expect_equal(rownames(result), genes)

  # Check values
  expect_equal(result[1, 1] |> as.numeric(), 5)   # GENE_1, cell 2
  expect_equal(result[1, 2] |> as.numeric(), 10)  # GENE_1, cell 4
  expect_equal(result[2, 1] |> as.numeric(), 0)   # GENE_3, cell 2
  expect_equal(result[3, 4] |> as.numeric(), 1)   # GENE_5, cell 8
})

test_that("odm_to_sparse_matrix handles all-zero genes", {
  test_mat <- matrix(c(
    0, 0, 0, 0,  # gene1: all zeros
    1, 2, 3, 4   # gene2: has values
  ), nrow = 2, byrow = TRUE)
  rownames(test_mat) <- c("GENE_1", "GENE_2")

  genes <- c("GENE_1", "GENE_2")
  cell_idx <- c(1, 2, 3, 4)

  result <- odm_to_sparse_matrix(test_mat, genes, cell_idx)

  # Should still create correct dimensions
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 4)

  # First row should be all zeros
  expect_equal(sum(result[1, ]), 0)
})

test_that("odm_to_sparse_matrix optionally skips rownames", {
  test_mat <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  rownames(test_mat) <- c("GENE_1", "GENE_2")

  result <- odm_to_sparse_matrix(test_mat, c("GENE_1", "GENE_2"), c(1, 2), set_rownames = FALSE)

  expect_null(rownames(result))
})
