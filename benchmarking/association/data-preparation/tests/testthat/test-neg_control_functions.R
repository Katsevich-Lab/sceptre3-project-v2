# test-make_neg_control.R
# Unit tests for neg-control specific functions

library(testthat)
library(Matrix)
library(dplyr)

# Source utility functions and neg-control functions
source("../../utils_data_prep.R")
source("../../neg_control_functions.R")


# Test select_genes_random() ----------------------------------------------

test_that("select_genes_random samples correct number", {
  # Create mock matrix with rownames (simulates ODM)
  mock_odm <- matrix(0, nrow = 100, ncol = 10)
  rownames(mock_odm) <- paste0("GENE_", 1:100)

  result <- select_genes_random(mock_odm, num_genes = 10, random_seed = 42)

  expect_length(result, 10)
  expect_true(all(result %in% paste0("GENE_", 1:100)))
})

test_that("select_genes_random is reproducible", {
  mock_odm <- matrix(0, nrow = 50, ncol = 5)
  rownames(mock_odm) <- paste0("GENE_", 1:50)

  result1 <- select_genes_random(mock_odm, num_genes = 5, random_seed = 123)
  result2 <- select_genes_random(mock_odm, num_genes = 5, random_seed = 123)

  expect_equal(result1, result2)
})

test_that("select_genes_random handles requesting more than available", {
  mock_odm <- matrix(0, nrow = 10, ncol = 5)
  rownames(mock_odm) <- paste0("GENE_", 1:10)

  expect_warning(
    result <- select_genes_random(mock_odm, num_genes = 20, random_seed = 42),
    "Requested .* genes but only .* available"
  )

  expect_length(result, 10)  # Should return all available
})


# Test create_pseudo_targets() --------------------------------------------

test_that("create_pseudo_targets creates sequential names", {
  guides <- c("NT_guide_1", "NT_guide_2", "NT_guide_3")
  result <- create_pseudo_targets(guides)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_equal(result$grna_id, guides)
  expect_equal(result$grna_target, c("non-targeting1", "non-targeting2", "non-targeting3"))
})

test_that("create_pseudo_targets has no duplicates", {
  guides <- c("guide_A", "guide_B", "guide_C", "guide_D")
  result <- create_pseudo_targets(guides)

  expect_false(any(duplicated(result$grna_target)))
})


# Test identify_nt_only_cells() -------------------------------------------

test_that("identify_nt_only_cells finds correct cells", {
  # Create test assignment matrix: 5 guides × 6 cells
  # Guides: 2 NT, 3 targeting
  # Cells: 1-2 have NT only, 3-4 have targeting, 5-6 have mix
  scep_assn_mat <- as(Matrix::Matrix(c(
    1, 0, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1,
    0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0
  ), nrow = 5, byrow = TRUE, sparse = TRUE), "nMatrix")
  rownames(scep_assn_mat) <- c("NT_guide_1", "NT_guide_2",
                                 "targeting_guide_1", "targeting_guide_2", "targeting_guide_3")

  grna_target_df <- data.frame(
    grna_id = c("NT_guide_1", "NT_guide_2", "targeting_guide_1", "targeting_guide_2", "targeting_guide_3"),
    grna_target = c("non-targeting", "non-targeting", "gene1", "gene2", "gene3"),
    stringsAsFactors = FALSE
  )

  result <- identify_nt_only_cells(scep_assn_mat, grna_target_df)

  # Should return cells 1 and 2 only (NT only)
  expect_equal(result, c(1, 2, 6))
})

test_that("identify_nt_only_cells excludes ignored targets", {
  scep_assn_mat <- Matrix::Matrix(c(
    1, 0, 1,  # NT_guide (non-targeting)
    0, 1, 0,  # off_target (should be ignored)
    0, 0, 1   # targeting_guide (gene1)
  ), nrow = 3, byrow = TRUE, sparse = TRUE) |>
    as("generalMatrix") |>
    as("nMatrix")
  rownames(scep_assn_mat) <- c("NT_guide", "off_target", "targeting_guide")

  grna_target_df <- data.frame(
    grna_id = c("NT_guide", "off_target", "targeting_guide"),
    grna_target = c("non-targeting", "nt_off_target", "gene1"),
    stringsAsFactors = FALSE
  )

  result <- identify_nt_only_cells(
    scep_assn_mat,
    grna_target_df,
    targets_to_ignore = "nt_off_target"
  )

  # Should return only cell 1 (cell 2 has ignored target)
  expect_equal(result, 1)
})


# Test filter_all_objects() -----------------------------------------------

test_that("filter_all_objects removes empty cells and genes", {
  # Create test data with some all-zero rows/columns
  response_matrix <- Matrix::Matrix(c(
    10, 0, 5, 0,   # gene1: has expression
    0,  0, 0, 0,   # gene2: all zeros (should be removed)
    6,  0, 3, 0    # gene3: has expression
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  rownames(response_matrix) <- c("gene1", "gene2", "gene3")

  grna_indicator <- Matrix::Matrix(c(
    1, 0, 1, 0,  # guide1
    0, 0, 0, 0,  # guide2: all zeros (should be removed)
    0, 0, 0, 1   # guide3
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  rownames(grna_indicator) <- c("guide1", "guide2", "guide3")

  grna_target_df <- data.frame(
    grna_id = c("guide1", "guide2", "guide3"),
    grna_target = c("target1", "target2", "target3"),
    stringsAsFactors = FALSE
  )

  cell_idx <- c(10, 20, 30, 40)

  result <- filter_all_objects(response_matrix, grna_indicator, grna_target_df, cell_idx)

  # Should remove gene2, guide2, and cells with zero UMI (cells 2 and 4)
  expected_response_matrix <- Matrix::Matrix(c(10,5,6,3), nrow = 2, byrow = TRUE, sparse = TRUE) |>
    `rownames<-`(c("gene1", "gene3")) |>
    as("dgCMatrix")
  expected_grna_indicator  <- Matrix::Matrix(c(
    1, 1
  ), nrow = 1, byrow = TRUE, sparse = TRUE) |>
    `rownames<-`(c("guide1"))
  expected_grna_target_df <- data.frame(
    grna_id = c("guide1"),
    grna_target = c("target1"),
    stringsAsFactors = FALSE
  )
  
  expect_equal(result$response_matrix, expected_response_matrix)  # gene1, gene3
  expect_equal(result$grna_indicator, expected_grna_indicator)   # guide1, guide3
  expect_equal(result$grna_target_df, expected_grna_target_df)   # target1, target3
  expect_equal(result$cell_idx, c(10, 30))       # cells 10, 30
  expect_equal(result$genes, c("gene1", "gene3"))
})


# Test prepare_cell_metadata() --------------------------------------------

test_that("prepare_cell_metadata creates correct structure", {
  grna_indicator <- Matrix::Matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  rownames(grna_indicator) <- c("guide1", "guide2", "guide3")

  grna_target_df <- data.frame(
    grna_id = c("guide1", "guide2", "guide3"),
    grna_target = c("target1", "target2", "target3"),
    stringsAsFactors = FALSE
  )

  cell_covariates <- data.frame(
    response_n_nonzero = c(50, 60, 70, 80, 90),
    response_n_umis = c(500, 600, 700, 800, 900),
    grna_n_nonzero = c(1, 1, 1, 1, 1),
    grna_n_umis = c(10, 12, 14, 16, 18),
    batch = c("A", "A", "B", "B", "A")
  )
  rownames(cell_covariates) <- paste0("cell_", 1:5)

  cell_idx <- c(1, 3, 5)

  result <- prepare_cell_metadata(grna_indicator, grna_target_df, cell_covariates, cell_idx)

  # Check cell_info structure
  expect_equal(nrow(result$cell_info), 3)
  expect_true(all(c("cell_idx", "expressed_guide", "grna_target", "cell_name") %in% names(result$cell_info)))
  expect_equal(result$cell_info$cell_idx, c(1, 3, 5))
  expect_equal(result$cell_info$grna_target, c("target1", "target2", "target3"))

  # Check cell_covariates_subset structure
  expect_equal(nrow(result$cell_covariates_subset), 3)
  expect_true(all(c("response_n_nonzero_full", "response_n_umis_full",
                    "grna_n_nonzero_full", "grna_n_umis_full", "batch") %in%
                    names(result$cell_covariates_subset)))
})

