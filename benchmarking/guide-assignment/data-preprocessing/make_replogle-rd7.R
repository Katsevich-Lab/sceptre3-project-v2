this_file <- normalizePath(sys.frame(1)$ofile)
script_dir <- dirname(this_file)
source(file.path(script_dir, "convert_odm_to_h5ad.R"))
source("~/.Rprofile")



make_data <- function(data_dir, output_base_dir,
                     make_crispat=FALSE, make_pertpy=FALSE,
                     make_sceptre=FALSE, make_cleanser=FALSE,
                     sizes = NULL, use_real_responses = FALSE) {
  
  # Setup Python/reticulate if needed for h5ad
  if(make_crispat || make_pertpy) {
    library(reticulate)
    conda_create("r-anndata", packages = c("python=3.12"))
    use_condaenv("r-anndata", required = TRUE)
    py_install(c("numpy", "scipy", "h5py", "anndata"), 
               envname = "r-anndata", pip = TRUE)
    py_config()
  }
  
  # PHASE 1: Load and process gRNA matrix
  grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(data_dir, "grna.odm"))
  cell_names <- paste0("CELL_", 1:ncol(grna_odm))
  print("Loading gRNA matrix...")
  grna_mat <- odm_to_R(grna_odm)  # SLOW - do once
  print("gRNA matrix loaded")
  colnames(grna_mat) <- cell_names
  
  # Store original dimensions before clearing from memory
  original_num_grnas <- nrow(grna_mat)
  original_num_cells <- ncol(grna_mat)
  
  # Loop through sizes for gRNA-based outputs
  for(i in 1:nrow(sizes)) {
    name <- sizes[i, "name"]
    num_grnas <- min(sizes[i, "num_grnas"], nrow(grna_mat))
    num_cells <- min(sizes[i, "num_cells"], ncol(grna_mat))
    
    print(paste("Processing", name, "- gRNAs:", num_grnas, "cells:", num_cells))
    
    # Subset gRNA matrix
    grna_subset <- grna_mat[1:num_grnas, 1:num_cells]
    output_dir <- file.path(output_base_dir, name)
    
    # Write crispat
    if(make_crispat) {
      crispat_dir <- file.path(output_dir, "crispat")
      dir.create(crispat_dir, recursive=TRUE, showWarnings=FALSE)
      print(paste("  Writing crispat to", crispat_dir))
      R_to_h5ad(grna_subset, file.path(crispat_dir, "grna_matrix.h5ad"))
    }
    
    # Write pertpy
    if(make_pertpy) {
      pertpy_dir <- file.path(output_dir, "pertpy")
      dir.create(pertpy_dir, recursive=TRUE, showWarnings=FALSE)
      print(paste("  Writing pertpy to", pertpy_dir))
      R_to_h5ad(grna_subset, file.path(pertpy_dir, "grna_matrix.h5ad"))
    }
    
    # Write cleanser
    if(make_cleanser) {
      cleanser_dir <- file.path(output_dir, "cleanser")
      dir.create(cleanser_dir, recursive=TRUE, showWarnings=FALSE)
      print(paste("  Writing cleanser to", cleanser_dir))
      Matrix::writeMM(grna_subset, file.path(cleanser_dir, "grna_matrix.mtx"))
    }
    
    # Write sceptre gRNA matrix
    if(make_sceptre) {
      sceptre_dir <- file.path(output_dir, "sceptre")
      dir.create(sceptre_dir, recursive=TRUE, showWarnings=FALSE)
      print(paste("  Writing sceptre gRNA matrix to", sceptre_dir))
      saveRDS(grna_subset, file.path(sceptre_dir, "grna_matrix.rds"))
    }
  }
  
  # Clear gRNA matrix from memory
  rm(grna_mat)
  gc()
  print("gRNA matrix cleared from memory")
  
  # PHASE 2: If sceptre needed, process response matrix
  if(make_sceptre) {
    # Load grna_target_data_frame once (needed for both real and fake responses)
    scep <- readRDS(file.path(data_dir, "sceptre_object.rds"))
    grna_target_df <- scep@grna_target_data_frame

    if(use_real_responses) {
      print("Processing with real response matrix...")
      # Load response ODM metadata (not the data itself)
      response_odm <- ondisc::initialize_odm_from_backing_file(file.path(data_dir, "gene.odm"))

      # Find the maximum number of responses we'll need across all sizes
      max_responses_needed <- max(sizes$num_responses[!is.infinite(sizes$num_responses)])
      if(length(max_responses_needed) == 0 || is.na(max_responses_needed)) {
          # All are Inf, we need to load everything (but this shouldn't happen in practice)
          max_responses_needed <- nrow(response_odm)
          warning("All num_responses are Inf - loading entire response matrix!")
      }

      # Load only the rows we need
      print(paste("Loading response matrix - first", max_responses_needed, "rows only..."))
      response_mat <- odm_to_R(response_odm, num_rows = max_responses_needed)
      print(paste("Response matrix loaded:", nrow(response_mat), "x", ncol(response_mat)))
      colnames(response_mat) <- cell_names

      # Loop through sizes for response outputs
      for(i in 1:nrow(sizes)) {
        name <- sizes[i, "name"]
        num_grnas <- min(sizes[i, "num_grnas"], original_num_grnas)  # Using stored dimension
        num_responses <- min(sizes[i, "num_responses"], nrow(response_mat))
        num_cells <- min(sizes[i, "num_cells"], original_num_cells)

        print(paste("Processing", name, "- responses:", num_responses, "cells:", num_cells))

        # Subset response matrix
        response_subset <- response_mat[1:num_responses, 1:num_cells]

        # Subset grna_target_data_frame to match the gRNAs
        grna_target_subset <- grna_target_df[1:num_grnas, ]

        # Write to sceptre directory
        sceptre_dir <- file.path(output_base_dir, name, "sceptre")
        print(paste("  Writing sceptre response matrix to", sceptre_dir))
        saveRDS(response_subset, file.path(sceptre_dir, "response_matrix.rds"))
        write.csv(grna_target_subset, file.path(sceptre_dir, "grna_target_data_frame.csv"), row.names=FALSE)
      }

      # Clear response matrix from memory
      rm(response_mat)
      gc()
      print("Response matrix cleared from memory")

    } else {
      print("Creating minimal fake response matrix for guide assignment only...")

      # Loop through sizes and create fake response matrices
      for(i in 1:nrow(sizes)) {
        name <- sizes[i, "name"]
        num_grnas <- min(sizes[i, "num_grnas"], original_num_grnas)
        num_cells <- min(sizes[i, "num_cells"], original_num_cells)

        print(paste("Processing", name, "- fake responses: 1 gene x", num_cells, "cells"))

        # Create minimal sparse matrix: 1 gene x num_cells, all zeros
        fake_response_mat <- Matrix::sparseMatrix(
            i = integer(0),  # No non-zero entries
            j = integer(0),
            x = numeric(0),
            dims = c(1L, num_cells),
            dimnames = list(
                "FAKE_GENE_1",  # Single fake gene
                cell_names[1:num_cells]  # Use same cell names as grna_matrix
            )
        )

        # Subset grna_target_data_frame to match the gRNAs
        grna_target_subset <- grna_target_df[1:num_grnas, ]

        # Write to sceptre directory
        sceptre_dir <- file.path(output_base_dir, name, "sceptre")
        print(paste("  Writing fake sceptre response matrix to", sceptre_dir))
        saveRDS(fake_response_mat, file.path(sceptre_dir, "response_matrix.rds"))
        write.csv(grna_target_subset, file.path(sceptre_dir, "grna_target_data_frame.csv"), row.names=FALSE)
      }
    }
  }
  
  print("Data generation complete!")
}



sizes <- data.frame(
  name = c("replogle-rd7_small", "replogle-rd7_medium", "replogle-rd7"),
  num_grnas = c(500, 2500,  Inf),
  num_cells = c(5000, 50000, Inf),
  # num_responses = c(NA, NA,  Inf),  # not needed for fake data
  stringsAsFactors = FALSE
)


# paths
gasperini_base_fp <- .get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR")
data_dir <- file.path(gasperini_base_fp, "processed/rd7")
output_base_dir <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data")

make_data(
  data_dir = data_dir,
  output_base_dir = output_base_dir,
  make_crispat = FALSE,
  make_pertpy = FALSE,
  make_cleanser = FALSE,
  make_sceptre = TRUE,
  sizes = sizes,
  use_real_responses = FALSE  
)
