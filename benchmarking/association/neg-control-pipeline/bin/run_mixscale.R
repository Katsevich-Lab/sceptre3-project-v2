#!/usr/bin/env Rscript
# Method: mixscale
# Task: negative control analysis in association testing
# NOTE: low MOI only!
#
# Strategy: Complement-style analysis - for each pseudo-target, test it against all other NT guides
# For each target i: assign cells as either target_i or "nt" (all others), then test target_i vs "nt"
# Try Run_wmvRegDE first (with mixscale scores), fall back to Run_stdDE if it fails
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - assignments.rds: Named vector mapping cell IDs to original pseudo-target names
#   (the pseudo-targets to loop over are just unique(assignments) -- we test the
#    full cartesian product, so no separate target list file is needed)
#
# Outputs:
# - association_neg_control_mixscale.csv: All perturbation-gene test results

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id <- args[2]

cat("Running mixscale negative control analysis\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Dataset ID:", dataset_id, "\n")

# Load required packages
cat("Loading packages...\n")
library(Seurat)
library(Mixscale)
library(future)
library(future.apply)

# Set Seurat options
options(Seurat.object.assay.version = 'v3')

# Remove size limit for future globals (for large datasets with parallel processing)
options(future.globals.maxSize = +Inf)

# Configure parallelization based on NCPUS
n_processors <- as.integer(Sys.getenv("NCPUS", "1"))
cat("Processors allocated:", n_processors, "\n")

if(n_processors > 1) {
  cat("Setting up parallel processing with", n_processors, "workers\n")
  plan(multicore, workers = n_processors)
} else {
  cat("Running in serial mode\n")
  plan(sequential)
}

# Load data
cat("Loading data...\n")
response_mat <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
assign_vec_original <- readRDS(file.path(dataset_dir, "assignments.rds"))
# We always test every pseudo-target (full cartesian product), so the loop set
# is exactly the distinct targets present in the assignments.
all_targets <- sort(unique(assign_vec_original))

cat("  Response matrix:", nrow(response_mat), "genes x", ncol(response_mat), "cells\n")
cat("  Original assignments:", length(assign_vec_original), "cells\n")
cat("  Total pseudo-targets to test:", length(all_targets), "\n")

# Validate that all targets are present in assignments
stopifnot(all(all_targets %in% assign_vec_original))

# Critical: Remove cells with zero total counts
cell_umi_totals <- Matrix::colSums(response_mat)
if(any(cell_umi_totals == 0)) {
  n_zero <- sum(cell_umi_totals == 0)
  cat("WARNING: Removing", n_zero, "cells with zero UMI counts\n")
  response_mat <- response_mat[, cell_umi_totals > 0]
  assign_vec_original <- assign_vec_original[colnames(response_mat)]
  cat("  After filtering:", ncol(response_mat), "cells remain\n")
}

# Create base Seurat object with preprocessing (will update $gene in each loop iteration)
cat("Creating base Seurat object and preprocessing...\n")
seu_base <- CreateSeuratObject(
  counts = response_mat,
  assay  = "RNA",
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu_base) <- "RNA"

seu_base <- seu_base |>
  NormalizeData() |>
  FindVariableFeatures() |> # only 2000 most variable genes get scaled and used for perturb scores
  ScaleData(features = VariableFeatures(seu_base)) |>
  RunPCA(verbose = FALSE)

cat("  Preprocessing complete\n")

# Loop over each pseudo-target: test target_i vs all other NT guides (parallelized)
cat("Looping over", length(all_targets), "pseudo-targets (complement-style analysis)...\n")
start_time <- Sys.time()

all_results <- future_lapply(seq_along(all_targets), function(target_idx) {
  curr_target <- all_targets[target_idx]
  cat("  Processing pseudo-target", target_idx, "/", length(all_targets), ":", curr_target, "\n")

  # Create binary assignment: current target vs all others ("nt")
  assign_vec_binary <- ifelse(assign_vec_original == curr_target, curr_target, "nt")
  names(assign_vec_binary) <- names(assign_vec_original)

  # Clone base Seurat object and set binary assignments
  seu <- seu_base
  seu$gene <- factor(assign_vec_binary[colnames(seu)])

  # Calculate perturbation signatures (curr_target vs "nt")
  seu <- CalcPerturbSig(
    object        = seu,
    assay         = "RNA",
    slot          = "data",
    gd.class      = "gene",
    nt.cell.class = "nt",
    reduction     = "pca",
    num.neighbors = 20,
    new.assay.name = "PRTB"
  )

  # Run Mixscale
  seu <- RunMixscale(
    object          = seu,
    assay           = "PRTB",
    slot            = "scale.data",
    labels          = "gene",
    nt.class.name   = "nt",
    de.assay        = "RNA",
    max.de.genes    = 100,
    new.class.name  = "mixscale_score",
    fine.mode       = FALSE,
    verbose         = FALSE
  )

  # Prepare arguments for DE testing (only test curr_target, not "nt")
  de_args <- list(
    object          = seu,
    assay           = "RNA",
    slot            = "counts",
    labels          = "gene",
    nt.class.name   = "nt",
    PRTB_list       = curr_target,  # Only test the single target
    logfc.threshold = 0,
    min.cells.group = 0,
    min.pct         = 0,
    verbose         = FALSE
  )

  # Try Run_wmvRegDE first, fall back to Run_stdDE
  de_res <- tryCatch({
    cat("    Attempting Run_wmvRegDE...\n")
    do.call(Run_wmvRegDE, de_args)
  }, error = function(e) {
    cat("    Run_wmvRegDE failed, falling back to Run_stdDE\n")
    cat("    Error was:", conditionMessage(e), "\n")
    do.call(Mixscale:::Run_stdDE, de_args)
  })

  # Extract results for curr_target only (de_res should have only one element)
  if(curr_target %in% names(de_res)) {
    de_res[[curr_target]] |>
      dplyr::mutate(
        grna_target = curr_target,
        response_id = rownames(de_res[[curr_target]])
      )
  } else {
    warning("No results returned for target: ", curr_target)
    NULL
  }
}, future.seed = TRUE)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total loop time:", round(elapsed, 2), "seconds\n")

# Combine all results (remove NULL entries if any)
cat("Combining results...\n")
all_results <- all_results[!sapply(all_results, is.null)]
results <- do.call(rbind, all_results)
cat("  Total pairs:", nrow(results), "\n")

# Save results
cat("Saving results...\n")
write.csv(results, "association_neg_control_mixscale.csv", row.names = FALSE)

cat("Mixscale negative control analysis complete!\n")
cat("Results saved to: association_neg_control_mixscale.csv\n")
