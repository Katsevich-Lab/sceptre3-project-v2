#!/usr/bin/env Rscript
# Method: mixscale
# Task: negative control analysis in association testing
# NOTE: low MOI only!
#
# Strategy: Loop over mixscale_nt_guides, using each in turn as the NT control group
# Try Run_wmvRegDE first (with mixscale scores), fall back to Run_stdDE if it fails
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - assignments.rds: Named vector mapping cell IDs to pseudo-target names
# - mixscale_nt_guides.rds: Vector of pseudo-target names to loop over as NT controls
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
assign_vec <- readRDS(file.path(dataset_dir, "assignments.rds"))
mixscale_nt_targets <- readRDS(file.path(dataset_dir, "mixscale_nt_targets.rds"))

cat("  Response matrix:", nrow(response_mat), "genes x", ncol(response_mat), "cells\n")
cat("  Assignments:", length(assign_vec), "cells\n")
cat("  Mixscale pseudo-targets to loop over:", length(mixscale_nt_targets), "\n")

# Critical: Remove cells with zero total counts
cell_umi_totals <- Matrix::colSums(response_mat)
if(any(cell_umi_totals == 0)) {
  n_zero <- sum(cell_umi_totals == 0)
  cat("WARNING: Removing", n_zero, "cells with zero UMI counts\n")
  response_mat <- response_mat[, cell_umi_totals > 0]
  assign_vec <- assign_vec[colnames(response_mat)]
  cat("  After filtering:", ncol(response_mat), "cells remain\n")
}

# Create base Seurat object with preprocessing
cat("Creating base Seurat object and preprocessing...\n")
seu_base <- CreateSeuratObject(
  counts = response_mat,
  assay  = "RNA",
  min.cells = 0,
  min.features = 0
)

seu_base$gene <- factor(assign_vec[colnames(seu_base)])
DefaultAssay(seu_base) <- "RNA"

seu_base <- seu_base |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA(verbose = FALSE)

cat("  Preprocessing complete\n")

# Loop over each pseudo-target, using it as the control group (parallelized)
cat("Looping over", length(mixscale_nt_targets), "pseudo-targets...\n")
start_time <- Sys.time()

all_results <- future_lapply(seq_along(mixscale_nt_targets), function(nt_idx) {
  curr_nt_target <- mixscale_nt_targets[nt_idx]
  cat("  Processing pseudo-target", nt_idx, "/", length(mixscale_nt_targets), ":", curr_nt_target, "\n")

  seu <- seu_base
  my_perturbations <- setdiff(levels(seu$gene), curr_nt_target)

  # Calculate perturbation signatures (always works)
  seu <- CalcPerturbSig(
    object        = seu,
    assay         = "RNA",
    slot          = "data",
    gd.class      = "gene",
    nt.cell.class = curr_nt_target,
    reduction     = "pca",
    num.neighbors = 20,
    new.assay.name = "PRTB"
  )

  # Run Mixscale (always works)
  seu <- RunMixscale(
    object          = seu,
    assay           = "PRTB",
    slot            = "scale.data",
    labels          = "gene",
    nt.class.name   = curr_nt_target,
    de.assay        = "RNA",
    max.de.genes    = 100,
    new.class.name  = "mixscale_score",
    fine.mode       = FALSE,
    verbose         = FALSE
  )

  # Prepare arguments for DE testing (same for both functions)
  de_args <- list(
    object          = seu,
    assay           = "RNA",
    slot            = "counts",
    labels          = "gene",
    nt.class.name   = curr_nt_target,
    PRTB_list       = my_perturbations,
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

  # Extract all results for this control group
  lapply(names(de_res), function(target_name) {
    de_res[[target_name]] |>
      dplyr::mutate(
        pseudo_target = target_name,
        nt_control = curr_nt_target,
        gene = rownames(de_res[[target_name]])
      )
  }) |> do.call(what = rbind)
}, future.seed = TRUE) |>
  setNames(mixscale_nt_targets)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total loop time:", round(elapsed, 2), "seconds\n")

# Combine all results
cat("Combining results...\n")
results <- do.call(rbind, all_results)
cat("  Total pairs:", nrow(results), "\n")

# Save results
cat("Saving results...\n")
write.csv(results, "association_neg_control_mixscale.csv", row.names = FALSE)

cat("Mixscale negative control analysis complete!\n")
cat("Results saved to: association_neg_control_mixscale.csv\n")
