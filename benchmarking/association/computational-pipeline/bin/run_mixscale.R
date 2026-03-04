#!/usr/bin/env Rscript
# Method: mixscale
# Task: computational benchmarking in association testing
# NOTE: low MOI only!
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - assignments.rds: Named vector mapping cell IDs to original pseudo-target names
#
# Outputs:
# - association_computational_mixscale.csv: All perturbation-gene test results

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id <- args[2]

cat("Running mixscale computational analysis\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Dataset ID:", dataset_id, "\n")

# Load required packages
cat("Loading packages...\n")
library(Seurat)
library(Mixscale)
library(dplyr)

# Set Seurat options
options(Seurat.object.assay.version = 'v3')

# Load data
cat("Loading data...\n")
response_mat <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
assign_vec <- readRDS(file.path(dataset_dir, "assignments.rds"))

NT_NAME = "non-targeting"
if(! NT_NAME %in% assign_vec) {
  stop("'non-targeting' is not present in 'assignments.rds'!")
}



cat("  Response matrix:", nrow(response_mat), "genes x", ncol(response_mat), "cells\n")
cat("  assignments:", length(assign_vec), "cells\n")

# Critical: Remove cells with zero total counts
cell_umi_totals <- Matrix::colSums(response_mat)
if(any(cell_umi_totals == 0)) {
  n_zero <- sum(cell_umi_totals == 0)
  cat("WARNING: Removing", n_zero, "cells with zero UMI counts\n")
  response_mat <- response_mat[, cell_umi_totals > 0]
  assign_vec <- assign_vec[colnames(response_mat)]
  cat("  After filtering:", ncol(response_mat), "cells remain\n")
}
all_targets = setdiff(as.character(assign_vec), NT_NAME)
start_time <- Sys.time()

# Create base Seurat object with preprocessing (will update $gene in each loop iteration)
cat("Creating base Seurat object and preprocessing...\n")
seu <- CreateSeuratObject(
  counts = response_mat,
  assay  = "RNA",
  min.cells = 0,
  min.features = 0
)

DefaultAssay(seu) <- "RNA"

seu <- seu |>
  NormalizeData() |>
  FindVariableFeatures() |> # only 2000 most variable genes get scaled and used for perturb scores
  ScaleData() |>
  RunPCA(verbose = FALSE)

cat("  Preprocessing complete\n")

  seu$gene <- factor(assign_vec[colnames(seu)])

  # Calculate perturbation signatures (curr_target vs "nt")
  seu <- CalcPerturbSig(
    object        = seu,
    assay         = "RNA",
    slot          = "data",
    gd.class      = "gene",
    nt.cell.class = NT_NAME,
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
    nt.class.name   = NT_NAME,
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
    nt.class.name   = NT_NAME,
    PRTB_list       = all_targets, 
    logfc.threshold = 0,
    min.cells.group = 0,
    min.pct         = 0,
    verbose         = FALSE
  )

  de_res <- do.call(Run_wmvRegDE, de_args)

  for(j in seq_along(de_res)) {
    de_res[[j]] <- mutate(
      de_res[[j]],
        grna_target = names(de_res)[j],
        response_id = rownames(de_res[[j]])
    )
  }


end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total time:", round(elapsed, 2), "seconds\n")

# Combine all results (remove NULL entries if any)
cat("Combining results...\n")
results <- do.call(rbind, de_res)
cat("  Total pairs:", nrow(results), "\n")

# Save results
cat("Saving results...\n")
write.csv(results, "association_computational_mixscale.csv", row.names = FALSE)

cat("Mixscale computational benchmarking complete!\n")
cat("Results saved to: association_computational_mixscale.csv\n")
