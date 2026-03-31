#!/usr/bin/env Rscript
# Method: mixscale
# Task: positive control analysis in association testing
# NOTE: low MOI only!
#
# Inputs (from dataset_dir):
# - response_matrix.rds: Gene expression sparse matrix (genes x cells)
# - assignments.rds: Named vector mapping cell IDs to target gene names
#
# Outputs:
# - results_mixscale.rds: Raw output from Run_wmvRegDE (list of data frames)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id <- args[2]

cat("Running mixscale positive control analysis\n")
cat("Dataset directory:", dataset_dir, "\n")
cat("Dataset ID:", dataset_id, "\n")

# Load required packages
cat("Loading packages...\n")
# library(Matrix)
library(Seurat)
library(Mixscale)

# Set Seurat options
options(Seurat.object.assay.version = 'v3')

# Step 1: Load data and filter zero-count cells
cat("Loading data...\n")
response_mat <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
assign_vec <- readRDS(file.path(dataset_dir, "assignments.rds"))

NT_NAME = "non-targeting"
if(! NT_NAME %in% assign_vec) {
  stop("'non-targeting' is not present in 'assignments.rds'!")
}

cat("  Response matrix:", nrow(response_mat), "genes x", ncol(response_mat), "cells\n")
cat("  Assignments:", length(assign_vec), "cells\n")

# Critical: Remove cells with zero total counts (causes Run_wmvRegDE to fail)
cell_umi_totals <- Matrix::colSums(response_mat)
if(any(cell_umi_totals == 0)) {
  n_zero <- sum(cell_umi_totals == 0)
  cat("WARNING: Removing", n_zero, "cells with zero UMI counts\n")
  response_mat <- response_mat[, cell_umi_totals > 0]
  assign_vec <- assign_vec[colnames(response_mat)]
  cat("  After filtering:", ncol(response_mat), "cells remain\n")
}

# Step 2: Create Seurat object with assignments
cat("Creating Seurat object...\n")
seu <- CreateSeuratObject(
  counts = response_mat,
  assay  = "RNA",
  min.cells = 0,
  min.features = 0
)

# Add perturbation assignments as metadata
seu$gene <- factor(assign_vec[colnames(seu)])
DefaultAssay(seu) <- "RNA"

cat("  Unique perturbations:", length(unique(seu$gene)), "\n")

# Step 3: Standard Seurat preprocessing
cat("Running preprocessing...\n")
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, features = VariableFeatures(seu))
seu <- RunPCA(seu, verbose=FALSE)

# Step 4: Calculate perturbation signatures
cat("Calculating perturbation signatures...\n")
seu <- CalcPerturbSig(
  object        = seu,
  assay         = "RNA",
  slot          = "data",
  gd.class      = "gene",
  nt.cell.class = "non-targeting",
  reduction     = "pca",
  num.neighbors = 20,
  new.assay.name = "PRTB"
)

# Step 5: Run Mixscale
cat("Running Mixscale...\n")
seu <- RunMixscale(
  object          = seu,
  assay           = "PRTB",
  slot            = "scale.data",
  labels          = "gene",
  nt.class.name   = "non-targeting",
  de.assay        = "RNA",
  max.de.genes    = 100,
  new.class.name  = "mixscale_score",
  fine.mode       = FALSE,
  verbose         = FALSE
)

# Check how many targets got mixscale scores
n_with_scores <- sum(tapply(seu$mixscale_score, seu$gene,
                            function(x) sum(x != 0 & x != 1)) > 0)
cat("  Targets with non-trivial mixscale scores:", n_with_scores, "\n")

# Step 6: Run weighted regression DE test
cat("Running weighted regression DE test...\n")
my_perturbations <- setdiff(levels(seu$gene), "non-targeting")
cat("  Testing", length(my_perturbations), "perturbations\n")

# TODO this very inefficient! Can we do on-targets only?
# it seems that massive logfc.threshold may skip all except
# the initially DE genes for each target, or at least for
# the targets with at least 5 initially DE genes
# But i'm not totally sure on the behavior of that so 
# I'm hesitant to try
warning("NOTE: currently computing every target-gene pair, not just on-targest!!!")

# prtb_scores <- Tool(seu, slot = "RunMixscale")
# targets_with_scores <- names(prtb_scores)

# lapply(my_perturbations, function(gene) {
#   # for the genes that got mixscale scores, a massive logfc.threshold
#   # filters out every gene from the test EXCEPT the on-target, which
#   # is guaranteed no matter logfc.threshold
#   # this saves us from computing many p-values when we only need on-targets
#   logfc_thresh <- ifelse(gene %in% targets_with_scores, 9999, 0)
# })


de_res <- Run_wmvRegDE(
  object          = seu,
  assay           = "RNA",
  slot            = "counts",
  labels          = "gene",
  nt.class.name   = "non-targeting",
  PRTB_list       = my_perturbations,
  # logfc.threshold = 0,
  # min.cells.group = 0,
  # min.pct         = 0,
  verbose         = FALSE
)

results <- lapply(my_perturbations, function(gene) de_res[[gene]][gene,]) |>
  do.call(what=rbind)

# Step 7: Save raw output
cat("Saving results...\n")
write.csv(results, "association_on_target_mixscale.csv")

cat("Mixscale on-target association testing complete!\n")
cat("Results saved to: association_on_target_mixscale.csv\n")
