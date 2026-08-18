#!/usr/bin/env Rscript
# Computational-pipeline runner for the ORIGINAL 2021 SCEPTRE method
# (Katsevich-Lab/sceptre-manuscript). Runs the published distilled-CRT + skew-t
# method over the discovery pairs and reports its runtime.
#
# The statistical procedure is unchanged from the published version: a dCRT test
# with a NegBin test statistic, calibrated by fitting a skew-t distribution to the
# resampled null statistics. It proceeds in five steps, with the first 4 being precomputations:
#   (1) gene precompute, round 1  For each gene, fit NB(expression ~ covariates) and
#                                 record the dispersion and the log
#                                 geometric mean of expression.
#
#   (2) regularize the sizes      Shrink the per-gene dispersions (the sctransform approach).
#
#   (3) gene precompute, round 2  Refit each gene at its regularized dispersion to get
#                                 the per-cell fitted-mean offsets
#
#   (4) gRNA precompute           For each target, fit logistic(perturbed ~ covariates)
#                                 to get each cell's probability of being perturbed;
#                                 this is the null distribution the test resamples from.
#
#   (5) pairwise test             For each (target, gene) pair, run the distilled CRT
#                                 and return a p-value
#
# For computational benchmarking, we do not parallelize, so this is run with the
# precomputations as a sequential pure-R loop. The total time, including precomputations,
# is the reported method runtime.
#
# WHERE THE MANUSCRIPT FUNCTIONS COME FROM
#   The manuscript package is itself named `sceptre` (version 1.0.0) and exports the
#   functions used below (run_gRNA_precomputation, run_gene_precomputation,
#   run_sceptre_using_precomp, regularize_thetas, log_geom_mean). This script runs
#   inside the sceptre-manuscript apptainer image, which installs ONLY that package,
#   so `library(sceptre)` loads the manuscript functions and nothing else.
#
# Two modeling choices, made so results are comparable to modern SCEPTRE and the
# other methods:
#  - A union gRNA strategy is used rather than singleton.
#  - The test is two-sided.
#
# Inputs (from dataset_dir):
# - response_matrix.rds          genes x cells sparse dgCMatrix (rownames = gene ids)
# - grna_matrix.rds              gRNAs x cells binary matrix (rownames = grna ids)
# - grna_target_data_frame.csv   grna_id -> grna_target mapping
# - cell_covariates.csv          per-cell covariates (rows aligned to matrix columns)
# - formula_object.rds           covariate formula string
# - discovery_pairs.rds          grna_target x response_id pairs to test
#
# Output:
# - association_computational_sceptre_manuscript.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre) })

# Build the covariate design matrix that the manuscript code consumes via `~ .`.
# We expand the covariate formula with model.matrix() so the manuscript method
# adjusts for exactly the same covariates as modern sceptre, then drop the intercept
# column (glm/glm.nb/vglm add their own).
build_covariate_matrix <- function(cell_covariates, formula_object) {
  mm <- model.matrix(formula_object, data = cell_covariates)
  keep <- colnames(mm) != "(Intercept)"
  as.data.frame(mm[, keep, drop = FALSE])
}

# For each unique grna_target, form the UNION indicator across the gRNAs that
# target it, plus the logistic gRNA precomputation (per-cell probabilities).
build_target_precomp <- function(targets, grna_matrix, grna_target_df, covariate_matrix) {
  target_to_grnas <- split(as.character(grna_target_df$grna_id),
                           as.character(grna_target_df$grna_target))
  indicators <- gRNA_precomps <- vector("list", length(targets)) |> setNames(targets)
  for (t in targets) {
    ids <- intersect(target_to_grnas[[t]], rownames(grna_matrix))
    if (length(ids) == 0L) stop("No gRNAs in grna_matrix for target ", t)
    ind <- as.integer(Matrix::colSums(grna_matrix[ids, , drop = FALSE]) > 0)
    indicators[[t]]    <- ind
    gRNA_precomps[[t]] <- run_gRNA_precomputation(ind, covariate_matrix)
  }
  list(indicators = indicators, gRNA_precomps = gRNA_precomps)
}

# Run the full manuscript pipeline (steps 1-5) over the (grna_target, response_id)
# pairs, in memory and serially. Returns a data.frame of results.
run_old_sceptre_pairs <- function(pairs, response_matrix, grna_matrix, grna_target_df,
                                  covariate_matrix, B, seed, side, regularization_amount) {
  genes   <- unique(as.character(pairs$response_id))
  targets <- unique(as.character(pairs$grna_target))
  stopifnot(all(genes %in% rownames(response_matrix)))

  # (4) gRNA precompute for every target (amortized across genes)
  cat("gRNA precomputation over", length(targets), "targets...\n")
  tp <- build_target_precomp(targets, grna_matrix, grna_target_df, covariate_matrix)

  # (1) gene precompute round 1: unregularized size + log geom mean per gene
  cat("Gene precomputation round 1 over", length(genes), "genes...\n")
  sizes_unreg <- numeric(length(genes)); names(sizes_unreg) <- genes
  log_gmeans  <- numeric(length(genes)); names(log_gmeans)  <- genes
  for (g in genes) {
    expr <- as.numeric(response_matrix[g, ])
    sizes_unreg[[g]] <- run_gene_precomputation(expr, covariate_matrix, gene_precomp_size = NULL)[["gene_precomp_size"]]
    if (regularization_amount > 0) log_gmeans[[g]] <- log_geom_mean(expr)
  }

  # (2) regularize sizes across genes (sctransform-style shrinkage)
  if (regularization_amount > 0) {
    cat("Regularizing gene sizes (amount =", regularization_amount, ")...\n")
    sizes_reg <- regularize_thetas(genes_log_gmean = log_gmeans, theta = sizes_unreg,
                                   bw_adjust = regularization_amount, plot_me = FALSE)
    names(sizes_reg) <- names(sizes_unreg)
  } else {
    sizes_reg <- sizes_unreg
  }

  # (3) gene precompute round 2 + (5) pair dCRT, gene-outer / target-inner so each
  # gene's offsets are computed once and reused across its paired targets.
  pairs_by_gene <- split(as.character(pairs$grna_target), as.character(pairs$response_id))
  cat("Gene precomputation round 2 + pairwise dCRT (B =", B, ", side =", side, ")...\n")
  results <- vector("list", nrow(pairs)); ri <- 0L
  for (g in genes) {
    expr <- as.numeric(response_matrix[g, ])
    offsets_g <- run_gene_precomputation(expr, covariate_matrix,
                                         gene_precomp_size = sizes_reg[[g]])[["gene_precomp_offsets"]]
    size_g <- sizes_reg[[g]]
    for (t in pairs_by_gene[[g]]) {
      res <- run_sceptre_using_precomp(
        expressions          = expr,
        gRNA_indicators      = tp$indicators[[t]],
        gRNA_precomp         = tp$gRNA_precomps[[t]],
        gene_precomp_size    = size_g,
        gene_precomp_offsets = offsets_g,
        B = B, seed = seed, side = side, reduced_output = TRUE)
      ri <- ri + 1L
      results[[ri]] <- data.frame(grna_target = t, response_id = g, res,
                                  stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, results)
}

# ----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]

# Tunables (env-overridable so the Nextflow module can set them)
B                     <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_B", "500"))
regularization_amount <- as.numeric(Sys.getenv("SCEPTRE_MANUSCRIPT_REG", "3"))
side                  <- Sys.getenv("SCEPTRE_MANUSCRIPT_SIDE", "both")
seed                  <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_SEED", "1234"))

cat("Running MANUSCRIPT sceptre (computational)\n")
cat("  dataset_dir:", dataset_dir, "\n  dataset_id:", dataset_id, "\n")
cat("  B:", B, " regularization_amount:", regularization_amount, " side:", side, " seed:", seed, "\n")

# --- load inputs ---
response_matrix <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna_matrix     <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
grna_target_df  <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(dataset_dir, "cell_covariates.csv"), stringsAsFactors = FALSE)
formula_object  <- as.formula(readRDS(file.path(dataset_dir, "formula_object.rds")))
discovery_pairs <- readRDS(file.path(dataset_dir, "discovery_pairs.rds"))

cat("  formula:", deparse(formula_object), "\n")
covariate_matrix <- build_covariate_matrix(cell_covariates, formula_object)
pairs <- data.frame(grna_target = as.character(discovery_pairs$grna_target),
                    response_id = as.character(discovery_pairs$response_id),
                    stringsAsFactors = FALSE)

cat("Data loaded:\n")
cat("  response_matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  grna_matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  covariate_matrix cols:", paste(colnames(covariate_matrix), collapse = ", "), "\n")
cat("  pairs:", nrow(pairs), "\n")
stopifnot(ncol(response_matrix) == ncol(grna_matrix),
          ncol(response_matrix) == nrow(cell_covariates),
          nrow(covariate_matrix) == ncol(response_matrix))

# --- run (timed: entire precompute + test, for fair runtime accounting) ---
start_time <- Sys.time()
results <- run_old_sceptre_pairs(pairs, response_matrix, grna_matrix, grna_target_df,
                                 covariate_matrix, B = B, seed = seed, side = side,
                                 regularization_amount = regularization_amount)
end_time <- Sys.time()

write.csv(results, "association_computational_sceptre_manuscript.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Manuscript sceptre complete!\n")
cat("Results saved to: association_computational_sceptre_manuscript.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
