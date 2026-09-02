#!/usr/bin/env Rscript
# Negative-control runner for the ORIGINAL 2021 SCEPTRE method
# (Katsevich-Lab/sceptre-manuscript). Tests the pipeline's negative-control pairs,
# whose p-values should be approximately uniform if the method is calibrated.
#
# The statistical procedure is unchanged from the published version: a distilled
# conditional randomization test with a negative-binomial test statistic, calibrated
# by fitting a skew-t distribution to the resampled null statistics. It proceeds in
# five steps; the first four are "precomputations" -- quantities computed once and
# reused across many pairs:
#   (1) gene precompute, round 1  for each gene, fit NB(expression ~ covariates) and
#                                 record the dispersion (size / theta) and the log
#                                 geometric mean of expression
#   (2) regularize the sizes      shrink the per-gene dispersions toward a smooth
#                                 trend across all genes (the sctransform approach)
#   (3) gene precompute, round 2  refit each gene at its regularized dispersion to get
#                                 the per-cell fitted-mean offsets
#   (4) gRNA precompute           for each target, fit logistic(perturbed ~ covariates)
#                                 to get each cell's probability of being perturbed
#   (5) pairwise test             for each (target, gene) pair, run the distilled CRT
#
# PARALLELISM. Unlike the computational pipeline (which is deliberately single-CPU so
# runtimes are comparable), this pipeline is free to use all allocated cores. Steps
# (1), (4) and (3)+(5) are embarrassingly parallel and are run with mclapply over
# NCPUS cores; step (2) is a reduce across all genes and is inherently serial.
# Parallelising here is safe and does NOT change the results: the published test
# function sets its own seed at the start of every pair, so each pair's p-value is
# independent of execution order. Output is therefore identical to a serial run.
#
# HIGH MOI ONLY. Low-MOI data is handled by sceptre_v030; this errors out on a
# low-MOI dataset rather than silently running the wrong analysis.
#
# NO QC: the 2021 package performs no cell- or gene-level filtering of any kind, so
# every method in the benchmark sees exactly the same cells and genes.
#
# Inputs (from dataset_dir): response_matrix.rds, grna_matrix.rds,
# grna_target_data_frame.csv, cell_covariates.csv, formula_object.rds,
# discovery_pairs.rds (the negative-control pairs).
#
# Output: association_neg_control_sceptre_manuscript.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]

B                     <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_B", "500"))
regularization_amount <- as.numeric(Sys.getenv("SCEPTRE_MANUSCRIPT_REG", "3"))
side                  <- Sys.getenv("SCEPTRE_MANUSCRIPT_SIDE", "both")
seed                  <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_SEED", "1234"))
n_cores               <- as.integer(Sys.getenv("NCPUS", "1"))

cat("Running MANUSCRIPT sceptre (negative control)\n")
cat("  dataset_dir:", dataset_dir, "\n  dataset_id:", dataset_id, "\n")
cat("  B:", B, " regularization_amount:", regularization_amount, " side:", side,
    " seed:", seed, " cores:", n_cores, "\n")

# Parallel map over independent units of work; falls back to lapply on one core.
par_lapply <- function(X, FUN) {
  if (n_cores <= 1) return(lapply(X, FUN))
  res <- parallel::mclapply(X, FUN, mc.cores = n_cores, mc.preschedule = TRUE)
  failed <- vapply(res, function(z) inherits(z, "try-error"), logical(1))
  if (any(failed)) {
    stop("parallel worker failed: ",
         conditionMessage(attr(res[[which(failed)[1]]], "condition")))
  }
  res
}

# --- MOI guard: this method is high-MOI only ---
DATASET_NAMES <- c("gasperini", "replogle")
dataset_name <- DATASET_NAMES[sapply(DATASET_NAMES, function(n) grepl(n, dataset_id, ignore.case = TRUE))]
if (length(dataset_name) != 1) stop("Could not determine dataset from dataset_id: ", dataset_id)
moi <- list(gasperini = "high", replogle = "low")[[dataset_name]]
if (moi != "high") {
  stop("sceptre_manuscript is a HIGH-MOI method, but dataset '", dataset_id, "' is ", moi,
       " MOI. Use sceptre_v030 for low-MOI data.")
}
cat("  detected dataset:", dataset_name, " MOI:", moi, "\n")

# --- load inputs ---
response_matrix <- readRDS(file.path(dataset_dir, "response_matrix.rds"))
grna_matrix     <- readRDS(file.path(dataset_dir, "grna_matrix.rds"))
grna_target_df  <- read.csv(file.path(dataset_dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
cell_covariates <- read.csv(file.path(dataset_dir, "cell_covariates.csv"), stringsAsFactors = FALSE)
formula_object  <- as.formula(readRDS(file.path(dataset_dir, "formula_object.rds")))
discovery_pairs <- readRDS(file.path(dataset_dir, "discovery_pairs.rds"))

# Design matrix from the covariate formula (log transforms live in the formula), with
# the intercept dropped since glm/glm.nb/vglm add their own.
mm <- model.matrix(formula_object, data = cell_covariates)
covariate_matrix <- as.data.frame(mm[, colnames(mm) != "(Intercept)", drop = FALSE])

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

genes   <- unique(pairs$response_id)
targets <- unique(pairs$grna_target)
stopifnot(all(genes %in% rownames(response_matrix)))

start_time <- Sys.time()

# (4) gRNA precompute: union indicator per target + logistic assignment probabilities
cat("gRNA precomputation over", length(targets), "targets...\n")
target_to_grnas <- split(as.character(grna_target_df$grna_id),
                         as.character(grna_target_df$grna_target))
tp <- par_lapply(targets, function(t) {
  ids <- intersect(target_to_grnas[[t]], rownames(grna_matrix))
  if (length(ids) == 0L) stop("No gRNAs in grna_matrix for target ", t)
  ind <- as.integer(Matrix::colSums(grna_matrix[ids, , drop = FALSE]) > 0)
  list(ind = ind, prob = run_gRNA_precomputation(ind, covariate_matrix))
})
names(tp) <- targets

# (1) gene precompute round 1: unregularized dispersion + log geometric mean
cat("Gene precomputation round 1 over", length(genes), "genes...\n")
r1 <- par_lapply(genes, function(g) {
  expr <- as.numeric(response_matrix[g, ])
  list(size = run_gene_precomputation(expr, covariate_matrix, gene_precomp_size = NULL)[["gene_precomp_size"]],
       gmean = if (regularization_amount > 0) log_geom_mean(expr) else NA_real_)
})
sizes_unreg <- vapply(r1, function(z) z$size, numeric(1)); names(sizes_unreg) <- genes
log_gmeans  <- vapply(r1, function(z) z$gmean, numeric(1)); names(log_gmeans)  <- genes

# (2) regularize dispersions across genes -- a reduce, inherently serial
if (regularization_amount > 0) {
  cat("Regularizing gene sizes (amount =", regularization_amount, ")...\n")
  sizes_reg <- regularize_thetas(genes_log_gmean = log_gmeans, theta = sizes_unreg,
                                 bw_adjust = regularization_amount, plot_me = FALSE)
  names(sizes_reg) <- names(sizes_unreg)
} else {
  sizes_reg <- sizes_unreg
}

# (3) gene precompute round 2 + (5) pairwise dCRT, parallel over genes so each gene's
# offsets are computed once and reused across the targets it is paired with.
pairs_by_gene <- split(pairs$grna_target, pairs$response_id)
cat("Gene precomputation round 2 + pairwise dCRT (B =", B, ", side =", side, ")...\n")
res_list <- par_lapply(genes, function(g) {
  expr <- as.numeric(response_matrix[g, ])
  offsets_g <- run_gene_precomputation(expr, covariate_matrix,
                                       gene_precomp_size = sizes_reg[[g]])[["gene_precomp_offsets"]]
  do.call(rbind, lapply(pairs_by_gene[[g]], function(t) {
    res <- run_sceptre_using_precomp(
      expressions = expr, gRNA_indicators = tp[[t]]$ind, gRNA_precomp = tp[[t]]$prob,
      gene_precomp_size = sizes_reg[[g]], gene_precomp_offsets = offsets_g,
      B = B, seed = seed, side = side, reduced_output = TRUE)
    data.frame(grna_target = t, response_id = g, res, stringsAsFactors = FALSE)
  }))
})
results <- do.call(rbind, res_list)

end_time <- Sys.time()

write.csv(results, "association_neg_control_sceptre_manuscript.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Manuscript sceptre negative control complete!\n")
cat("Results saved to: association_neg_control_sceptre_manuscript.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
