#!/usr/bin/env Rscript
# Positive-control runner for the ORIGINAL 2021 SCEPTRE method
# (Katsevich-Lab/sceptre-manuscript). Tests each target against its OWN gene, so the
# p-values measure the method's power to detect a true effect.
#
# The 2021 package has no dedicated power-check routine, but a power check is simply
# the ordinary association test applied to on-target pairs: we construct the
# positive-control pairs (every target that is itself a measured gene, paired with
# that gene) and run the standard distilled CRT on them. This mirrors what the modern
# sceptre positive-control run does with construct_positive_control_pairs().
#
# The five-step procedure (four precomputations then the pairwise test) is unchanged
# from the published version; see the negative-control script for the full outline.
#
# PARALLELISM. This pipeline is free to use all allocated cores (unlike the
# computational pipeline, which is deliberately single-CPU for comparable runtimes).
# The gene precompute, gRNA precompute and pair loop run under mclapply over NCPUS
# cores; the dispersion regularization is a reduce and stays serial. This does not
# change results: the published test sets its own seed at the start of each pair, so
# output is identical to a serial run (verified byte-for-byte).
#
# HIGH MOI ONLY. Low-MOI data is handled by sceptre_v030.
#
# NO QC: the 2021 package performs no filtering of any kind, and the modern
# positive-control run likewise skips QC, so all methods see the same cells and genes.
#
# Inputs (from dataset_dir): response_matrix.rds, grna_matrix.rds,
# grna_target_data_frame.csv, cell_covariates.csv. (No formula_object.rds or
# discovery_pairs.rds in this pipeline -- the formula is built here and the pairs are
# constructed below, matching the modern positive-control run.)
#
# Output: association_on_target_sceptre_manuscript.csv

suppressPackageStartupMessages({ library(Matrix); library(sceptre); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
dataset_dir <- args[1]
dataset_id  <- args[2]

B                     <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_B", "500"))
regularization_amount <- as.numeric(Sys.getenv("SCEPTRE_MANUSCRIPT_REG", "3"))
side                  <- Sys.getenv("SCEPTRE_MANUSCRIPT_SIDE", "both")
seed                  <- as.integer(Sys.getenv("SCEPTRE_MANUSCRIPT_SEED", "1234"))
n_cores               <- as.integer(Sys.getenv("NCPUS", "1"))

cat("Running MANUSCRIPT sceptre (positive control)\n")
cat("  dataset_dir:", dataset_dir, "\n  dataset_id:", dataset_id, "\n")
cat("  B:", B, " regularization_amount:", regularization_amount, " side:", side,
    " seed:", seed, " cores:", n_cores, "\n")

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

# Covariate formula, matching the modern positive-control run for this dataset.
formula_object <- ~ log(response_n_nonzero_full + 1) + log(response_n_umis_full + 1) +
                    log(grna_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + prep_batch
cat("  formula:", deparse(formula_object), "\n")
mm <- model.matrix(formula_object, data = cell_covariates)
covariate_matrix <- as.data.frame(mm[, colnames(mm) != "(Intercept)", drop = FALSE])

# Positive-control pairs: every target that is itself a measured gene, paired with it.
pc_targets <- intersect(unique(as.character(grna_target_df$grna_target)), rownames(response_matrix))
pairs <- data.frame(grna_target = pc_targets, response_id = pc_targets, stringsAsFactors = FALSE)

cat("Data loaded:\n")
cat("  response_matrix:", nrow(response_matrix), "genes x", ncol(response_matrix), "cells\n")
cat("  grna_matrix:", nrow(grna_matrix), "gRNAs x", ncol(grna_matrix), "cells\n")
cat("  positive control pairs:", nrow(pairs), "\n")
if (nrow(pairs) == 0L) stop("No positive control pairs: no grna_target is a row of response_matrix.")
stopifnot(ncol(response_matrix) == ncol(grna_matrix),
          ncol(response_matrix) == nrow(cell_covariates),
          nrow(covariate_matrix) == ncol(response_matrix))

genes   <- unique(pairs$response_id)
targets <- unique(pairs$grna_target)

start_time <- Sys.time()

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

cat("Gene precomputation round 1 over", length(genes), "genes...\n")
r1 <- par_lapply(genes, function(g) {
  expr <- as.numeric(response_matrix[g, ])
  list(size = run_gene_precomputation(expr, covariate_matrix, gene_precomp_size = NULL)[["gene_precomp_size"]],
       gmean = if (regularization_amount > 0) log_geom_mean(expr) else NA_real_)
})
sizes_unreg <- vapply(r1, function(z) z$size, numeric(1)); names(sizes_unreg) <- genes
log_gmeans  <- vapply(r1, function(z) z$gmean, numeric(1)); names(log_gmeans)  <- genes

if (regularization_amount > 0) {
  cat("Regularizing gene sizes (amount =", regularization_amount, ")...\n")
  sizes_reg <- regularize_thetas(genes_log_gmean = log_gmeans, theta = sizes_unreg,
                                 bw_adjust = regularization_amount, plot_me = FALSE)
  names(sizes_reg) <- names(sizes_unreg)
} else {
  sizes_reg <- sizes_unreg
}

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

write.csv(results, "association_on_target_sceptre_manuscript.csv", row.names = FALSE)

elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Manuscript sceptre positive control complete!\n")
cat("Results saved to: association_on_target_sceptre_manuscript.csv\n")
cat("Total pairs analyzed:", nrow(results), "\n")
cat("Total time:", round(elapsed, 2), "seconds\n")
