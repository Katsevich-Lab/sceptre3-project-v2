#!/usr/bin/env Rscript
# Read-only check of the full-Gasperini inputs for the assignment -> positive-
# control analysis, before any dataset is built. Writes nothing but its log.
#
#   SIF=<pipeline-images>/sceptre.sif \
#     sbatch -c 2 -t 1:00:00 ../../guide-assignment/grna-simulator/run-in-image.sh check-inputs.R
#
# Checks (each reported PASS/FAIL; all run even if one fails):
#   1. both ODMs have 207,324 cells, in the same order as the sceptre object's
#      covariate rows (the ODMs store no cell names, so order is checked by
#      recomputing per-cell UMI and nonzero totals from each ODM)
#   2. grna.odm's rows are the sceptre object's target table guides; that table
#      matches the sceptre/grna_target_data_frame.csv copy; 377 ENSG targets, 754 guides
#   3. all 377 ENSG targets are rows of response.odm
#   4. cleanser/grna_matrix.mtx (fishash's input) equals grna.odm entry for entry,
#      in the same row and column order
#   5. fishashplus/grna_matrix.rds (fishash+'s input) equals grna.odm entry for
#      entry; reports its dimnames
#   6. the sceptre object's covariates equal sceptre/_full_cov_mat.csv
#   7. the fishash output's indices are within the .mtx dimensions

suppressPackageStartupMessages({
  library(Matrix); library(sceptre); library(ondisc)
})
source("~/.Rprofile")                                    # .get_config_path()

in_dir   <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "input_data", "gasperini")
pipe_dir <- file.path(in_dir, "sceptre-pipeline")
fishash_csv <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment", "outputs",
                         "gb_c1_full", "assignments_fishash_gasperini.csv")

results <- character(0)
check <- function(name, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) { cat("  error:", conditionMessage(e), "\n"); FALSE })
  results[name] <<- if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", results[name], name))
}
as_dgc <- function(m) drop0(as(as(m, "CsparseMatrix"), "generalMatrix"))
same_matrix <- function(a, b) {
  identical(dim(a), dim(b)) && identical(a@p, b@p) && identical(a@i, b@i) &&
    identical(as.numeric(a@x), as.numeric(b@x))
}

# ---- load ------------------------------------------------------------------
cat("reading sceptre object + ODMs from", pipe_dir, "\n")
scep <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp    = file.path(pipe_dir, "sceptre_object.rds"),
  response_odm_file_fp = file.path(pipe_dir, "response.odm"),
  grna_odm_file_fp     = file.path(pipe_dir, "grna.odm")
)
response_odm <- initialize_odm_from_backing_file(file.path(pipe_dir, "response.odm"))
grna_odm     <- initialize_odm_from_backing_file(file.path(pipe_dir, "grna.odm"))
cov_df    <- scep@covariate_data_frame
target_df <- scep@grna_target_data_frame
# The ODMs store no cell names (colnames() is empty), so cell counts come from
# the length of one row, and cell order is checked in 1b/1c by recomputing the
# per-cell totals from the ODMs and comparing them with the covariate rows.
n_c_grna <- length(grna_odm[1, ]); n_c_resp <- length(response_odm[1, ])
cat("  response.odm:", length(rownames(response_odm)), "x", n_c_resp,
    "; colnames stored:", length(colnames(response_odm)), "\n")
cat("  grna.odm:    ", length(rownames(grna_odm)), "x", n_c_grna,
    "; colnames stored:", length(colnames(grna_odm)), "\n")
cat("  covariates:  ", nrow(cov_df), "rows; columns:", paste(names(cov_df), collapse = ", "), "\n")
cat("  first barcodes:", paste(head(rownames(cov_df), 3), collapse = ", "), "\n")

# ---- 1a. cell counts ---------------------------------------------------------
cat("\n1. cells\n")
check("1a both ODMs and the covariates have 207,324 cells",
      n_c_grna == 207324L && n_c_resp == 207324L && nrow(cov_df) == 207324L)

# ---- 2. guides and targets ---------------------------------------------------
cat("\n2. guides and targets\n")
cat("  target table columns:", paste(names(target_df), collapse = ", "), "\n")
check("2a grna.odm rows == target table guides (same order)", identical(rownames(grna_odm), target_df$grna_id))
target_csv <- read.csv(file.path(in_dir, "sceptre", "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
check("2b target table == sceptre/grna_target_data_frame.csv (grna_id, grna_target)",
      identical(target_csv$grna_id, target_df$grna_id) &&
        identical(target_csv$grna_target, target_df$grna_target))
ensg_guides  <- target_df$grna_id[startsWith(target_df$grna_target, "ENSG")]
ensg_targets <- unique(target_df$grna_target[startsWith(target_df$grna_target, "ENSG")])
cat("  ENSG targets:", length(ensg_targets), "; their guides:", length(ensg_guides), "\n")
check("2c 377 ENSG targets with 754 guides", length(ensg_targets) == 377L && length(ensg_guides) == 754L)

# ---- 3. target genes in the response matrix ---------------------------------
cat("\n3. target genes\n")
check("3a all 377 ENSG targets are rows of response.odm", all(ensg_targets %in% rownames(response_odm)))

# ---- grna.odm as a sparse matrix (row by row) --------------------------------
cat("\nreading grna.odm into a sparse matrix, row by row\n")
n_g <- length(rownames(grna_odm)); n_c <- n_c_grna
ilist <- jlist <- xlist <- vector("list", n_g)
for (g in seq_len(n_g)) {
  v  <- grna_odm[g, ]
  nz <- which(v != 0)
  ilist[[g]] <- rep.int(g, length(nz)); jlist[[g]] <- nz; xlist[[g]] <- v[nz]
}
odm_mat <- as_dgc(sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = as.numeric(unlist(xlist)),
                               dims = c(n_g, n_c)))
rm(ilist, jlist, xlist)
cat("  ", nrow(odm_mat), "x", ncol(odm_mat), ";", length(odm_mat@x), "nonzeros\n")

# ---- 1b/1c. cell order, via per-cell totals -----------------------------------
report_match <- function(label, a, b) {
  cat(sprintf("  %s: %d of %d cells equal\n", label, sum(a == b), length(a)))
  identical(as.numeric(a), as.numeric(b))
}
check("1b grna.odm column sums == grna_n_umis (cell order)",
      report_match("grna_n_umis", Matrix::colSums(odm_mat), cov_df$grna_n_umis))
check("1b grna.odm column nonzeros == grna_n_nonzero (cell order)",
      report_match("grna_n_nonzero", Matrix::colSums(odm_mat != 0), cov_df$grna_n_nonzero))

cat("\nreading response.odm row by row for per-cell totals\n")
resp_umis <- numeric(n_c_resp); resp_nonzero <- numeric(n_c_resp)
for (r in seq_along(rownames(response_odm))) {
  v <- response_odm[r, ]
  resp_umis <- resp_umis + v; resp_nonzero <- resp_nonzero + (v != 0)
}
check("1c response.odm column sums == response_n_umis (cell order)",
      report_match("response_n_umis", resp_umis, cov_df$response_n_umis))
check("1c response.odm column nonzeros == response_n_nonzero (cell order)",
      report_match("response_n_nonzero", resp_nonzero, cov_df$response_n_nonzero))
rm(resp_umis, resp_nonzero)

# ---- 4. fishash input (.mtx) --------------------------------------------------
cat("\n4. cleanser/grna_matrix.mtx\n")
mtx <- as_dgc(readMM(file.path(in_dir, "cleanser", "grna_matrix.mtx")))
cat("  ", nrow(mtx), "x", ncol(mtx), ";", length(mtx@x), "nonzeros\n")
check("4a .mtx == grna.odm (dims, order, every entry)", same_matrix(mtx, odm_mat))
rm(mtx)

# ---- 5. fishash+ input (.rds) -------------------------------------------------
cat("\n5. fishashplus/grna_matrix.rds\n")
rds <- readRDS(file.path(in_dir, "fishashplus", "grna_matrix.rds"))
cat("  class:", class(rds), ";", nrow(rds), "x", ncol(rds), "\n")
cat("  first rownames:", paste(head(rownames(rds), 3), collapse = ", "), "\n")
cat("  first colnames:", paste(head(colnames(rds), 3), collapse = ", "), "\n")
check("5a .rds rownames == grna.odm rows (same order)", identical(rownames(rds), rownames(grna_odm)))
check("5b .rds colnames are CELL_1..CELL_207324 or == covariate barcodes",
      identical(colnames(rds), paste0("CELL_", seq_len(n_c))) || identical(colnames(rds), rownames(cov_df)))
check("5c .rds == grna.odm (dims, order, every entry)", same_matrix(as_dgc(unname(rds)), odm_mat))
rm(rds)

# ---- 6. covariates -----------------------------------------------------------
cat("\n6. covariates\n")
full_cov <- read.csv(file.path(in_dir, "sceptre", "_full_cov_mat.csv"), row.names = 1, stringsAsFactors = FALSE)
num_cols <- c("response_n_nonzero", "response_n_umis", "response_p_mito", "grna_n_nonzero", "grna_n_umis")
check("6a _full_cov_mat.csv rows == covariate rows (same order)", identical(rownames(full_cov), rownames(cov_df)))
for (k in num_cols) {
  check(sprintf("6b %s equal", k), isTRUE(all.equal(as.numeric(full_cov[[k]]), as.numeric(cov_df[[k]]))))
}
check("6c prep_batch equal", identical(as.character(full_cov$prep_batch), as.character(cov_df$prep_batch)))

# ---- 7. fishash output indices -----------------------------------------------
cat("\n7. fishash output (gb_c1_full)\n")
fa <- read.csv(fishash_csv, colClasses = c("integer", "integer"))
cat("  ", nrow(fa), "assignments; cell_id range", paste(range(fa$cell_id), collapse = "-"),
    "; grna_id range", paste(range(fa$grna_id), collapse = "-"), "\n")
check("7a fishash indices within 1..13,077 guides and 1..207,324 cells",
      min(fa$grna_id) >= 1L && max(fa$grna_id) <= n_g && min(fa$cell_id) >= 1L && max(fa$cell_id) <= n_c)
check("7b fishash has no duplicate (cell, guide) rows", !anyDuplicated(fa))

# ---- summary -----------------------------------------------------------------
cat("\n==== summary ====\n")
for (k in names(results)) cat(sprintf("  %s  %s\n", results[k], k))
if (any(results == "FAIL")) quit(status = 1)
cat("all checks passed\n")
