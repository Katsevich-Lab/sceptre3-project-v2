#!/usr/bin/env Rscript
# Build the full-Gasperini positive-control dataset for comparing guide-assignment
# methods: one shared response matrix + covariates, and one binary gRNA matrix
# per assignment method. Cells: all 207,324 (no assignment-based selection).
# Targets: the 377 TSS (ENSG) targets, 754 guides; no NT guides (high MOI uses
# the complement control group, so NT rows would not change any p-value).
#
#   SIF=<pipeline-images>/sceptre.sif \
#     sbatch -c 2 -t 1:00:00 ../../guide-assignment/grna-simulator/run-in-image.sh build-dataset.R
#
# Inputs (checked by check-inputs.R, job 9173649, all PASS):
#   guide_assignment/input_data/gasperini/sceptre-pipeline/{sceptre_object.rds, response.odm, grna.odm}
#     ODM column j == covariate row j (the ODMs store no cell names; barcodes
#     are the covariate rownames)
#   guide_assignment/outputs/gb_c1_full/assignments_fishash_gasperini.csv
#     1-based (cell, guide) indices into cleanser/grna_matrix.mtx == grna.odm order
#   guide_assignment/outputs/gasperini_fishashplus_refit100/assignments_fishashplus_gasperini.csv
#     cell CELL_j == column j; guide by ID; converged at REFIT = 100
#
# Writes <LOCAL_BENCHMARKING_DIR>/association/assign-pos-control/input_data/gasperini_allcells/:
#   shared/response_matrix.rds        377 genes x 207,324 cells, raw counts
#   shared/cell_covariates.csv        cell_id + the 4 *_full covariates + prep_batch
#   shared/grna_target_data_frame.csv 754 guides: grna_id, grna_target, chr, start, end
#   shared/qc_excluded_cells.rds      cells sceptre's default cell-QC rule removes,
#                                     applied to the FULL-transcriptome covariates
#   grna/<method>.rds                 binary 754 guides x 207,324 cells, one per method
#   assigned_cells_per_target.csv     cells assigned per target and method (all / after QC)
# Refuses to run if the output directory already exists.

suppressPackageStartupMessages({
  library(Matrix); library(sceptre); library(ondisc); library(data.table)
})
source("~/.Rprofile")                                    # .get_config_path()

bench    <- .get_config_path("LOCAL_BENCHMARKING_DIR")
in_dir   <- file.path(bench, "guide_assignment", "input_data", "gasperini")
pipe_dir <- file.path(in_dir, "sceptre-pipeline")
fishash_csv     <- file.path(bench, "guide_assignment", "outputs", "gb_c1_full",
                             "assignments_fishash_gasperini.csv")
fishashplus_csv <- file.path(bench, "guide_assignment", "outputs", "gasperini_fishashplus_refit100",
                             "assignments_fishashplus_gasperini.csv")
out_root <- file.path(bench, "association", "assign-pos-control", "input_data", "gasperini_allcells")

THRESHOLDS <- c(3L, 5L, 10L)          # thresholding: assigned iff raw count >= t
N_CELLS    <- 207324L
# sceptre 0.99.0 run_qc() defaults, applied here to the full-transcriptome covariates
QC_UMIS_RANGE    <- c(0.01, 0.99)
QC_NONZERO_RANGE <- c(0.01, 0.99)
QC_P_MITO        <- 0.2
# expected totals in the two assignment files (checked before subsetting)
N_FISHASH     <- 6971335L
N_FISHASHPLUS <- 7405793L

if (dir.exists(out_root)) stop("output directory already exists, not overwriting: ", out_root)

# ---- load --------------------------------------------------------------------
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
barcodes  <- rownames(cov_df)
stopifnot(length(barcodes) == N_CELLS, !anyDuplicated(barcodes),
          identical(rownames(grna_odm), target_df$grna_id))

# ---- targets and guides ------------------------------------------------------
tss_df <- target_df[startsWith(target_df$grna_target, "ENSG"),
                    c("grna_id", "grna_target", "chr", "start", "end")]
rownames(tss_df) <- NULL
genes <- unique(tss_df$grna_target)
stopifnot(nrow(tss_df) == 754L, length(genes) == 377L,
          all(genes %in% rownames(response_odm)))
tss_rows <- match(tss_df$grna_id, rownames(grna_odm))    # rows of grna.odm, in tss_df order
cat("targets:", length(genes), "; guides:", nrow(tss_df), "\n")

# ---- response matrix (377 genes x all cells) ---------------------------------
cat("\nreading response.odm rows for the 377 target genes\n")
read_rows <- function(odm, rows, n_cells) {
  ilist <- jlist <- xlist <- vector("list", length(rows))
  for (k in seq_along(rows)) {
    v  <- odm[rows[k], ]
    stopifnot(length(v) == n_cells)
    nz <- which(v != 0)
    ilist[[k]] <- rep.int(k, length(nz)); jlist[[k]] <- nz; xlist[[k]] <- v[nz]
  }
  sparseMatrix(i = unlist(ilist), j = unlist(jlist), x = as.numeric(unlist(xlist)),
               dims = c(length(rows), n_cells))
}
response_matrix <- read_rows(response_odm, genes, N_CELLS)
dimnames(response_matrix) <- list(genes, barcodes)
# cross-check each gene against gene_summary_stats.csv (computed from the same data)
gss <- read.csv(file.path(in_dir, "gene_summary_stats.csv"), stringsAsFactors = FALSE)
gss <- gss[match(genes, gss$gene), ]
stopifnot(identical(gss$gene, genes),
          all(Matrix::rowSums(response_matrix) == gss$gene_total),
          all(Matrix::rowSums(response_matrix != 0) == gss$gene_n_nonzero))
cat("  ", nrow(response_matrix), "x", ncol(response_matrix), ";", length(response_matrix@x),
    "nonzeros; per-gene totals and nonzero counts match gene_summary_stats.csv\n")

# ---- covariates (full data, assignment-independent) --------------------------
cell_covariates <- data.frame(
  cell_id                 = barcodes,
  response_n_nonzero_full = cov_df$response_n_nonzero,
  response_n_umis_full    = cov_df$response_n_umis,
  grna_n_nonzero_full     = cov_df$grna_n_nonzero,
  grna_n_umis_full        = cov_df$grna_n_umis,
  prep_batch              = as.character(cov_df$prep_batch),
  stringsAsFactors = FALSE
)

# ---- cell QC on the full-transcriptome covariates ----------------------------
# Same rule as sceptre 0.99.0 determine_cells_to_retain(): drop cells outside the
# [lo, hi] quantiles (stats::quantile, strict < and >) of response_n_umis and of
# response_n_nonzero, and cells with response_p_mito > threshold.
clip <- function(v, probs) {
  cutoffs <- stats::quantile(x = v, probs = probs)
  which(v < cutoffs[1] | v > cutoffs[2])
}
excl_umis    <- clip(cov_df$response_n_umis, QC_UMIS_RANGE)
excl_nonzero <- clip(cov_df$response_n_nonzero, QC_NONZERO_RANGE)
excl_p_mito  <- which(cov_df$response_p_mito > QC_P_MITO)
qc_excluded  <- sort(unique(c(excl_umis, excl_nonzero, excl_p_mito)))
qc_kept      <- setdiff(seq_len(N_CELLS), qc_excluded)
cat("\ncell QC (full covariates): n_umis removes", length(excl_umis),
    "; n_nonzero removes", length(excl_nonzero), "; p_mito removes", length(excl_p_mito),
    "; total removed", length(qc_excluded), "; kept", length(qc_kept), "\n")

# ---- raw gRNA counts for the 754 guides --------------------------------------
cat("\nreading grna.odm rows for the 754 guides\n")
grna_counts <- read_rows(grna_odm, tss_rows, N_CELLS)    # rows in tss_df order
cat("  ", nrow(grna_counts), "x", ncol(grna_counts), ";", length(grna_counts@x), "nonzeros\n")

binary <- function(i, j) {
  m <- sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(tss_df), N_CELLS))
  dimnames(m) <- list(tss_df$grna_id, barcodes)
  m
}
grna <- list()

# ---- fishash: 1-based (cell, guide) indices in grna.odm order ----------------
cat("\nfishash:", fishash_csv, "\n")
fa <- fread(fishash_csv, colClasses = c("integer", "integer"))
stopifnot(identical(names(fa), c("cell_id", "grna_id")), nrow(fa) == N_FISHASH,
          !anyDuplicated(fa), all(fa$cell_id %in% seq_len(N_CELLS)),
          all(fa$grna_id %in% seq_along(rownames(grna_odm))))
fa <- fa[grna_id %in% tss_rows]
grna$fishash <- binary(match(fa$grna_id, tss_rows), fa$cell_id)
fa_counts <- grna_counts[cbind(match(fa$grna_id, tss_rows), fa$cell_id)]
cat("  ", N_FISHASH, "assignments in total;", nrow(fa), "on the 754 TSS guides;",
    "raw count at those entries: min", min(fa_counts), ", zero in", sum(fa_counts == 0), "\n")
rm(fa)

# ---- fishash+: cell CELL_j == column j; guide by ID --------------------------
cat("\nfishash+:", fishashplus_csv, "\n")
fp <- fread(fishashplus_csv, colClasses = c("character", "character"))
stopifnot(identical(names(fp), c("cell_id", "grna_id")), nrow(fp) == N_FISHASHPLUS,
          !anyDuplicated(fp), all(grepl("^CELL_[0-9]+$", fp$cell_id)),
          all(fp$grna_id %in% rownames(grna_odm)))
fp[, j := as.integer(sub("^CELL_", "", cell_id))]
stopifnot(all(fp$j %in% seq_len(N_CELLS)))
fp <- fp[grna_id %in% tss_df$grna_id]
grna$fishashplus <- binary(match(fp$grna_id, tss_df$grna_id), fp$j)
fp_counts <- grna_counts[cbind(match(fp$grna_id, tss_df$grna_id), fp$j)]
cat("  ", N_FISHASHPLUS, "assignments in total;", nrow(fp), "on the 754 TSS guides;",
    "raw count at those entries: min", min(fp_counts), "\n")
# run_fishashplus.R sets min_count = 2, so every assigned entry must have count >= 2;
# a wrong cell or guide mapping would put assignments on entries below that
stopifnot(min(fp_counts) >= 2)
rm(fp)

# ---- thresholding: raw count >= t --------------------------------------------
nz <- summary(grna_counts)                               # i, j, x of nonzero entries
for (t in THRESHOLDS) {
  keep <- nz$x >= t
  grna[[paste0("thresh", t)]] <- binary(nz$i[keep], nz$j[keep])
}

# ---- per-target assigned cells ----------------------------------------------
per_target <- rbindlist(lapply(names(grna), function(m) {
  by_target <- lapply(genes, function(g) {
    rows <- which(tss_df$grna_target == g)
    any_guide <- Matrix::colSums(grna[[m]][rows, , drop = FALSE]) > 0
    data.table(method = m, grna_target = g,
               n_cells_all = sum(any_guide), n_cells_qc = sum(any_guide[qc_kept]))
  })
  rbindlist(by_target)
}))
cat("\nassigned cells per target (union of its 2 guides):\n")
print(per_target[, .(total_assignments_754 = sum(grna[[method[1]]]),
                     min_cells = min(n_cells_all), median_cells = median(n_cells_all),
                     targets_with_0 = sum(n_cells_all == 0),
                     targets_with_0_after_qc = sum(n_cells_qc == 0)), by = method])

# ---- write -------------------------------------------------------------------
for (d in c("shared", "grna")) dir.create(file.path(out_root, d), recursive = TRUE)
saveRDS(response_matrix, file.path(out_root, "shared", "response_matrix.rds"))
write.csv(cell_covariates, file.path(out_root, "shared", "cell_covariates.csv"), row.names = FALSE)
write.csv(tss_df, file.path(out_root, "shared", "grna_target_data_frame.csv"), row.names = FALSE)
saveRDS(list(excluded_idx = qc_excluded, excluded_barcodes = barcodes[qc_excluded],
             n_removed = c(n_umis = length(excl_umis), n_nonzero = length(excl_nonzero),
                           p_mito = length(excl_p_mito), total = length(qc_excluded)),
             rule = list(response_n_umis_range = QC_UMIS_RANGE,
                         response_n_nonzero_range = QC_NONZERO_RANGE,
                         p_mito_threshold = QC_P_MITO)),
        file.path(out_root, "shared", "qc_excluded_cells.rds"))
for (m in names(grna)) saveRDS(grna[[m]], file.path(out_root, "grna", paste0(m, ".rds")))
fwrite(per_target, file.path(out_root, "assigned_cells_per_target.csv"))
cat("\nwrote", out_root, "\n")
