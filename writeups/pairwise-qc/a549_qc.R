#!/usr/bin/env Rscript
# Real-data test on a549 CRISPRi (GSE236304): every targeting guide is a TSS
# knockdown of its own gene (real positive control). We subsample each target's
# treatment cells to a controlled n_trt (mimicking realistic per-perturbation
# cell counts), keep the full data for model fitting, and ask: how many REAL
# knockdowns does the current filter (observed n_nonzero_trt>=7) discard that the
# margin filter (expected m = n_trt * p_ctrl >=7) keeps? Plus null calibration.

suppressMessages({ library(sceptre); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
C    <- ifelse(length(args) >= 1, as.integer(args[1]), 80L)   # treatment cells per positive control
dir  <- "/Users/ekatsevi/data/external/perturbseq-survey/a549_crispri_perturbseq_GSE236304/sceptre"
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
out <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")

g   <- readRDS(file.path(dir, "grna_matrix_aligned.rds"))
r   <- readRDS(file.path(dir, "response_matrix.rds"))
tdf <- read.csv(file.path(dir, "grna_target_data_frame.csv"), stringsAsFactors = FALSE)
N   <- ncol(r)
tgts <- unique(tdf$grna_target[tdf$grna_target != "non-targeting"])
assigned <- g >= 5
p_all <- Matrix::rowSums(r > 0) / N
SEEDS <- 1:5

## per-target treatment cells (union over its guides)
target_cells <- lapply(tgts, function(t) {
  gids <- tdf$grna_id[tdf$grna_target == t]; gids <- gids[gids %in% rownames(g)]
  which(Matrix::colSums(assigned[gids, , drop = FALSE]) > 0)
}); names(target_cells) <- tgts

## build pseudo-guides: POS (subsampled real knockdown) + NULL (random cells)
pg_cells <- list(); meta <- list()
for (t in tgts) {
  cells <- target_cells[[t]]; if (length(cells) < C) next
  for (s in SEEDS) {
    set.seed(1000 * s + which(tgts == t)); sub <- sample(cells, C)
    id <- sprintf("POS_%s_s%d", t, s)
    pg_cells[[id]] <- sub
    meta[[id]] <- data.frame(pg = id, type = "pos", gene = t, seed = s, p_base = p_all[t])
  }
}
cand <- rownames(r)[p_all > 0.03 & p_all < 0.97]
n_null <- 1600
for (i in seq_len(n_null)) {
  set.seed(500000 + i); sub <- sample(N, C); gene <- sample(cand, 1)
  id <- sprintf("NULL_%06d", i)
  pg_cells[[id]] <- sub
  meta[[id]] <- data.frame(pg = id, type = "null", gene = gene, seed = NA, p_base = p_all[gene])
}
meta <- do.call(rbind, meta); rownames(meta) <- NULL
pg_ids <- names(pg_cells)
cat(sprintf("C=%d | pos pseudo-guides=%d (%d targets) | null=%d\n",
            C, sum(meta$type=="pos"), length(unique(meta$gene[meta$type=="pos"])), sum(meta$type=="null")))

## pseudo-guide count matrix
ii <- unlist(lapply(seq_along(pg_ids), function(k) rep(k, length(pg_cells[[k]]))))
jj <- unlist(pg_cells)
gm <- sparseMatrix(i = ii, j = jj, x = 20L, dims = c(length(pg_ids), N),
                   dimnames = list(pg_ids, NULL))
grna_tdf2 <- data.frame(grna_id = pg_ids, grna_target = pg_ids, stringsAsFactors = FALSE)
disc <- data.frame(grna_target = meta$pg, response_id = meta$gene, stringsAsFactors = FALSE)

## run sceptre (high-MOI, complement control, left-sided for knockdown)
so <- import_data(r, gm, grna_tdf2, moi = "high") |>
  set_analysis_parameters(
    discovery_pairs = disc, side = "left",
    formula_object = formula(~ log(response_n_umis + 1) + log(response_n_nonzero + 1)),
    resampling_mechanism = "permutations") |>
  assign_grnas(method = "thresholding", threshold = 5) |>
  run_qc(n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L,
         response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1), p_mito_threshold = 1)
so <- run_discovery_analysis(so)
res <- as.data.frame(get_result(so, "run_discovery_analysis"))

## join meta, compute margin filter
res <- merge(res, meta, by.x = c("grna_target", "response_id"), by.y = c("pg", "gene"))
n_cntrl <- N - C
res$p_hat_cntrl <- res$n_nonzero_cntrl / (N - C)     # complement control fraction
res$m_expected  <- C * res$p_hat_cntrl
res$pass_old <- res$n_nonzero_trt >= 7 & res$n_nonzero_cntrl >= 7
res$pass_new <- res$m_expected    >= 7 & res$n_nonzero_cntrl >= 7
saveRDS(res, file.path(out, sprintf("a549_res_C%d.rds", C)))

## ===== summary =====
cat("\n==== a549 real CRISPRi, C =", C, "====\n")
nl <- res[res$type == "null" & is.finite(res$p_value), ]
po <- res[res$type == "pos"  & is.finite(res$p_value), ]
cat(sprintf("NULL calibration: n=%d  P(p<.05)=%.3f  (pass_old=%.3f, pass_new=%.3f)\n",
            nrow(nl), mean(nl$p_value < .05),
            mean(nl$p_value[nl$pass_old] < .05), mean(nl$p_value[nl$pass_new] < .05)))
zsl <- qnorm(pmin(pmax(nl$p_value, 1e-6), 1 - 1e-6))   # left-sided signed z
cat(sprintf("Entanglement (null): cor(n_nonzero_trt, z)=%+.3f  cor(m, z)=%+.3f\n",
            cor(nl$n_nonzero_trt, zsl), cor(nl$m_expected, zsl)))

cat("\nPOSITIVE CONTROLS (real knockdowns):\n")
po$sig <- po$p_value < 0.05
cat(sprintf("  significant (p<.05, left): %d / %d ; median log2FC = %.2f\n",
            sum(po$sig), nrow(po), median(po$log_2_fold_change, na.rm = TRUE)))
lost <- po[po$sig & !po$pass_old, ]
cat(sprintf("  real detectable knockdowns DROPPED by current filter: %d (%.1f%% of detectable)\n",
            nrow(lost), 100 * nrow(lost) / max(1, sum(po$sig))))
cat(sprintf("     of those RESCUED by margin filter: %d (%.1f%%)\n",
            sum(lost$pass_new), 100 * mean(lost$pass_new)))
cat(sprintf("  median log2FC of dropped-but-real: %.2f (strong knockdowns)\n",
            median(lost$log_2_fold_change, na.rm = TRUE)))

## BH discoveries among positive controls under each filter
bh <- function(sub){ sub<-sub[is.finite(sub$p_value),]; sum(p.adjust(sub$p_value,"BH")<0.1) }
cat(sprintf("\n  positive-control discoveries (BH<0.1): current=%d  margin=%d\n",
            bh(po[po$pass_old,]), bh(po[po$pass_new,])))
cat("DONE ->", sprintf("a549_res_C%d.rds", C), "\n")
