#!/usr/bin/env Rscript
# Controlled simulation: sceptre's current pairwise QC (observed n_nonzero_trt>=7)
# vs a margin-based independent filter (expected nonzero trt m = n_trt*p_ctrl>=7).
# Ground truth known per pair. Real sceptre is run with permissive QC so every
# pair is tested; the two filters are applied post hoc.

suppressMessages({ library(sceptre); library(Matrix) })
set.seed(1)
args   <- commandArgs(trailingOnly = TRUE)
SCALE  <- ifelse(length(args) >= 1, as.numeric(args[1]), 1)
SIDE   <- ifelse(length(args) >= 2, args[2], "both")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
outdir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else normalizePath(".")

theta       <- 5
mu_of_p     <- function(p, theta) theta * ((1 - p)^(-1 / theta) - 1)

## ---- design: power grid (true effects + kappa=0 controls) ----
p_levels     <- c(0.15, 0.30, 0.50)
ntrt_levels  <- c(40, 80)
kappa_levels <- c(0, 0.50, 0.70, 0.85, 0.95)
R_grid       <- round(14 * SCALE)
grid <- expand.grid(p = p_levels, n_trt = ntrt_levels, kappa = kappa_levels,
                    rep = seq_len(R_grid), KEEP.OUT.ATTRS = FALSE)

## ---- design: large null panel for calibration (kappa=0, spread over spectrum) ----
n_null <- round(500 * SCALE)
nullp  <- data.frame(
  p     = sample(p_levels, n_null, replace = TRUE),
  n_trt = sample(ntrt_levels, n_null, replace = TRUE),
  kappa = 0, rep = NA)
grid <- rbind(grid, nullp)

grid$panel   <- c(rep("grid", nrow(grid) - n_null), rep("null", n_null))
grid$gene_id <- sprintf("gene_%04d", seq_len(nrow(grid)))
grid$grna_id <- sprintf("g_%04d",   seq_len(nrow(grid)))
grid$mu0     <- mu_of_p(grid$p, theta)
G_test <- nrow(grid)

n_NT_cells <- round(10000 * SCALE)
G_bg       <- 150
N <- n_NT_cells + sum(grid$n_trt)
cat(sprintf("Cells: %d | test guides: %d (grid %d + null %d) | NT cells: %d\n",
            N, G_test, sum(grid$panel=="grid"), n_null, n_NT_cells))

s <- rlnorm(N, 0, 0.3)   # library size factor

## ---- assign cells to guides (low-MOI) ----
cell_guide <- character(N)
n_NT_grnas <- 20
cell_guide[seq_len(n_NT_cells)] <- sprintf("nt_%02d", sample.int(n_NT_grnas, n_NT_cells, replace = TRUE))
ptr <- n_NT_cells
trt_cell_list <- vector("list", G_test)
for (i in seq_len(G_test)) {
  idx <- (ptr + 1):(ptr + grid$n_trt[i]); cell_guide[idx] <- grid$grna_id[i]
  trt_cell_list[[i]] <- idx; ptr <- ptr + grid$n_trt[i]
}

## ---- gRNA matrix ----
all_grnas <- c(sprintf("nt_%02d", seq_len(n_NT_grnas)), grid$grna_id)
grna_mat  <- matrix(0L, length(all_grnas), N, dimnames = list(all_grnas, NULL))
carrier_row <- match(cell_guide, all_grnas)
grna_mat[cbind(carrier_row, seq_len(N))] <- rpois(N, 50) + 20L
grna_mat <- grna_mat + matrix(rpois(length(all_grnas) * N, 0.02), length(all_grnas))
grna_tdf <- data.frame(grna_id = all_grnas,
                       grna_target = c(rep("non-targeting", n_NT_grnas), grid$gene_id),
                       stringsAsFactors = FALSE)

## ---- response matrix ----
all_genes <- c(grid$gene_id, sprintf("bg_%03d", seq_len(G_bg)))
resp <- matrix(0L, length(all_genes), N, dimnames = list(all_genes, NULL))
for (i in seq_len(G_test)) {
  eff <- rep(1, N); eff[trt_cell_list[[i]]] <- (1 - grid$kappa[i])
  resp[i, ] <- rnbinom(N, mu = grid$mu0[i] * s * eff, size = theta)
}
bg_mu0 <- runif(G_bg, 0.2, 3)
for (b in seq_len(G_bg)) resp[G_test + b, ] <- rnbinom(N, mu = bg_mu0[b] * s, size = theta)

## ---- run sceptre, permissive QC ----
disc_pairs <- data.frame(grna_target = grid$gene_id, response_id = grid$gene_id,
                         stringsAsFactors = FALSE)
so <- import_data(resp, grna_mat, grna_tdf, moi = "low") |>
  set_analysis_parameters(discovery_pairs = disc_pairs, side = SIDE,
                          resampling_mechanism = "permutations") |>
  assign_grnas(method = "thresholding", threshold = 5) |>
  run_qc(n_nonzero_trt_thresh = 0L, n_nonzero_cntrl_thresh = 0L, p_mito_threshold = 1)
so <- run_discovery_analysis(so)
res <- as.data.frame(get_result(so, "run_discovery_analysis"))

## ---- join ground truth; compute margin filter ----
res <- merge(res, grid[, c("gene_id","p","n_trt","kappa","mu0","panel")],
             by.x = "response_id", by.y = "gene_id", all.x = TRUE)
n_cntrl_cells    <- sum(grepl("^nt_", cell_guide))
res$p_hat_cntrl  <- res$n_nonzero_cntrl / n_cntrl_cells
res$m_expected   <- res$n_trt * res$p_hat_cntrl            # margin-based filter statistic
res$pass_old     <- res$n_nonzero_trt >= 7 & res$n_nonzero_cntrl >= 7
res$pass_new     <- res$m_expected    >= 7 & res$n_nonzero_cntrl >= 7
saveRDS(res, file.path(outdir, sprintf("sim_qc_res_%s.rds", SIDE)))

## ================= SUMMARY =================
cat("\n===== SUMMARY (side =", SIDE, ") =====\nTested pairs:", nrow(res), "\n")

cat("\n--- Calibration on NULL pairs (kappa=0) ---\n")
nl <- res[res$kappa == 0 & !is.na(res$p_value), ]
cal <- function(x) sprintf("n=%d  P(p<.05)=%.3f  P(p<.01)=%.3f", length(x),
                           mean(x < .05), mean(x < .01))
cat(" all nulls        :", cal(nl$p_value), "\n")
cat(" pass_old nulls   :", cal(nl$p_value[nl$pass_old]), "\n")
cat(" pass_new nulls   :", cal(nl$p_value[nl$pass_new]), "\n")
# nulls the NEW filter keeps but OLD drops (the newly-admitted region): still uniform?
newonly <- nl[nl$pass_new & !nl$pass_old, ]
cat(" new-only nulls   :", cal(newonly$p_value), "  <-- must stay calibrated\n")

cat("\n--- corr(filter statistic, |z|) under null (want ~0) ---\n")
z_null <- abs(qnorm(pmax(nl$p_value, 1e-6) / (if (SIDE=="both") 2 else 1)))
cat(sprintf(" cor(n_nonzero_trt, |z|) = %+.3f   <-- current filter\n",
            cor(nl$n_nonzero_trt, z_null, use = "complete.obs")))
cat(sprintf(" cor(m_expected,    |z|) = %+.3f   <-- margin filter\n",
            cor(nl$m_expected,    z_null, use = "complete.obs")))

cat("\n--- Retention of TRUE effects by filter (mean pass rate) ---\n")
te <- res[res$panel == "grid" & res$kappa > 0, ]
agg <- aggregate(cbind(pass_old, pass_new) ~ kappa + p + n_trt, te, mean)
print(agg[order(agg$p, agg$n_trt, agg$kappa), ], row.names = FALSE)

cat("\n--- Real detectable signal discarded by OLD filter ---\n")
te$detectable <- te$p_value < 0.05
cat(sprintf(" true-effect pairs detectable (p<.05): %d / %d\n", sum(te$detectable), nrow(te)))
lost <- te[te$detectable & !te$pass_old, ]
cat(sprintf("   dropped by OLD filter : %d (%.1f%% of detectable)\n",
            nrow(lost), 100 * nrow(lost) / sum(te$detectable)))
cat(sprintf("   of those kept by NEW  : %d (%.1f%%)\n",
            sum(lost$pass_new), 100 * mean(lost$pass_new)))

cat("\n--- DESeq2-style bottom line: BH discoveries within each filtered set ---\n")
disc_count <- function(sub) {
  sub <- sub[!is.na(sub$p_value), ]
  padj <- p.adjust(sub$p_value, "BH")
  data.frame(tested = nrow(sub),
             true_disc  = sum(padj < 0.1 & sub$kappa > 0),
             false_disc = sum(padj < 0.1 & sub$kappa == 0))
}
cat(" OLD filter:\n"); print(disc_count(res[res$pass_old, ]), row.names = FALSE)
cat(" NEW filter:\n"); print(disc_count(res[res$pass_new, ]), row.names = FALSE)
cat("\nDONE ->", sprintf("sim_qc_res_%s.rds", SIDE), "\n")
