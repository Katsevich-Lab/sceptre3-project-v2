#!/usr/bin/env Rscript
# Extract from the real Replogle RD7 genome-wide CRISPRi trans-analysis (~93M pairs,
# permissive QC): (1) all cis self-targeting pairs = real knockdowns (positive
# controls); (2) a trans sample for the null/entanglement analysis. Compute the
# margin filter m = n_trt * (n_nonzero_cntrl / n_NT).
suppressMessages({ library(arrow); library(dplyr) })
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
scr <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else normalizePath(".")
dir <- "/Users/ekatsevi/data/projects/sceptre3/nf_pipelines/rd7_trans/sceptre_outputs/trans_results"
nt  <- readRDS(file.path(scr, "rd7_ntrt.rds"))
n_NT <- nt$n_nt; ntrt_tab <- nt$n_trt

ds <- open_dataset(dir)
cat("total rows:", ds$num_rows, "\n")

## (1) all cis (self-targeting) pairs = positive controls
t0 <- Sys.time()
cis <- ds |> filter(response_id == grna_target) |>
  select(response_id, grna_target, n_nonzero_trt, n_nonzero_cntrl, pass_qc, p_value, log_2_fold_change) |>
  collect() |> as.data.frame()
cat(sprintf("cis pairs: %d  (%.0fs)\n", nrow(cis), as.numeric(Sys.time()-t0, units="secs")))

## (2) trans sample for null/entanglement: read first 8 files fully
files <- list.files(dir, pattern="result_[0-9]+\\.parquet$", full.names=TRUE)
set.seed(1); some <- files[1:8]
trans <- open_dataset(some) |> filter(response_id != grna_target) |>
  select(response_id, grna_target, n_nonzero_trt, n_nonzero_cntrl, p_value, log_2_fold_change) |>
  collect() |> as.data.frame()
cat("trans sample rows:", nrow(trans), "\n")
# subsample trans to 800k for speed
if (nrow(trans) > 800000) trans <- trans[sample(nrow(trans), 800000), ]

addm <- function(d){
  d$n_trt <- ntrt_tab[d$grna_target]
  d$p_hat_cntrl <- d$n_nonzero_cntrl / n_NT
  d$m_expected  <- d$n_trt * d$p_hat_cntrl
  d$pass_old <- d$n_nonzero_trt >= 7 & d$n_nonzero_cntrl >= 7
  d$pass_new <- d$m_expected    >= 7 & d$n_nonzero_cntrl >= 7
  d
}
cis <- addm(cis); trans <- addm(trans)
saveRDS(list(cis=cis, trans=trans, n_NT=n_NT), file.path(scr, "rd7_extract.rds"))

## quick summary
cat("\n==== Replogle RD7 CIS positive controls (real knockdowns) ====\n")
cat("n cis pairs with n_trt known:", sum(!is.na(cis$n_trt)), "\n")
cis <- cis[!is.na(cis$n_trt) & is.finite(cis$p_value), ]
cat(sprintf("median log2FC = %.2f  (fraction with log2FC< -0.5: %.2f)\n",
            median(cis$log_2_fold_change), mean(cis$log_2_fold_change < -0.5)))
# real strong knockdowns: significant & downregulated
cis$strong <- cis$p_value < 1e-4 & cis$log_2_fold_change < -0.5
cat(sprintf("strong real knockdowns (p<1e-4 & log2FC<-0.5): %d\n", sum(cis$strong)))
sk <- cis[cis$strong, ]
cat(sprintf("  of these, n_nonzero_trt < 7 (DROPPED by current filter): %d (%.1f%%)\n",
            sum(sk$n_nonzero_trt < 7), 100*mean(sk$n_nonzero_trt < 7)))
cat(sprintf("  of these, KEPT by margin filter (m>=7): %d (%.1f%%)\n",
            sum(sk$pass_new), 100*mean(sk$pass_new)))
dropped <- sk[sk$n_nonzero_trt < 7, ]
cat(sprintf("  among dropped-by-current strong knockdowns: rescued by margin = %d (%.1f%%)\n",
            sum(dropped$pass_new), 100*mean(dropped$pass_new)))
cat(sprintf("  median n_nonzero_trt of dropped strong knockdowns: %.0f ; median m: %.0f\n",
            median(dropped$n_nonzero_trt), median(dropped$m_expected)))

cat("\n-- overall cis discovery counts (BH<0.1, left-tail knockdown) --\n")
cis$p_kd <- ifelse(cis$log_2_fold_change < 0, cis$p_value/2, 1 - cis$p_value/2) # approx one-sided
disc <- function(sub){ sub<-sub[is.finite(sub$p_value),]; sum(p.adjust(sub$p_value,"BH")<0.1)}
cat(sprintf("current filter keeps %d cis pairs -> %d discoveries\n", sum(cis$pass_old), disc(cis[cis$pass_old,])))
cat(sprintf("margin  filter keeps %d cis pairs -> %d discoveries\n", sum(cis$pass_new), disc(cis[cis$pass_new,])))

cat("\n==== TRANS (null-ish) entanglement ====\n")
tr <- trans[is.finite(trans$p_value) & !is.na(trans$n_trt), ]
tr$z <- sign(tr$log_2_fold_change) * qnorm(pmin(pmax(1 - tr$p_value/2, 1e-8), 1-1e-8))
cat(sprintf("cor(n_nonzero_trt, z) = %+.3f | cor(m, z) = %+.3f\n",
            cor(tr$n_nonzero_trt, tr$z, use="complete.obs"), cor(tr$m_expected, tr$z, use="complete.obs")))
cat(sprintf("mean z among n_nonzero_trt<7: %+.3f | >=7: %+.3f\n",
            mean(tr$z[tr$n_nonzero_trt<7]), mean(tr$z[tr$n_nonzero_trt>=7])))
cat("DONE -> rd7_extract.rds\n")
