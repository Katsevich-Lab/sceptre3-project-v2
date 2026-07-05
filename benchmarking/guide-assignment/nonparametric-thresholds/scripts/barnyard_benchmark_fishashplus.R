#!/usr/bin/env Rscript
# Authoritative barnyard head-to-head: fishash vs fishash+ (= depth_fix), on the authors'
# EXACT barnyard cohort + EXACT Table-2 accuracy metric (score_table2.R logic).
# fishash  = contingency_assign(cell_margin="observed")  [reproduces fishash]
# fishash+ = contingency_assign(cell_margin="ambient")   [our one change: denoised ambient depth]
# Everything else identical (min_count, hypergeometric tail, Guo-Sarkar FDR), so the ONLY
# difference is the cell-margin exposure -> isolates the depth fix.
suppressPackageStartupMessages({library(Matrix); library(fishash)})
source("scripts/contingency_method.R")
REPRO <- "external/repro_work"; OUT <- "results/collaborator_writeup"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

samples <- c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")
chem_of <- function(s) if (grepl("Cropseq", s)) "CROP-seq" else "direct-capture"

# authors' exact Table-2 accuracy: per QC cell, correct iff assigned guides are a single species
# matching the GEX species call; cells assigned both / neither species are NA (excluded).
table2_acc <- function(A, gt, orig_qc_pass, has_homo) {
  n_homo <- Matrix::colSums(A[gt == "homo_guide", , drop = FALSE] > 0)
  n_mus  <- Matrix::colSums(A[gt == "mus_guide",  , drop = FALSE] > 0)
  ah <- n_homo > 0; am <- n_mus > 0
  binz <- ah; binz[ah & am] <- NA; binz[!ah & !am] <- NA
  sel <- orig_qc_pass; pred <- binz[sel]; truth <- has_homo[sel]
  tp <- sum(pred & truth, na.rm=T); fp <- sum(pred & !truth, na.rm=T)
  tn <- sum(!pred & !truth, na.rm=T); fn <- sum(!pred & truth, na.rm=T)
  naT <- sum(is.na(pred) & truth); naF <- sum(is.na(pred) & !truth)
  (tp + tn) / (tp + fp + tn + fn + naT + naF)
}

out <- list()
for (s in samples) {
  meta   <- read.csv(file.path(REPRO, paste0(s, "_meta.csv")))
  guides <- read.csv(file.path(REPRO, paste0(s, "_guides.csv"))); gt <- guides$guide_type
  gm <- as(readMM(file.path(REPRO, paste0(s, "_grna_counts.mtx"))), "CsparseMatrix")
  rownames(gm) <- guides$guide; colnames(gm) <- paste0("c", seq_len(ncol(gm)))
  sum_gex <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_homo <- meta$homo_sum_gex / sum_gex; mito_frac <- meta$mito_sum / sum_gex
  orig_qc_pass <- (mito_frac < .15) & (meta$features_gex <= 6000) & (sum_gex <= 20000) &
                  (1500 <= meta$features_gex) & (3500 <= sum_gex) & (frac_homo < .1 | frac_homo > .9)
  orig_qc_pass[is.na(orig_qc_pass)] <- FALSE
  has_homo <- (frac_homo >= .1) & (sum_gex > 100)

  fa  <- contingency_assign(gm, q=0.05, refit=10, min_count=2, cell_margin="observed", tail="hyper", fdr="GS")
  fap <- contingency_assign(gm, q=0.05, refit=10, min_count=2, cell_margin="ambient",  tail="hyper", fdr="GS")
  out[[length(out)+1]] <- data.frame(sample=s, chemistry=chem_of(s),
    fishash  = round(table2_acc(fa$assigned,  gt, orig_qc_pass, has_homo), 4),
    `fishash+`= round(table2_acc(fap$assigned, gt, orig_qc_pass, has_homo), 4), check.names=FALSE)
  cat(sprintf("%-24s fishash=%.4f  fishash+=%.4f\n", s, out[[length(out)]]$fishash, out[[length(out)]]$`fishash+`))
}
df <- do.call(rbind, out)
df <- rbind(df, data.frame(sample="MEAN", chemistry="",
  fishash=round(mean(df$fishash),4), `fishash+`=round(mean(df$`fishash+`),4), check.names=FALSE))
write.csv(df, file.path(OUT, "barnyard_benchmark_fishashplus.csv"), row.names=FALSE)
cat("\n"); print(df, row.names=FALSE)
