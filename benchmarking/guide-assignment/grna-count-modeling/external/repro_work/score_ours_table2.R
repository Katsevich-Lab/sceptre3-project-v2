#!/usr/bin/env Rscript
# Score OUR methods (ambient test, Otsu, valley) on the authors' EXACT barnyard
# cohort + EXACT Table-2 metric, so they are directly comparable to the
# now-exactly-reproduced geomux/fishash numbers.
suppressPackageStartupMessages({library(Matrix)})
WORK <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
LIB  <- normalizePath(file.path(WORK, "..", "..", "..", "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
source(LIB)

binary_accuracy_metrics <- function(pred, truth) {
  tp <- sum(pred & truth, na.rm=T); fp <- sum(pred & !truth, na.rm=T)
  tn <- sum(!pred & !truth, na.rm=T); fn <- sum(!pred & truth, na.rm=T)
  naTrue <- sum(is.na(pred) & truth); naFalse <- sum(is.na(pred) & !truth)
  tot <- tp+fp+tn+fn+naTrue+naFalse
  c(Accuracy=(tp+tn)/tot, stderr=sqrt((tp+fn+naTrue)*(tp/(tp+fn+naTrue))*(1-tp/(tp+fn+naTrue))/tot^2))
}
samples <- c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")
methods <- list(
  `ours (ambient)` = function(gm) ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=1)$assignment_matrix,
  otsu   = function(gm) assign_by_threshold(gm, otsu_threshold_log1p)$assignment_matrix,
  valley = function(gm) assign_by_threshold(gm, smoothed_valley_threshold)$assignment_matrix)

out <- list()
for (s in samples) {
  meta   <- read.csv(file.path(WORK, paste0(s, "_meta.csv")))
  guides <- read.csv(file.path(WORK, paste0(s, "_guides.csv"))); gt <- guides$guide_type
  gm <- as(readMM(file.path(WORK, paste0(s, "_grna_counts.mtx"))), "CsparseMatrix")  # guides x cells
  rownames(gm) <- guides$guide
  sum_gex <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_homo <- meta$homo_sum_gex / sum_gex; mito_frac <- meta$mito_sum / sum_gex
  orig_qc_pass <- (mito_frac < .15) & (meta$features_gex <= 6000) & (sum_gex <= 20000) &
                  (1500 <= meta$features_gex) & (3500 <= sum_gex) & (frac_homo < .1 | frac_homo > .9)
  has_homo <- (frac_homo >= .1) & (sum_gex > 100)
  for (mn in names(methods)) {
    A <- methods[[mn]](gm)
    n_homo <- Matrix::colSums(A[gt=="homo_guide", , drop=FALSE] > 0)
    n_mus  <- Matrix::colSums(A[gt=="mus_guide",  , drop=FALSE] > 0)
    ah <- n_homo > 0; am <- n_mus > 0
    binz <- ah; binz[ah & am] <- NA; binz[!ah & !am] <- NA
    sel <- orig_qc_pass
    r <- binary_accuracy_metrics(binz[sel], has_homo[sel])
    out[[length(out)+1]] <- data.frame(sample=s, method=mn, Accuracy=round(r["Accuracy"],4))
  }
}
df <- do.call(rbind, out); rownames(df) <- NULL
w <- reshape(df, idvar="method", timevar="sample", direction="wide"); names(w) <- sub("Accuracy.","",names(w))
print(w, row.names=FALSE)
write.csv(w, file.path(WORK, "ours_table2_exact.csv"), row.names=FALSE)
