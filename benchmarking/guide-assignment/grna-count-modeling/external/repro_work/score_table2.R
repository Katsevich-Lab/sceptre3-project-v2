#!/usr/bin/env Rscript
# Faithful replication of fishash_analysis 82_plot_barnyard_results.R Table-2 logic.
suppressPackageStartupMessages({library(Matrix)})

WORK <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))

samples <- c("mix0hr_Cropseq","mix72hr_Cropseq","mix0hr_DirectCapture","mix72hr_DirectCapture")

# ---- authors' exact accuracy metric (binary_accuracy_metrics) ----
binary_accuracy_metrics <- function(pred, truth) {
  tp <- sum(pred & truth, na.rm=T); fp <- sum(pred & !truth, na.rm=T)
  tn <- sum(!pred & !truth, na.rm=T); fn <- sum(!pred & truth, na.rm=T)
  naTrue <- sum(is.na(pred) & truth); naFalse <- sum(is.na(pred) & !truth)
  tot <- tp+fp+tn+fn+naTrue+naFalse
  tot_true <- tp+fn+naTrue; tot_false <- tn+fp+naFalse
  stopifnot(tot == length(truth))
  acc_true <- tp/tot_true; acc_false <- tn/tot_false
  acc_stderr <- tot_true*acc_true*(1-acc_true) + tot_false*acc_false*(1-acc_false)
  acc_stderr <- sqrt(acc_stderr/tot^2)
  c(Accuracy=(tp+tn)/tot, Accuracy_stderr=acc_stderr,
    TP=tp, FP=fp, TN=tn, FN=fn, naTrue=naTrue, naFalse=naFalse, tot=tot)
}

# read a method assignment mtx (guides x cells) -> n_assigned per species per cell
read_assign <- function(path, guide_type) {
  m <- as(readMM(path), "CsparseMatrix")  # guides x cells, logical/numeric
  list(homo = Matrix::colSums(m[guide_type=="homo_guide", , drop=FALSE] > 0),
       mus  = Matrix::colSums(m[guide_type=="mus_guide",  , drop=FALSE] > 0))
}

out <- list()
for (s in samples) {
  meta   <- read.csv(file.path(WORK, paste0(s, "_meta.csv")))
  guides <- read.csv(file.path(WORK, paste0(s, "_guides.csv")))
  gt <- guides$guide_type

  # ---- QC fields (exactly as 82_plot) ----
  sum_gex   <- meta$homo_sum_gex + meta$mus_sum_gex
  frac_homo <- meta$homo_sum_gex / (meta$homo_sum_gex + meta$mus_sum_gex)
  features_gex <- meta$features_gex
  mito_frac <- meta$mito_sum / sum_gex
  is_empty  <- sum_gex <= 100
  orig_qc_pass <- (mito_frac < .15) & (features_gex <= 6000) & (sum_gex <= 20000) &
                  (1500 <= features_gex) & (3500 <= sum_gex) &
                  (frac_homo < .1 | frac_homo > .9)
  has_homo <- (frac_homo >= .1) & !is_empty   # GEX truth

  for (method in c("geomux","fishash")) {
    if (method == "fishash") {
      a <- read_assign(file.path(WORK, paste0(s, "_fishash.mtx")), gt)
      n_homo <- a$homo; n_mus <- a$mus
    } else {
      # geomux tsv -> per-cell assigned guides (guide_ids_original are 0-based row idx into 200 guides)
      gx <- read.delim(file.path(WORK, paste0(s, "_geomux.tsv")))
      n_homo <- integer(nrow(meta)); n_mus <- integer(nrow(meta))
      gx2 <- gx[gx$moi > 0 & nzchar(gx$guide_ids_original), , drop=FALSE]
      # cell_id is 0-based original cell index
      for (r in seq_len(nrow(gx2))) {
        ci <- gx2$cell_id[r] + 1L
        gids <- as.integer(strsplit(gx2$guide_ids_original[r], "\\|")[[1]]) + 1L
        th <- sum(gt[gids] == "homo_guide"); tm <- sum(gt[gids] == "mus_guide")
        n_homo[ci] <- n_homo[ci] + th; n_mus[ci] <- n_mus[ci] + tm
      }
    }
    ah <- n_homo > 0; am <- n_mus > 0
    atype <- ifelse(ah & am, "both", ifelse(ah, "homo", ifelse(am, "mus", "neither")))
    binz <- ah
    binz[atype == "both"] <- NA; binz[atype == "neither"] <- NA

    sel <- orig_qc_pass
    res <- binary_accuracy_metrics(binz[sel], has_homo[sel])
    out[[length(out)+1]] <- data.frame(sample=s, method=method,
        n_qc=sum(sel), Accuracy=round(res["Accuracy"],4), stderr=round(res["Accuracy_stderr"],4),
        both=sum(atype[sel]=="both"), neither=sum(atype[sel]=="neither"))
  }
}
df <- do.call(rbind, out); rownames(df) <- NULL
print(df, row.names=FALSE)
cat("\n--- Accuracy table (rows=method, cols=sample) ---\n")
library(stats)
wide <- reshape(df[,c("sample","method","Accuracy")], idvar="method", timevar="sample", direction="wide")
print(wide, row.names=FALSE)
