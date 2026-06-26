#!/usr/bin/env Rscript
# Evaluate the principled ambient-proportion test (hypergeometric / geomux-family,
# FDR-controlled per-(cell,guide) threshold) against Otsu and the smoothed valley.
#   (A) Barnyard: Fishash Table-2 accuracy + REALIZED ambient false-discovery
#       fraction vs the nominal FDR q (a calibration check the ground truth allows).
#   (B) Sims (ground truth): per-guide precision / recall / Jaccard.

suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

# ---------- (A) BARNYARD ------------------------------------------------------
load_qc <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub, "_features.tsv.gz"))), header = FALSE)
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub, "_matrix.mtx.gz")))), "CsparseMatrix")
  g <- feat$V3 == "Gene Expression"
  hg <- g & grepl("^GRCh38_", feat$V1); mg <- g & grepl("^mm10_", feat$V1)
  mito <- g & grepl("_MT-|_mt-", feat$V2)
  gu <- Matrix::colSums(m[g, ]); nf <- Matrix::colSums(m[g, ] > 0)
  mf <- Matrix::colSums(m[mito, ]) / pmax(gu, 1)
  h <- Matrix::colSums(m[hg, ]); mo <- Matrix::colSums(m[mg, ]); fh <- h / (h + mo)
  keep <- mf < .15 & nf >= 1500 & nf <= 6000 & gu >= 3500 & gu <= 20000 & pmax(fh, 1 - fh) >= .9
  grna <- m[feat$V3 == "CRISPR Guide Capture", keep]
  gidx <- as.integer(sub("nt_", "", feat$V1[feat$V3 == "CRISPR Guide Capture"]))
  list(grna = grna, species = ifelse(fh >= .9, "human", "mouse")[keep],
       native = ifelse(gidx <= 100, "human", "mouse"))
}
table2 <- function(A, species, native) {        # >=1 correct-species guide & 0 wrong
  A <- as.matrix(A); correct <- native[row(A)] == species[col(A)]
  nc <- Matrix::colSums(A & correct); nw <- Matrix::colSums(A & !correct)
  list(acc = mean(nc >= 1 & nw == 0), n = ncol(A),
       amb_fp = if (sum(A) > 0) sum(A & !correct) / sum(A) else NA_real_,  # realized ambient-FDR
       calls_per_cell = mean(Matrix::colSums(A)))
}
samples <- c(`Cropseq mix0hr` = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
             `Cropseq mix72hr` = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
             `DirectCapture mix0hr` = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
             `DirectCapture mix72hr` = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix")

barn <- list()
for (sn in names(samples)) {
  d <- load_qc(samples[[sn]])
  A_o <- assign_by_threshold(d$grna, otsu_threshold_log1p)$assignment_matrix
  A_v <- assign_by_threshold(d$grna, smoothed_valley_threshold)$assignment_matrix
  meth <- list(otsu = A_o, valley = A_v)
  for (qq in c(0.01, 0.05, 0.10))
    meth[[sprintf("ambient(q=%.2f)", qq)]] <- ambient_test_assign(d$grna, q = qq, model = "hypergeometric")$assignment_matrix
  for (mn in names(meth)) {
    r <- table2(meth[[mn]], d$species, d$native)
    barn[[length(barn)+1]] <- data.frame(sample = sn, method = mn,
      acc = round(r$acc, 4), se = round(sqrt(r$acc*(1-r$acc)/r$n), 4),
      realized_ambient_FDR = round(r$amb_fp, 4), calls_per_cell = round(r$calls_per_cell, 2))
  }
  cat("barnyard done:", sn, "\n")
}
barn_df <- do.call(rbind, barn)
cat("\n===== BARNYARD: Table-2 accuracy + realized ambient-FDR =====\n")
print(barn_df, row.names = FALSE)
write.csv(barn_df, file.path(HERE, "results", "principled_barnyard.csv"), row.names = FALSE)

# ---------- (B) SIMS ---------------------------------------------------------
sims <- list(
  `Repl 0.5%pert` = "replogle-rd7_simulated_100k_0.005-pert",
  `Repl 1.5%pert` = "replogle-rd7_simulated_100k_0.015-pert",
  `Gasp-calib`    = "sims_gasperini_calibrated")
pgrid <- function(A, truth, gm) {                # per-guide precision/recall/jaccard
  A <- as(A, "RsparseMatrix"); gnames <- rownames(gm)
  do.call(rbind, lapply(seq_along(gnames), function(g) {
    counts <- as.numeric(gm[g, ]); pred <- as.numeric(A[g, ]) > 0
    tr <- as.numeric(truth[gnames[g], ]) > 0 & counts > 0
    tp <- sum(pred & tr)
    data.frame(prec = if (sum(pred)) tp/sum(pred) else NA, rec = if (sum(tr)) tp/sum(tr) else NA,
               jac = if (sum(pred|tr)) tp/sum(pred|tr) else NA)
  }))
}
sim_rows <- list()
for (sn in names(sims)) {
  base <- file.path(DATA, sims[[sn]])
  gm <- as(readRDS(file.path(base, "sceptre", "grna_matrix.rds")), "RsparseMatrix")
  tf <- if (file.exists(file.path(base, "true_pert_matrix.rds"))) "true_pert_matrix.rds" else "true_perturbation_status.rds"
  truth <- as(readRDS(file.path(base, tf)), "RsparseMatrix")
  meth <- list(otsu = assign_by_threshold(gm, otsu_threshold_log1p)$assignment_matrix,
               valley = assign_by_threshold(gm, smoothed_valley_threshold)$assignment_matrix,
               `ambient(q=0.05)` = ambient_test_assign(gm, q = 0.05, model = "hypergeometric")$assignment_matrix)
  for (mn in names(meth)) {
    pg <- pgrid(meth[[mn]], truth, gm)
    sim_rows[[length(sim_rows)+1]] <- data.frame(dataset = sn, method = mn,
      precision = round(mean(pg$prec, na.rm=TRUE), 3), recall = round(mean(pg$rec, na.rm=TRUE), 3),
      jaccard = round(mean(pg$jac, na.rm=TRUE), 3))
  }
  cat("sim done:", sn, "\n")
}
sim_df <- do.call(rbind, sim_rows)
cat("\n===== SIMS: per-guide metrics vs ground truth =====\n")
print(sim_df, row.names = FALSE)
write.csv(sim_df, file.path(HERE, "results", "principled_sims.csv"), row.names = FALSE)

# ---------- plot: barnyard Table-2 accuracy by method ------------------------
pd <- barn_df; pd$method <- factor(pd$method, levels = unique(pd$method))
p <- ggplot(pd, aes(method, acc, fill = grepl("ambient", method))) +
  geom_col() + geom_errorbar(aes(ymin = acc - se, ymax = acc + se), width = .3) +
  facet_wrap(~sample, nrow = 1) + coord_cartesian(ylim = c(0.7, 1)) +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#2c7fb8"), guide = "none") +
  labs(title = "Principled ambient test vs Otsu/valley (barnyard, Table-2 metric)",
       subtitle = "blue = ambient-proportion test (hypergeometric, FDR q)", x = NULL, y = "accuracy") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(HERE, "results", "principled_barnyard.png"), p, width = 12, height = 4, dpi = 120)
cat("\nWrote results/principled_barnyard.{csv,png} and principled_sims.csv\n")
