#!/usr/bin/env Rscript
# Definitive same-cohort, same-metric comparison on the barnyard:
#   ours (plain BH-hypergeometric = geomux core) vs REAL geomux (default lor=10,
#   and lor=0) vs REAL fishash. Identical exported cohort, identical Table-2 metric.
suppressPackageStartupMessages({library(Matrix); library(fishash); library(SummarizedExperiment)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
EXP  <- file.path(HERE, "results", "barnyard_cohort_export")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

t2 <- function(cell_of, guide_of, cells, guides) {
  sp <- setNames(cells$species, cells$barcode); nat <- setNames(guides$native, guides$guide)
  correct <- as.logical(nat[as.character(guide_of)] == sp[as.character(cell_of)])
  ci <- match(as.character(cell_of), cells$barcode)
  ok <- !is.na(ci) & !is.na(correct)
  ci <- ci[ok]; correct <- correct[ok]
  ncor <- numeric(nrow(cells)); nwro <- numeric(nrow(cells))
  ac <- tapply(correct, ci, sum);  ncor[as.integer(names(ac))] <- ac
  aw <- tapply(!correct, ci, sum); nwro[as.integer(names(aw))] <- aw
  acc <- mean(ncor >= 1 & nwro == 0)
  amb <- if (length(correct)) mean(!correct) else NA
  c(acc = round(acc,4), amb_fdr = round(amb,4), se = round(sqrt(acc*(1-acc)/nrow(cells)),4),
    cells_assigned = sum(ncor+nwro>0))
}
order_samp <- c("Cropseq_mix0hr","Cropseq_mix72hr","DirectCapture_mix0hr","DirectCapture_mix72hr")
rows <- list()
for (sn in order_samp) {
  s <- file.path(EXP, sn)
  cells <- read.csv(file.path(s,"cells.csv")); guides <- read.csv(file.path(s,"guides.csv"))
  gm <- as(readMM(file.path(s,"guide_counts.mtx")), "CsparseMatrix"); rownames(gm)<-guides$guide; colnames(gm)<-cells$barcode

  # ours
  A <- as(ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=1)$assignment_matrix, "TsparseMatrix")
  rows[[length(rows)+1]] <- data.frame(sample=sn, method="ours (BH-hyper)",
    t(t2(cells$barcode[A@j+1L], guides$guide[A@i+1L], cells, guides)))
  # geomux default + nolor
  for (tag in c("default","nolor")) {
    g <- read.csv(file.path(s, paste0("geomux_",tag,"_calls.csv")))
    rows[[length(rows)+1]] <- data.frame(sample=sn, method=paste0("geomux ",tag),
      t(t2(as.character(g$cell_barcode), as.character(g$guide), cells, guides)))
  }
  # fishash
  res <- fishash(gm); cd <- as.data.frame(colData(res)); cd$barcode <- rownames(cd)
  asn <- strsplit(as.character(cd$assignment), ",")
  fc_cell <- rep(cd$barcode, lengths(asn)); fc_guide <- unlist(asn)
  keep <- !is.na(fc_guide) & fc_guide != ""
  rows[[length(rows)+1]] <- data.frame(sample=sn, method="fishash",
    t(t2(fc_cell[keep], fc_guide[keep], cells, guides)))
  cat("done", sn, "\n")
}
df <- do.call(rbind, rows)
cat("\n===== SAME-COHORT, SAME-METRIC: ours vs real geomux vs real fishash =====\n")
print(df, row.names = FALSE)
write.csv(df, file.path(HERE, "results", "final_same_cohort.csv"), row.names = FALSE)

# wide accuracy table
library(stats)
w <- reshape(df[,c("sample","method","acc")], idvar="sample", timevar="method", direction="wide")
names(w) <- sub("acc.","",names(w)); cat("\nAccuracy (wide):\n"); print(w, row.names=FALSE)
