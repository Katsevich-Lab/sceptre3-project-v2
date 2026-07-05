#!/usr/bin/env Rscript
# Compute the Table 2 metric of Kamm/Yeung/Forrest (Fishash, bioRxiv 701179v2)
# for the two nonparametric thresholds on the Liu et al. barnyard data:
#   accuracy = fraction of cells with >=1 guide from the CORRECT species and
#              0 guides from the INCORRECT species.
# Replicates their cell QC (mito<15%, 1500-6000 features, 3500-20000 gene UMIs,
# >=90% UMIs from one species) so numbers are comparable to their table.
# Columns: CROP-seq = MCH2 chemistry; DirectCapture = LRB100 chemistry.

suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

samples <- c(
  "Cropseq mix0hr"        = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
  "Cropseq mix72hr"       = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
  "DirectCapture mix0hr"  = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
  "DirectCapture mix72hr" = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix"
)
methods <- list(otsu = function(c) otsu_threshold_log1p(c)$t,
                valley = function(c) smoothed_valley_threshold(c)$t)

fmt <- function(p, n) sprintf("%.4f ± %.4f", p, sqrt(p * (1 - p) / n))
res <- list()
for (sn in names(samples)) {
  stub <- samples[[sn]]
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub, "_features.tsv.gz"))), header = FALSE)
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub, "_matrix.mtx.gz")))), "CsparseMatrix")
  is_guide <- feat$V3 == "CRISPR Guide Capture"
  hg <- feat$V3 == "Gene Expression" & grepl("^GRCh38_", feat$V1)
  mg <- feat$V3 == "Gene Expression" & grepl("^mm10_",   feat$V1)
  gene <- feat$V3 == "Gene Expression"
  mito <- gene & grepl("_MT-|_mt-", feat$V2)

  gene_umis <- Matrix::colSums(m[gene, , drop = FALSE])
  n_feat    <- Matrix::colSums(m[gene, , drop = FALSE] > 0)
  mito_frac <- Matrix::colSums(m[mito, , drop = FALSE]) / pmax(gene_umis, 1)
  h <- Matrix::colSums(m[hg, , drop = FALSE]); mo <- Matrix::colSums(m[mg, , drop = FALSE])
  fh <- h / (h + mo); purity <- pmax(fh, 1 - fh)
  keep <- mito_frac < 0.15 & n_feat >= 1500 & n_feat <= 6000 &
          gene_umis >= 3500 & gene_umis <= 20000 & purity >= 0.90
  species <- ifelse(fh >= 0.9, "human", "mouse")[keep]

  grna <- as(m[is_guide, keep, drop = FALSE], "RsparseMatrix")
  gidx <- as.integer(sub("nt_", "", feat$V1[is_guide]))
  guide_native <- ifelse(gidx <= 100, "human", "mouse")
  n_cell <- length(species)

  for (mn in names(methods)) {
    # per-guide threshold -> assignment matrix (guides x cells)
    A <- matrix(FALSE, nrow(grna), n_cell)
    for (g in seq_len(nrow(grna))) {
      counts <- as.numeric(grna[g, ]); t <- methods[[mn]](counts)
      if (is.finite(t)) A[g, ] <- counts >= t
    }
    correct_sp <- guide_native[row(A)] == species[col(A)]   # same-species pairs
    n_correct <- Matrix::colSums(A & correct_sp)             # per cell: # right-species guides
    n_wrong   <- Matrix::colSums(A & !correct_sp)            # per cell: # wrong-species guides
    acc <- mean(n_correct >= 1 & n_wrong == 0)
    res[[paste(mn, sn)]] <- data.frame(method = mn, sample = sn, n = n_cell,
                                       acc = acc, txt = fmt(acc, n_cell))
  }
  cat(sprintf("[%s] %d cells pass QC (%d human, %d mouse)\n",
              sn, n_cell, sum(species == "human"), sum(species == "mouse")))
}
df <- do.call(rbind, res)
wide <- reshape(df[, c("method", "sample", "txt")], idvar = "method",
                timevar = "sample", direction = "wide")
names(wide) <- sub("txt.", "", names(wide))
cat("\n===== Table-2 metric (accuracy of gRNA species assignment) =====\n")
print(wide[, c("method", names(samples))], row.names = FALSE)
write.csv(df, file.path(HERE, "results", "barnyard_table2_metric.csv"), row.names = FALSE)
