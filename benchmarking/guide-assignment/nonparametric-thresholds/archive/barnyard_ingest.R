#!/usr/bin/env Rscript
# Ingest the CLEANSER barnyard (species-mixing) data and score guide-assignment
# SPECIFICITY using its real ambient ground truth: a gRNA from one species'
# library detected in the OTHER species' cell is guaranteed ambient (a false
# positive). For each method we report the false-positive rate over all
# wrong-species (cell x gRNA) pairs -- the one thing no simulation can validate.
#
# Builds, per mix sample: sceptre/grna_matrix.rds (200 gRNAs x kept cells),
# meta.rds (per-cell species, per-gRNA native species), ambient_truth.rds.

suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
OUT  <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

samples <- c(
  lrb100_0hr  = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
  lrb100_72hr = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix",
  mch2_0hr    = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
  mch2_72hr   = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix"
)

ingest <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub, "_features.tsv.gz"))), header = FALSE)
  bc   <- readLines(gzfile(file.path(BARN, paste0(stub, "_barcodes.tsv.gz"))))
  m    <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub, "_matrix.mtx.gz")))), "CsparseMatrix")
  is_guide <- feat$V3 == "CRISPR Guide Capture"
  hg <- feat$V3 == "Gene Expression" & grepl("^GRCh38_", feat$V1)
  mg <- feat$V3 == "Gene Expression" & grepl("^mm10_",   feat$V1)
  # per-cell species from GENE rows only (never the gRNA rows)
  h  <- Matrix::colSums(m[hg, , drop = FALSE]); mo <- Matrix::colSums(m[mg, , drop = FALSE])
  fh <- h / (h + mo)
  species <- ifelse(fh > 0.9, "human", ifelse(fh < 0.1, "mouse", NA_character_))
  keep <- !is.na(species); species <- species[keep]
  grna <- m[is_guide, keep, drop = FALSE]
  rownames(grna) <- feat$V1[is_guide]; colnames(grna) <- bc[keep]
  gidx <- as.integer(sub("nt_", "", rownames(grna)))
  guide_native <- ifelse(gidx <= 100, "human", "mouse")   # lib1=human, lib2=mouse
  list(grna = grna, species = species, guide_native = guide_native)
}

thr_methods <- list(
  `otsu`       = function(c) otsu_threshold_log1p(c)$t,
  `valley`     = function(c) smoothed_valley_threshold(c)$t,
  `count>=1`   = function(c) 1, `count>=3` = function(c) 3, `count>=5` = function(c) 5
)

all_summ <- list(); diag_obj <- NULL
for (sn in names(samples)) {
  d <- ingest(samples[[sn]])
  ng <- nrow(d$grna); gm <- as(d$grna, "RsparseMatrix")
  is_mouse <- d$species == "mouse"
  # save framework artifacts
  od <- file.path(OUT, paste0("barnyard_", sn)); dir.create(file.path(od, "sceptre"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(d$grna, file.path(od, "sceptre", "grna_matrix.rds"))
  saveRDS(list(species = d$species, guide_native = d$guide_native), file.path(od, "meta.rds"))

  # per-guide: ambient cells = wrong-species cells; score each method's FP there
  rows <- list()
  for (g in seq_len(ng)) {
    counts  <- as.numeric(gm[g, ])
    ambient <- if (d$guide_native[g] == "human") is_mouse else !is_mouse   # wrong-species cells
    n_amb   <- sum(ambient); n_amb_nz <- sum(ambient & counts > 0)
    for (mn in names(thr_methods)) {
      t <- thr_methods[[mn]](counts); assigned <- is.finite(t) & counts >= t
      rows[[length(rows) + 1]] <- data.frame(method = mn,
        fp = sum(assigned & ambient), n_amb = n_amb,
        native_assigned = sum(assigned & !ambient), n_native = sum(!ambient))
    }
  }
  rr <- do.call(rbind, rows)
  summ <- aggregate(cbind(fp, n_amb, native_assigned, n_native) ~ method, rr, sum)
  summ$ambient_FPR        <- round(summ$fp / summ$n_amb, 4)
  summ$native_assign_rate <- round(summ$native_assigned / summ$n_native, 3)
  summ$sample <- sn
  all_summ[[sn]] <- summ
  cat(sprintf("\n[%s]  %d cells kept (%d human, %d mouse); ambient pairs/guide avg = %d\n",
              sn, length(d$species), sum(d$species == "human"), sum(is_mouse), round(mean(rr$n_amb[rr$method=="otsu"]))))
  print(summ[order(summ$ambient_FPR), c("method","ambient_FPR","native_assign_rate")], row.names = FALSE)
  if (sn == "lrb100_72hr") diag_obj <- d
}
summary_all <- do.call(rbind, all_summ)
write.csv(summary_all, file.path(HERE, "results", "barnyard_specificity.csv"), row.names = FALSE)

# ---- diagnostic: native vs ambient count distributions + thresholds ---------
d <- diag_obj; gm <- as(d$grna, "RsparseMatrix"); is_mouse <- d$species == "mouse"
sel <- c(which(d$guide_native == "human")[1:3], which(d$guide_native == "mouse")[1:3])
png(file.path(HERE, "results", "barnyard_native_vs_ambient.png"), width = 1300, height = 720, res = 115)
op <- par(mfrow = c(2, 3), mar = c(3.4, 3.4, 2.2, 0.6), mgp = c(2, 0.6, 0))
for (g in sel) {
  counts <- as.numeric(gm[g, ]); z <- log1p(counts)
  amb <- if (d$guide_native[g] == "human") is_mouse else !is_mouse
  br <- seq(0, max(z) + 1e-6, length.out = 40)
  hn <- hist(z[!amb], breaks = br, plot = FALSE); ha <- hist(z[amb], breaks = br, plot = FALSE)
  yn <- hn$counts; yn[yn == 0] <- NA; ya <- ha$counts; ya[ya == 0] <- NA
  plot(hn$mids, yn, type = "h", log = "y", lwd = 3, col = "grey55",
       xlim = c(0, max(z) + .3), ylim = c(1, max(hn$counts, ha$counts)),
       xlab = "log(1+count)", ylab = "cells", main = sprintf("%s (%s-native)", rownames(gm)[g], d$guide_native[g]), cex.main = .9)
  points(ha$mids, ya, type = "h", lwd = 3, col = adjustcolor("firebrick", .8))
  to <- otsu_threshold_log1p(counts)$t; tv <- smoothed_valley_threshold(counts)$t
  if (is.finite(to)) abline(v = log1p(to), col = "darkorange2", lwd = 2)
  if (is.finite(tv)) abline(v = log1p(tv), col = "royalblue", lwd = 2, lty = 2)
  legend("topright", bty = "n", cex = .7, legend = c("native cells", "ambient (wrong sp.)", "Otsu", "valley"),
         col = c("grey55", "firebrick", "darkorange2", "royalblue"), lwd = 2, lty = c(1,1,1,2))
}
par(op); dev.off()
cat("\nWrote results/barnyard_specificity.csv and results/barnyard_native_vs_ambient.png\n")
