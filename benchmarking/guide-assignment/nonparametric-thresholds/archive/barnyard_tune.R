#!/usr/bin/env Rscript
# Can the nonparametric methods under-assign less if we tune their threshold
# knob instead of using one fixed setting? We (1) decompose valley's failures
# (under- vs over-assignment, abstention rate), and (2) sweep two levers:
#   - fallback: when valley abstains, use Otsu instead of assigning nothing
#   - shift:    lower the integer cut by `shift` (floored at 1) -> assign more
# Metric = Fishash Table 2 accuracy (>=1 correct-species guide & 0 wrong).

suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

samples <- c(
  "Cropseq mix0hr"        = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
  "Cropseq mix72hr"       = "GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
  "DirectCapture mix0hr"  = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
  "DirectCapture mix72hr" = "GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix")

load_qc <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub, "_features.tsv.gz"))), header = FALSE)
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub, "_matrix.mtx.gz")))), "CsparseMatrix")
  g  <- feat$V3 == "Gene Expression"
  hg <- g & grepl("^GRCh38_", feat$V1); mg <- g & grepl("^mm10_", feat$V1)
  mito <- g & grepl("_MT-|_mt-", feat$V2)
  gu <- Matrix::colSums(m[g, ]); nf <- Matrix::colSums(m[g, ] > 0)
  mf <- Matrix::colSums(m[mito, ]) / pmax(gu, 1)
  h <- Matrix::colSums(m[hg, ]); mo <- Matrix::colSums(m[mg, ]); fh <- h / (h + mo)
  keep <- mf < .15 & nf >= 1500 & nf <= 6000 & gu >= 3500 & gu <= 20000 & pmax(fh, 1 - fh) >= .9
  grna <- as(m[feat$V3 == "CRISPR Guide Capture", keep], "RsparseMatrix")
  gidx <- as.integer(sub("nt_", "", feat$V1[feat$V3 == "CRISPR Guide Capture"]))
  list(grna = grna, species = ifelse(fh >= .9, "human", "mouse")[keep],
       guide_native = ifelse(gidx <= 100, "human", "mouse"))
}

# per-guide cut from base method; valley may abstain (Inf)
guide_cuts <- function(grna) {
  t(vapply(seq_len(nrow(grna)), function(g) {
    counts <- as.numeric(grna[g, ])
    c(otsu = otsu_threshold_log1p(counts)$t, valley = smoothed_valley_threshold(counts)$t)
  }, numeric(2)))
}

# accuracy given a per-guide threshold vector and a shift
acc_for <- function(grna, species, native, thr, shift) {
  n_cell <- length(species); A <- matrix(FALSE, nrow(grna), n_cell)
  for (g in seq_len(nrow(grna))) {
    t <- thr[g]; if (!is.finite(t)) next
    A[g, ] <- as.numeric(grna[g, ]) >= max(1, t - shift)
  }
  correct <- native[row(A)] == species[col(A)]
  ncor <- Matrix::colSums(A & correct); nwro <- Matrix::colSums(A & !correct)
  mean(ncor >= 1 & nwro == 0)
}

shifts <- 0:3
rows <- list(); decomp <- list()
for (sn in names(samples)) {
  d <- load_qc(samples[[sn]]); cuts <- guide_cuts(d$grna)
  t_otsu <- cuts[, "otsu"]; t_val <- cuts[, "valley"]
  t_val_fb <- ifelse(is.finite(t_val), t_val, t_otsu)         # fallback to Otsu on abstain
  abstain_rate <- mean(!is.finite(t_val))

  # failure decomposition for the ORIGINAL valley (shift 0, hard abstain)
  A <- matrix(FALSE, nrow(d$grna), length(d$species))
  for (g in seq_len(nrow(d$grna))) if (is.finite(t_val[g])) A[g, ] <- as.numeric(d$grna[g, ]) >= t_val[g]
  correct <- d$guide_native[row(A)] == d$species[col(A)]
  ncor <- Matrix::colSums(A & correct); nwro <- Matrix::colSums(A & !correct)
  fail <- !(ncor >= 1 & nwro == 0)
  decomp[[sn]] <- data.frame(sample = sn, abstain_rate = round(abstain_rate, 3),
    fail_total = round(mean(fail), 3),
    under_only = round(mean(fail & ncor == 0 & nwro == 0), 3),
    over_only  = round(mean(fail & ncor >= 1 & nwro > 0), 3),
    both       = round(mean(fail & ncor == 0 & nwro > 0), 3))

  for (s in shifts) {
    rows[[length(rows)+1]] <- data.frame(sample = sn, shift = s,
      otsu          = acc_for(d$grna, d$species, d$guide_native, t_otsu,   s),
      valley        = acc_for(d$grna, d$species, d$guide_native, t_val,    s),
      `valley+fb`   = acc_for(d$grna, d$species, d$guide_native, t_val_fb, s),
      check.names = FALSE)
  }
  cat("done", sn, "(abstain", round(abstain_rate,3), ")\n")
}
dec <- do.call(rbind, decomp); df <- do.call(rbind, rows)
cat("\n=== valley failure decomposition (shift 0, hard abstain) ===\n"); print(dec, row.names = FALSE)
cat("\n=== Table-2 accuracy vs shift ===\n")
print(transform(df, otsu=round(otsu,4), valley=round(valley,4), `valley+fb`=round(`valley+fb`,4),
                check.names=FALSE), row.names = FALSE)
write.csv(df, file.path(HERE, "results", "barnyard_tune.csv"), row.names = FALSE)

long <- reshape(df, varying = c("otsu","valley","valley+fb"), v.names = "acc",
                timevar = "method", times = c("otsu","valley","valley+fb"), direction = "long")
p <- ggplot(long, aes(shift, acc, color = method)) +
  geom_line() + geom_point(size = 1.8) + facet_wrap(~sample, nrow = 1) +
  labs(title = "Tuning the cut reduces under-assignment (barnyard, Fishash Table-2 metric)",
       subtitle = "shift = lower the integer cut by this much (floored at 1); +fb = Otsu fallback on abstain",
       x = "threshold downshift", y = "accuracy") +
  theme_bw(base_size = 9) + theme(legend.position = "top")
ggsave(file.path(HERE, "results", "barnyard_tune.png"), p, width = 11, height = 3.6, dpi = 120)
cat("\nWrote results/barnyard_tune.png\n")
