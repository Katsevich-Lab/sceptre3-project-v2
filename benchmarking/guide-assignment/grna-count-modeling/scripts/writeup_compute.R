#!/usr/bin/env Rscript
# ============================================================================
# Compute backbone for the fishash / fishash+ / CLEANSER method-comparison writeup.
# Runs fishash (cell_margin="observed") and fishash+ (cell_margin="ambient") on every real dataset,
# plus CLEANSER's cheap <=2 ambient DEPTH (no MCMC). Real CLEANSER ASSIGNMENTS (MCMC) are produced
# separately by run_cleanser_batch.sh and read at render time from results/ambient_ceiling/cleanser_calls/.
#
# Saves to results/ambient_ceiling/writeup/ :
#   assign/<ds>_{fishash,fishashplus}.rds   sparse assignment matrices (for three-way Jaccard later)
#   depths_sampled.csv     per-cell depth (fishash rank-1, fishash+ rank-1, CLEANSER <=2), subsampled
#   depth_summary.csv      per-dataset median depths + spearman(f+ , cleanser)
#   assign_counts.csv      per-dataset per-method assignment counts + MOI
#   ambient_frac_by_count.csv   per (dataset,k,method): #ambient entries at count k / #entries at k
#   jaccard_ff.csv         per-dataset fishash-vs-fishash+ Jaccard (overall) + by count-bin
#   ceilings.csv           per-guide ambient ceiling (fishash, fishash+)
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats) })
source("scripts/contingency_method.R")
set.seed(1)
OUT <- "results/ambient_ceiling/writeup"; AOUT <- file.path(OUT,"assign")
dir.create(AOUT, recursive=TRUE, showWarnings=FALSE)

source("scripts/datasets.R")
paths <- dataset_paths()

Kmax <- 40L
depths <- list(); dsum <- list(); acount <- list(); afrac <- list(); jacc <- list(); ceil <- list()
for (nm in names(paths)) {
  t0 <- proc.time()["elapsed"]
  counts <- load_grna_matrix(paths[[nm]])
  G <- nrow(counts); C <- ncol(counts)
  cat(sprintf("[%s] %s %d x %d ...\n", format(Sys.time(),"%H:%M:%S"), nm, G, C)); flush.console()
  Af <- as(fishash_assign(counts)$assigned, "lgCMatrix")
  Ap <- as(fishashplus_assign(counts)$assigned, "lgCMatrix")
  saveRDS(Af, file.path(AOUT, paste0(nm,"_fishash.rds")))
  saveRDS(Ap, file.path(AOUT, paste0(nm,"_fishashplus.rds")))

  # per-cell depths: rank-1 (each method's signal-masked completion) + CLEANSER <=2
  bgf <- fishash::impute_masked_counts(counts, Af); bgp <- fishash::impute_masked_counts(counts, Ap)
  df <- Matrix::colSums(bgf); dp <- Matrix::colSums(bgp); dc <- kappa_leq2(counts)
  keep <- df>0 | dp>0 | dc>0; idx <- which(keep); if (length(idx) > 6000) idx <- sample(idx, 6000)
  depths[[nm]] <- data.frame(dataset=nm, fishash=df[idx], fishashplus=dp[idx], cleanser=dc[idx])
  pospos <- dp>0 & dc>0
  dsum[[nm]] <- data.frame(dataset=nm, guides=G, cells=C, med_lib=median(Matrix::colSums(counts)),
    med_depth_fishash=median(df[df>0]), med_depth_fishashplus=median(dp[dp>0]),
    med_depth_cleanser=median(dc[dc>0]),
    spearman_fp_clns=suppressWarnings(cor(dp[pospos], dc[pospos], method="spearman")))

  acount[[nm]] <- data.frame(dataset=nm, fishash=sum(Af), fishashplus=sum(Ap),
    moi_fishash=median(Matrix::colSums(Af)), moi_fishashplus=median(Matrix::colSums(Ap)))

  # ambient fraction by count: over nonzero entries, fraction UNASSIGNED at each count k
  ct <- as(counts,"TsparseMatrix"); i <- ct@i+1L; j <- ct@j+1L; y <- as.numeric(ct@x)
  af <- as.logical(Af[cbind(i,j)]); af[is.na(af)] <- FALSE
  ap <- as.logical(Ap[cbind(i,j)]); ap[is.na(ap)] <- FALSE
  kk <- pmin(y, Kmax)
  tot <- tapply(rep(1,length(kk)), kk, sum)
  amb_f <- tapply(!af, kk, sum); amb_p <- tapply(!ap, kk, sum)
  ks <- as.integer(names(tot))
  afrac[[nm]] <- rbind(
    data.frame(dataset=nm, k=ks, method="fishash",  n_ambient=as.numeric(amb_f), n_total=as.numeric(tot)),
    data.frame(dataset=nm, k=ks, method="fishash+", n_ambient=as.numeric(amb_p), n_total=as.numeric(tot)))

  # fishash vs fishash+ Jaccard (overall)
  inter <- sum(Af & Ap); uni <- sum(Af | Ap)
  jacc[[nm]] <- data.frame(dataset=nm, pair="fishash|fishash+", jaccard=ifelse(uni>0, inter/uni, 1),
    n_fishash=sum(Af), n_fishashplus=sum(Ap), n_intersect=inter, n_union=uni)

  # per-guide ambient ceilings
  ceil_f <- tapply(y[!af], factor(i[!af], levels=seq_len(G)), max)
  ceil_p <- tapply(y[!ap], factor(i[!ap], levels=seq_len(G)), max)
  ceil[[nm]] <- data.frame(dataset=nm, ceil_fishash=as.numeric(ceil_f), ceil_fishashplus=as.numeric(ceil_p))

  cat(sprintf("   fishash=%d fishash+=%d  jaccard=%.3f  (%.1fs)\n", sum(Af), sum(Ap),
    jacc[[nm]]$jaccard, proc.time()["elapsed"]-t0)); flush.console()
  saveRDS(list(depths=depths, dsum=dsum, acount=acount, afrac=afrac, jacc=jacc, ceil=ceil),
    file.path(OUT,"_progress.rds"))
  rm(counts, Af, Ap, bgf, bgp, ct); gc()
}
write.csv(do.call(rbind,depths), file.path(OUT,"depths_sampled.csv"), row.names=FALSE)
write.csv(do.call(rbind,dsum),   file.path(OUT,"depth_summary.csv"), row.names=FALSE)
write.csv(do.call(rbind,acount), file.path(OUT,"assign_counts.csv"), row.names=FALSE)
write.csv(do.call(rbind,afrac),  file.path(OUT,"ambient_frac_by_count.csv"), row.names=FALSE)
write.csv(do.call(rbind,jacc),   file.path(OUT,"jaccard_ff.csv"), row.names=FALSE)
write.csv(do.call(rbind,ceil),   file.path(OUT,"ceilings.csv"), row.names=FALSE)
cat("\nwriteup_compute done -> ", OUT, "\n")
