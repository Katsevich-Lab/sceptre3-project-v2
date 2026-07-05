#!/usr/bin/env Rscript
# Finalize the ambient-test form: does the simple Poisson approx match the exact
# hypergeometric? does iterative background re-estimation (n_iter) help? Check on
# sims (Jaccard vs truth) and 2 barnyard samples (Table-2 acc). Also sanity-run
# on the real 2024-2026 survey datasets (no truth: assignment-rate / separation).

suppressPackageStartupMessages({library(Matrix)})
HERE <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
SURV <- path.expand("~/data/external/perturbseq-survey")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

jac_vs_truth <- function(A, truth, gm) {
  A <- as(A, "RsparseMatrix"); gn <- rownames(gm)
  mean(vapply(seq_along(gn), function(g) {
    pred <- as.numeric(A[g, ]) > 0; tr <- as.numeric(truth[gn[g], ]) > 0 & as.numeric(gm[g, ]) > 0
    u <- sum(pred | tr); if (u == 0) NA_real_ else sum(pred & tr) / u
  }, numeric(1)), na.rm = TRUE)
}
configs <- list(
  `hyper n1 q.01` = list(model="hypergeometric", n_iter=1, q=0.01),
  `hyper n2 q.01` = list(model="hypergeometric", n_iter=2, q=0.01),
  `pois  n2 q.01` = list(model="poisson",        n_iter=2, q=0.01),
  `hyper n2 q.05` = list(model="hypergeometric", n_iter=2, q=0.05))

cat("===== SIMS: Jaccard vs ground truth =====\n")
for (ds in c("sims_gasperini_calibrated", "sims_sum_2np_3p")) {
  gm <- as(readRDS(file.path(DATA, ds, "sceptre/grna_matrix.rds")), "RsparseMatrix")
  truth <- as(readRDS(file.path(DATA, ds, "true_pert_matrix.rds")), "RsparseMatrix")
  cat(sprintf("\n%s (%d guides):\n", ds, nrow(gm)))
  for (cn in names(configs)) { cf <- configs[[cn]]
    A <- ambient_test_assign(gm, q=cf$q, model=cf$model, n_iter=cf$n_iter)$assignment_matrix
    cat(sprintf("  %-14s jaccard=%.3f\n", cn, jac_vs_truth(A, truth, gm))) }
}

# barnyard (2 samples)
load_qc <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub,"_features.tsv.gz"))), header=FALSE)
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub,"_matrix.mtx.gz")))), "CsparseMatrix")
  g <- feat$V3=="Gene Expression"; hg <- g & grepl("^GRCh38_",feat$V1); mg <- g & grepl("^mm10_",feat$V1)
  mito <- g & grepl("_MT-|_mt-",feat$V2)
  gu <- Matrix::colSums(m[g,]); nf <- Matrix::colSums(m[g,]>0); mf <- Matrix::colSums(m[mito,])/pmax(gu,1)
  h <- Matrix::colSums(m[hg,]); mo <- Matrix::colSums(m[mg,]); fh <- h/(h+mo)
  keep <- mf<.15 & nf>=1500 & nf<=6000 & gu>=3500 & gu<=20000 & pmax(fh,1-fh)>=.9
  gidx <- as.integer(sub("nt_","",feat$V1[feat$V3=="CRISPR Guide Capture"]))
  list(grna=m[feat$V3=="CRISPR Guide Capture",keep], species=ifelse(fh>=.9,"human","mouse")[keep],
       native=ifelse(gidx<=100,"human","mouse")) }
t2 <- function(A,sp,nat){A<-as.matrix(A);co<-nat[row(A)]==sp[col(A)];mean(Matrix::colSums(A&co)>=1 & Matrix::colSums(A&!co)==0)}
cat("\n===== BARNYARD: Table-2 accuracy =====\n")
for (st in c(`Cropseq mix0hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
             `DirectCapture mix0hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix")) {
  d <- load_qc(st); cat(sprintf("\n%s:\n", names(which(c(`Cropseq mix0hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",`DirectCapture mix0hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix")==st))))
  for (cn in names(configs)) { cf <- configs[[cn]]
    A <- ambient_test_assign(d$grna, q=cf$q, model=cf$model, n_iter=cf$n_iter)$assignment_matrix
    cat(sprintf("  %-14s acc=%.4f\n", cn, t2(A, d$species, d$native))) }
}

# survey datasets sanity (no truth): n cells, n guides, assignment rate, median guides/cell
cat("\n===== SURVEY (real 2024-2026, no ground truth): ambient test q=0.01 sanity =====\n")
for (dd in list.dirs(SURV, recursive=FALSE)) {
  f <- file.path(dd, "grna_matrix.rds"); if (!file.exists(f)) next
  gm <- readRDS(f); A <- ambient_test_assign(gm, q=0.01, model="poisson")$assignment_matrix
  cat(sprintf("  %-34s %5d guides x %6d cells | assigned/cell median=%.1f mean=%.2f | %% cells w/>=1 guide=%.1f\n",
      basename(dd), nrow(gm), ncol(gm), median(Matrix::colSums(A)), mean(Matrix::colSums(A)),
      100*mean(Matrix::colSums(A) >= 1)))
}
cat("\ndone\n")
