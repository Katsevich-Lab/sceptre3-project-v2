#!/usr/bin/env Rscript
# Faithful(er) fishash: the noise-conditioned Fisher test that corrects Simpson's
# paradox. For pair (g,c) the 2x2 table replaces the "other cells" counts with
# NOISE counts (uncalled background), iterated per Algorithm 1. Derivation gives
#   white = noise_row_g,  draws = n_c,  N* = n_c + T_noise - noise_col_c
#   p = P(Hypergeom(white, N*-white, n_c) >= y).
# Compare to the plain geomux-core test (ambient_test_assign, n_iter=1).

suppressPackageStartupMessages(library(Matrix))
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

fishash_test_assign <- function(grna, q = 0.05, niter = 6) {
  tg <- as(grna, "TsparseMatrix"); i <- tg@i + 1L; j <- tg@j + 1L; y <- tg@x
  n_c <- Matrix::colSums(grna); G <- nrow(grna); C <- ncol(grna)
  called <- rep(FALSE, length(y)); Omega <- rep(TRUE, length(y)); p <- rep(1, length(y))
  for (it in seq_len(niter)) {
    noise <- Omega
    nr <- as.numeric(tapply(y[noise], factor(i[noise], levels = seq_len(G)), sum)); nr[is.na(nr)] <- 0
    nc <- as.numeric(tapply(y[noise], factor(j[noise], levels = seq_len(C)), sum)); nc[is.na(nc)] <- 0
    Tn <- sum(y[noise])
    white <- nr[i]; Nstar <- n_c[j] + Tn - nc[j]
    p <- stats::phyper(y - 1, white, pmax(Nstar - white, 0), n_c[j], lower.tail = FALSE)
    called <- stats::p.adjust(p, "BH") < q; called[is.na(called)] <- FALSE
    Omega <- if (it <= 3) !called else Omega & !called
  }
  Matrix::sparseMatrix(i = i[called], j = j[called], x = TRUE, dims = dim(grna), dimnames = dimnames(grna))
}

# barnyard QC + metric
load_qc <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub,"_features.tsv.gz"))), header=FALSE)
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub,"_matrix.mtx.gz")))), "CsparseMatrix")
  g <- feat$V3=="Gene Expression"; hg<-g&grepl("^GRCh38_",feat$V1); mg<-g&grepl("^mm10_",feat$V1); mito<-g&grepl("_MT-|_mt-",feat$V2)
  gu<-Matrix::colSums(m[g,]); nf<-Matrix::colSums(m[g,]>0); mf<-Matrix::colSums(m[mito,])/pmax(gu,1)
  h<-Matrix::colSums(m[hg,]); mo<-Matrix::colSums(m[mg,]); fh<-h/(h+mo)
  keep<-mf<.15&nf>=1500&nf<=6000&gu>=3500&gu<=20000&pmax(fh,1-fh)>=.9
  gid<-feat$V1[feat$V3=="CRISPR Guide Capture"]
  list(grna=m[feat$V3=="CRISPR Guide Capture",keep], species=ifelse(fh>=.9,"human","mouse")[keep],
       native=ifelse(as.integer(sub("nt_","",gid))<=100,"human","mouse")) }
t2 <- function(A,sp,nat){A<-as.matrix(A);co<-nat[row(A)]==sp[col(A)]
  list(acc=mean(Matrix::colSums(A&co)>=1 & Matrix::colSums(A&!co)==0),
       amb_fdr=if(sum(A)>0) sum(A&!co)/sum(A) else NA)}
jac <- function(A,truth,gm){A<-as(A,"RsparseMatrix");gn<-rownames(gm)
  v<-vapply(seq_along(gn),function(g){pr<-as.numeric(A[g,])>0;tr<-as.numeric(truth[gn[g],])>0&as.numeric(gm[g,])>0
    pu<-sum(pr|tr);c(jac=if(pu==0)NA else sum(pr&tr)/pu, prec=if(sum(pr))sum(pr&tr)/sum(pr) else NA, rec=if(sum(tr))sum(pr&tr)/sum(tr) else NA)},numeric(3))
  rowMeans(v,na.rm=TRUE)}

samples <- c(`Cropseq 0hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
             `Cropseq 72hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
             `DirectCap 0hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
             `DirectCap 72hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix")
cat("===== BARNYARD: plain (geomux-core) vs faithful fishash (q=0.05) =====\n")
cat(sprintf("%-15s | %-18s | %-18s\n","sample","plain acc (ambFDR)","fishash acc (ambFDR)"))
for (sn in names(samples)) { d <- load_qc(samples[[sn]])
  pl <- t2(ambient_test_assign(d$grna,q=0.05,n_iter=1)$assignment_matrix, d$species,d$native)
  fa <- t2(fishash_test_assign(d$grna,q=0.05), d$species,d$native)
  cat(sprintf("%-15s | %.4f (%.4f)    | %.4f (%.4f)\n", sn, pl$acc, pl$amb_fdr, fa$acc, fa$amb_fdr)) }

cat("\n===== SIMS: plain vs fishash (Jaccard / precision / recall, q=0.05) =====\n")
for (ds in c("sims_sum_2np_3p","sims_gasperini_calibrated","sims_sum_repeat_old")) {
  gm <- as(readRDS(file.path(DATA,ds,"sceptre/grna_matrix.rds")),"RsparseMatrix")
  truth <- as(readRDS(file.path(DATA,ds,"true_pert_matrix.rds")),"RsparseMatrix")
  pl <- jac(ambient_test_assign(gm,q=0.05,n_iter=1)$assignment_matrix, truth, gm)
  fa <- jac(fishash_test_assign(gm,q=0.05), truth, gm)
  cat(sprintf("%-26s plain J=%.3f(p%.2f/r%.2f)  fishash J=%.3f(p%.2f/r%.2f)\n", ds,
      pl["jac"],pl["prec"],pl["rec"], fa["jac"],fa["prec"],fa["rec"])) }
cat("\ndone\n")
