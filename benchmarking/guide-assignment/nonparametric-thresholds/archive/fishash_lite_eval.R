#!/usr/bin/env Rscript
# Does fishash's key idea -- iteratively re-estimating the NOISE to correct
# Simpson's paradox (Algorithm 1) -- help over the plain geomux-core Fisher test
# on our cohort? Compare ambient_test_assign at n_iter = 1 (geomux core),
# 3, 5 (-> simplified fishash) on barnyard (Table-2 metric) and sims (Jaccard).
# Also EXPORT each QC'd barnyard cohort (mtx + species/native metadata) so the
# REAL geomux/fishash can be run on the identical cohort.

suppressPackageStartupMessages({library(Matrix)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
GA   <- normalizePath(file.path(HERE, ".."))
BARN <- path.expand("~/data/external/liu-2025-cleanser/GSE272457")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
EXP  <- file.path(HERE, "results", "barnyard_cohort_export"); dir.create(EXP, recursive = TRUE, showWarnings = FALSE)
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

load_qc <- function(stub) {
  feat <- read.delim(gzfile(file.path(BARN, paste0(stub,"_features.tsv.gz"))), header=FALSE)
  bc   <- readLines(gzfile(file.path(BARN, paste0(stub,"_barcodes.tsv.gz"))))
  m <- as(Matrix::readMM(gzfile(file.path(BARN, paste0(stub,"_matrix.mtx.gz")))), "CsparseMatrix")
  g <- feat$V3=="Gene Expression"; hg <- g & grepl("^GRCh38_",feat$V1); mg <- g & grepl("^mm10_",feat$V1)
  mito <- g & grepl("_MT-|_mt-",feat$V2)
  gu <- Matrix::colSums(m[g,]); nf <- Matrix::colSums(m[g,]>0); mf <- Matrix::colSums(m[mito,])/pmax(gu,1)
  h <- Matrix::colSums(m[hg,]); mo <- Matrix::colSums(m[mg,]); fh <- h/(h+mo)
  keep <- mf<.15 & nf>=1500 & nf<=6000 & gu>=3500 & gu<=20000 & pmax(fh,1-fh)>=.9
  grna <- m[feat$V3=="CRISPR Guide Capture", keep]
  gid  <- feat$V1[feat$V3=="CRISPR Guide Capture"]; rownames(grna) <- gid; colnames(grna) <- bc[keep]
  gidx <- as.integer(sub("nt_","",gid))
  list(grna=grna, species=ifelse(fh>=.9,"human","mouse")[keep], native=ifelse(gidx<=100,"human","mouse")) }
t2 <- function(A,sp,nat){A<-as.matrix(A);co<-nat[row(A)]==sp[col(A)];mean(Matrix::colSums(A&co)>=1 & Matrix::colSums(A&!co)==0)}
jac <- function(A, truth, gm){A<-as(A,"RsparseMatrix");gn<-rownames(gm)
  mean(vapply(seq_along(gn),function(g){p<-as.numeric(A[g,])>0;tr<-as.numeric(truth[gn[g],])>0&as.numeric(gm[g,])>0
    u<-sum(p|tr);if(u==0)NA_real_ else sum(p&tr)/u},numeric(1)),na.rm=TRUE)}

samples <- c(`Cropseq mix0hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix",
             `Cropseq mix72hr`="GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix",
             `DirectCapture mix0hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix",
             `DirectCapture mix72hr`="GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix")

cat("===== BARNYARD: Table-2 accuracy vs n_iter (noise correction) at q=0.05 =====\n")
cat(sprintf("%-22s %8s %8s %8s\n","sample","iter1","iter3","iter5"))
for (sn in names(samples)) {
  d <- load_qc(samples[[sn]])
  accs <- sapply(c(1,3,5), function(ni) t2(ambient_test_assign(d$grna, q=0.05, model="hypergeometric", n_iter=ni)$assignment_matrix, d$species, d$native))
  cat(sprintf("%-22s %8.4f %8.4f %8.4f\n", sn, accs[1], accs[2], accs[3]))
  # export QC'd cohort for the real packages
  tag <- gsub("[^a-zA-Z0-9]","_", sn); od <- file.path(EXP, tag); dir.create(od, showWarnings=FALSE)
  Matrix::writeMM(d$grna, file.path(od, "guide_counts.mtx"))            # guides x cells
  write.csv(data.frame(barcode=colnames(d$grna), species=d$species), file.path(od,"cells.csv"), row.names=FALSE)
  write.csv(data.frame(guide=rownames(d$grna), native=d$native), file.path(od,"guides.csv"), row.names=FALSE)
}

cat("\n===== SIMS: Jaccard vs n_iter at q=0.05 =====\n")
cat(sprintf("%-26s %8s %8s %8s\n","sim","iter1","iter3","iter5"))
for (ds in c("sims_sum_2np_3p","sims_gasperini_calibrated","sims_sum_repeat_old")) {
  gm <- as(readRDS(file.path(DATA, ds, "sceptre/grna_matrix.rds")), "RsparseMatrix")
  truth <- as(readRDS(file.path(DATA, ds, "true_pert_matrix.rds")), "RsparseMatrix")
  js <- sapply(c(1,3,5), function(ni) jac(ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=ni)$assignment_matrix, truth, gm))
  cat(sprintf("%-26s %8.3f %8.3f %8.3f\n", ds, js[1], js[2], js[3]))
}
cat("\nExported QC'd barnyard cohorts to", EXP, "\n")
