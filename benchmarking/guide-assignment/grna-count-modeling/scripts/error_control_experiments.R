#!/usr/bin/env Rscript
# ============================================================================
# error_control_experiments.R  --  evidence behind error_control_grna_assignment.qmd
#
# Five self-contained experiments on the error-control question for gRNA assignment.
# All run the SAME method (fishash+ = contingency_assign, cell_margin="ambient",
# hypergeometric tail) and vary only the multiple-testing / QC treatment.
# Run from the grna-count-modeling/ folder root:
#     Rscript scripts/error_control_experiments.R all          # everything (~10 min)
#     Rscript scripts/error_control_experiments.R bh_vs_gs      # one experiment
#     ... {bh_vs_gs | poscontrol | perguide | calibration | singleton_qc}
# Outputs -> results/error_control/.
#
# Experiment map (writeup section in brackets):
#   bh_vs_gs      [3] realized FDR of BH vs GS vs BY on fishash+ p-values (sim truth)
#   poscontrol    [3] on-target power BH vs GS on DC-TAP, knockdown-validated
#   perguide      [4] pooled-BH vs GS vs per-guide-BH: which error UNIT (sim truth)
#   calibration   [6] does the assignment threshold move gRNA-level DE Type-I? (DC-TAP NTC)
#   singleton_qc  [7] low-MOI singleton QC: throw-out vs FDR, collision decomposition (sim)
# ============================================================================
suppressPackageStartupMessages({
  library(fishash); library(Matrix); library(SummarizedExperiment)
  library(sparseMatrixStats); library(matrixStats); library(extraDistr)
})
source("scripts/contingency_method.R"); source("scripts/datasets.R")
OUT <- "results/error_control"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Q <- 0.05

## ---- instrumented fishash+ : BH/GS/BY branches; returns the final p-values -----------------
# BH and GS branches are byte-identical to scripts/contingency_method.R; adds "BY" + returns logp.
cassign <- function(counts, q = 0.05, refit = 10, min_count = 2,
                    cell_margin = "ambient", tail = "hyper", fdr = "GS", init_assigned = NULL) {
  counts <- as(counts, "CsparseMatrix"); obs_col <- Matrix::colSums(counts)
  n_rows <- sum(Matrix::rowSums(counts) > 0); n_cols <- sum(obs_col > 0)
  n_entries <- as.numeric(n_rows) * as.numeric(n_cols)
  ct <- as(counts, "TsparseMatrix"); i <- ct@i + 1L; j <- ct@j + 1L; y <- as.numeric(ct@x)
  mask <- init_assigned; prev <- NULL; assigned <- NULL; logp <- NULL
  for (it in seq_len(refit + 1)) {
    background <- if (it == 1 && is.null(init_assigned)) counts else fishash::impute_masked_counts(counts, mask)
    nr <- Matrix::rowSums(background); nc <- Matrix::colSums(background); Tn <- sum(background)
    bgc <- background[cbind(i, j)]; K <- nr[i] - bgc + y
    if (cell_margin == "ambient") { draws <- nc[j] - bgc + y; pop <- Tn - bgc + y }
    else                          { draws <- obs_col[j];      pop <- Tn - nc[j] + obs_col[j] }
    logp <- stats::phyper(y - 1, K, pop - K, draws, lower.tail = FALSE, log.p = TRUE)
    if (fdr == "BH")      A <- stats::p.adjust(exp(logp), "BH") < q
    else if (fdr == "BY") A <- stats::p.adjust(exp(logp), "BY") < q
    else {                                                     # Guo-Sarkar block (per-cell)
      ml <- sparseMatrix(i = i, j = j, x = logp, dims = dim(counts)); colmin <- sparseMatrixStats::colMins(ml)
      B <- sum(stats::p.adjust(pmin(exp(colmin) * n_rows, 1), "BH") <= q)
      A <- logp <= (log(q) - log(n_entries) + log(max(B, 1)))
    }
    A <- A & !is.na(A) & (y >= min_count)
    assigned <- sparseMatrix(i = i[A], j = j[A], x = TRUE, dims = dim(counts), dimnames = dimnames(counts))
    mask <- if (it > 3 && !is.null(mask)) (mask | assigned) else assigned
    if (it > 1 && sum(abs(prev - assigned)) == 0) break
    prev <- assigned
  }
  list(assigned = assigned, logp = logp, i = i, j = j, y = y,
       n_rows = n_rows, n_entries = n_entries, dims = dim(counts))
}
# apply a correction to a FIXED p-value set (isolation): identical p-values in, correction varies
apply_corr <- function(fit, fdr, q = 0.05, min_count = 2) {
  logp <- fit$logp; i <- fit$i; j <- fit$j; y <- fit$y
  if (fdr == "BH")      A <- stats::p.adjust(exp(logp), "BH") < q
  else if (fdr == "BY") A <- stats::p.adjust(exp(logp), "BY") < q
  else { ml <- sparseMatrix(i = i, j = j, x = logp, dims = fit$dims); colmin <- sparseMatrixStats::colMins(ml)
    B <- sum(stats::p.adjust(pmin(exp(colmin) * fit$n_rows, 1), "BH") <= q)
    A <- logp <= (log(q) - log(fit$n_entries) + log(max(B, 1))) }
  A <- A & !is.na(A) & (y >= min_count)
  sparseMatrix(i = i[A], j = j[A], x = TRUE, dims = fit$dims)
}
# per-guide BH: BH applied SEPARATELY within each guide's cells (across the ~independent cell axis)
apply_perguide <- function(fit, q = 0.05, min_count = 2) {
  p <- exp(fit$logp); i <- fit$i; A <- logical(length(p))
  for (g in unique(i)) { idx <- which(i == g); A[idx] <- p.adjust(p[idx], "BH") < q }
  A <- A & (fit$y >= min_count)
  sparseMatrix(i = fit$i[A], j = fit$j[A], x = TRUE, dims = fit$dims)
}
fdp <- function(assigned, truth) {
  A <- as(assigned, "lgCMatrix"); Tr <- as(truth > 0, "lgCMatrix"); tpm <- A & Tr
  na <- sum(A); tp <- sum(tpm); nt <- sum(Tr)
  c(fdp = if (na > 0) (na - tp)/na else 0, recall = if (nt > 0) tp/nt else NA, n_assigned = na)
}
sim1 <- function(n_guides, n_cells, moi, snr = 4, seed = 1) {
  set.seed(seed)
  s <- simulate_guidebender2(n_guides = n_guides, n_cells = n_cells, moi = moi, hurdle_prob = 0.1,
        snr = snr, count_per_cell = 100, frac_noise_endo = 1, return_sparse_only = TRUE)
  list(counts = assay(s, "counts"), truth = assay(s, "ground_truth"))
}

## ===================== [3] BH vs GS vs BY realized FDR (sim) =====================
exp_bh_vs_gs <- function(N_REP = 12) {
  regs <- list(G20_moi0.3=c(20,.3,4), G20_moi1=c(20,1,4), G20_moi3=c(20,3,4), G50_moi5=c(50,5,4),
               fewcell_n300=c(20,3,4), lowSNR=c(20,1,2), fewcell_lowSNR=c(20,3,2), verydense=c(50,8,4))
  ncell <- c(2000,2000,2000,2000, 300,2000,300,1500)
  rows <- list()
  for (k in seq_along(regs)) { p <- regs[[k]]; rg <- names(regs)[k]
    for (rep in seq_len(N_REP)) {
      d <- sim1(p[1], ncell[k], p[2], p[3], seed = 7000*k+rep)
      for (cc in c("BH","GS","BY")) { e <- fdp(cassign(d$counts, q=Q, fdr=cc)$assigned, d$truth)
        rows[[length(rows)+1]] <- data.frame(regime=rg, guides_per_cell=round(mean(Matrix::colSums(d$truth>0)),2),
          rep=rep, corr=cc, fdp=e["fdp"], recall=e["recall"]) } }
    cat("  done", rg, "\n") }
  df <- do.call(rbind, rows); write.csv(df, file.path(OUT,"bh_vs_gs_fdr.csv"), row.names=FALSE)
  agg <- aggregate(cbind(fdp,recall)~regime+guides_per_cell+corr, df, mean)
  print(agg[order(agg$guides_per_cell, agg$corr),], row.names=FALSE, digits=3); invisible(df)
}

## ===================== [3] on-target power BH vs GS (DC-TAP, knockdown-validated) ==========
exp_poscontrol <- function() {
  run <- function(which) {
    DS <- dctap_sceptre_dir(which); gm <- load_grna_matrix(file.path(DS,"grna_matrix_aligned.rds"))
    rm_ <- as(readRDS(file.path(DS,"response_matrix.rds")),"CsparseMatrix")
    tdf <- read.csv(file.path(DS,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
    gene_of <- setNames(tdf$grna_target,tdf$grna_id); sym_of <- setNames(tdf$grna_target_symbol,tdf$grna_id)
    pos <- tdf$grna_id[tdf$grna_target_symbol!="non-targeting"]
    A_GS <- as(cassign(gm,q=Q,fdr="GS")$assigned,"lgCMatrix"); A_BH <- as(cassign(gm,q=Q,fdr="BH")$assigned,"lgCMatrix")
    sf <- Matrix::colSums(rm_)/median(Matrix::colSums(rm_)); ne <- function(G) log1p(as.numeric(rm_[G,])/sf)
    genes <- unique(unname(gene_of[pos]))
    ctrl <- lapply(setNames(genes,genes), function(G) as.numeric(sparseMatrixStats::colMaxs(gm[names(gene_of)[gene_of==G],,drop=FALSE]))<=2)
    kd <- function(tr,ct,e){ if(sum(tr)<3||sum(ct)<20) return(c(NA,NA)); tt<-tryCatch(t.test(e[tr],e[ct]),error=function(x)NULL)
      c(log2((mean(expm1(e[tr]))+1e-6)/(mean(expm1(e[ct]))+1e-6)), if(!is.null(tt))-log10(tt$p.value) else NA) }
    do.call(rbind, lapply(pos, function(g){ G<-gene_of[g]; e<-ne(G); tGS<-as.logical(A_GS[g,]); tBH<-as.logical(A_BH[g,])
      on<-tBH&!tGS; sGS<-kd(tGS,ctrl[[G]]&!tGS,e); sBH<-kd(tBH,ctrl[[G]]&!tBH,e); sON<-kd(on,ctrl[[G]]&!tBH,e)
      data.frame(dataset=which,guide=g,gene=sym_of[g],n_GS=sum(tGS),n_BH=sum(tBH),n_BHonly=sum(on),
        lfc_GS=sGS[1],lfc_BH=sBH[1],lfc_BHonly=sON[1],nlp_BHonly=sON[2]) }))
  }
  R <- do.call(rbind, lapply(c("highmoi","lowmoi"), run)); write.csv(R, file.path(OUT,"poscontrol_bh_vs_gs.csv"), row.names=FALSE)
  for (w in unique(R$dataset)) { d<-R[R$dataset==w,]; ok<-d[is.finite(d$lfc_BHonly),]
    cat(sprintf("[%s] pos-control cells GS=%d BH=%d (+%.1f%%) | median kd log2FC GS=%.3f BH=%.3f | BH-only knockdown %d/%d guides\n",
      w,sum(d$n_GS),sum(d$n_BH),100*(sum(d$n_BH)/sum(d$n_GS)-1),median(d$lfc_GS,na.rm=T),median(d$lfc_BH,na.rm=T),
      sum(ok$lfc_BHonly<0),nrow(ok))) }
  invisible(R)
}

## ===================== [4] which error UNIT: pooled-BH vs GS vs per-guide-BH (sim) =========
exp_perguide <- function(N_REP = 10) {
  metrics <- function(assigned, truth){ A<-as(assigned,"lgCMatrix"); Tr<-as(truth>0,"lgCMatrix"); tpm<-A&Tr
    na<-sum(A); tp<-sum(tpm); nt<-sum(Tr); na_g<-Matrix::rowSums(A); fp_g<-na_g-Matrix::rowSums(tpm)
    fp_c<-Matrix::colSums(A)-Matrix::colSums(tpm); pg<-ifelse(na_g>0,fp_g/na_g,NA)
    c(entry_fdr=if(na>0)(na-tp)/na else 0, perguide_fdr=mean(pg,na.rm=T), perguide_fdr_max=max(pg,na.rm=T),
      percell_fwer=mean(fp_c>0), recall=if(nt>0)tp/nt else NA) }
  regs <- list(G20_moi1=c(20,1), G20_moi3=c(20,3), G50_moi5=c(50,5)); rows <- list()
  for (k in seq_along(regs)) { p<-regs[[k]]
    for (rep in seq_len(N_REP)) { d<-sim1(p[1],2000,p[2],4,seed=7000*k+rep); fit<-cassign(d$counts,q=Q,fdr="GS")
      A<-list(poolBH=apply_corr(fit,"BH",Q), GS=apply_corr(fit,"GS",Q), perguideBH=apply_perguide(fit,Q))
      for (cc in names(A)) rows[[length(rows)+1]]<-data.frame(regime=names(regs)[k],rep=rep,corr=cc,t(metrics(A[[cc]],d$truth))) }
    cat("  done", names(regs)[k], "\n") }
  df <- do.call(rbind, rows); write.csv(df, file.path(OUT,"perguide_fdr.csv"), row.names=FALSE)
  print(aggregate(cbind(entry_fdr,perguide_fdr,perguide_fdr_max,percell_fwer,recall)~regime+corr, df, mean), row.names=FALSE, digits=3)
  invisible(df)
}

## ===================== [6] does assignment threshold move gRNA-level Type-I? (DC-TAP NTC) ===
exp_calibration <- function() {
  DS <- dctap_sceptre_dir("highmoi"); gm <- load_grna_matrix(file.path(DS,"grna_matrix_aligned.rds"))
  rm_ <- as(readRDS(file.path(DS,"response_matrix.rds")),"CsparseMatrix"); tdf <- read.csv(file.path(DS,"grna_target_data_frame.csv"),stringsAsFactors=FALSE)
  targ <- unique(tdf$grna_target[tdf$grna_target_symbol!="non-targeting"]); nt <- intersect(tdf$grna_id[tdf$grna_target_symbol=="non-targeting"], rownames(gm))
  E <- log1p(sweep(as.matrix(rm_),2,Matrix::colSums(rm_)/median(Matrix::colSums(rm_)),"/"))
  En <- E[setdiff(rownames(rm_),targ),,drop=FALSE]; Et <- E[targ,,drop=FALSE]; noG <- (gm<=2)
  rej <- function(M,tr,ct){ if(sum(tr)<3)return(NULL); Mt<-M[,tr,drop=F];Mc<-M[,ct,drop=F]
    mt<-rowMeans2(Mt);mc<-rowMeans2(Mc);vt<-rowVars(Mt);vc<-rowVars(Mc);se<-sqrt(vt/sum(tr)+vc/sum(ct));se[se==0]<-NA
    df<-(vt/sum(tr)+vc/sum(ct))^2/((vt/sum(tr))^2/(sum(tr)-1)+(vc/sum(ct))^2/(sum(ct)-1)); 2*pt(-abs((mt-mc)/se),df) }
  fit <- cassign(gm,q=Q,fdr="GS")
  A <- list(`GS(.05)`=apply_corr(fit,"GS",.05), `perguideBH(.05)`=apply_perguide(fit,.05),
            `BH(.05)`=apply_corr(fit,"BH",.05), `BH(.30)`=apply_corr(fit,"BH",.30))
  A <- lapply(A, function(m){rownames(m)<-rownames(gm); as(m,"lgCMatrix")})
  for (mn in names(A)) { Pn<-c(); Pt<-c()
    for (g in nt){ tr<-as.logical(A[[mn]][g,]); ct<-as.logical(noG[g,])&!tr; p<-rej(En,tr,ct); if(is.null(p))next
      Pn<-c(Pn,p); Pt<-c(Pt,rej(Et,tr,ct)) }
    cat(sprintf("%-16s clean-null frac<.05=%.3f | structure-sensor(target genes) frac<.05=%.3f\n",mn,mean(Pn<.05,na.rm=T),mean(Pt<.05,na.rm=T))) }
  # sanity bracket: random exchangeable vs depth-biased treatment group (identifies the confounder)
  set.seed(1); glib<-Matrix::colSums(gm); C<-ncol(gm); Pr<-c(); Pd<-c(); ord<-order(glib,decreasing=TRUE)
  for (k in 1:98){ tr<-logical(C); tr[sample(C,982)]<-TRUE; Pr<-c(Pr,rej(En,tr,!tr))
    td<-logical(C); td[ord[((k-1)*10+1):((k-1)*10+982)]]<-TRUE; Pd<-c(Pd,rej(En,td,!td)) }
  cat(sprintf("[sanity] RANDOM group frac<.05=%.3f (calibration floor) | DEPTH-biased group frac<.05=%.3f (confounder)\n",
      mean(Pr<.05,na.rm=T), mean(Pd<.05,na.rm=T)))
}

## ===================== [7] singleton QC: throw-out vs FDR (sim) =====================
exp_singleton_qc <- function(N_REP = 8) {
  analyze <- function(counts, truth, q){ fit<-cassign(counts,q=q,fdr="BH"); A<-as(fit$assigned,"lgCMatrix"); Tr<-as(truth>0,"lgCMatrix")
    na<-Matrix::colSums(A); nt<-Matrix::colSums(Tr); ntp<-Matrix::colSums(A&Tr); nfp<-na-ntp
    ass<-na>=1; thr<-na>=2; kept<-na==1
    c(fdr=sum(nfp)/max(sum(na),1), throwout=mean(thr[ass]), frac_FP=if(sum(thr)>0)mean((nt<2)[thr])else NA,
      collision=sum(nfp[nt>=1])/max(sum(nfp),1), goodloss=sum(thr&nt==1&ntp>=1)/max(sum(nt==1),1),
      false_singleton=sum(kept&ntp==0&nfp==1)/max(sum(kept),1)) }
  for (moi in c(0.3,0.6)) { cat(sprintf("--- MOI %.1f ---\n",moi))
    for (q in c(0.05,0.10,0.20,0.35)) { rows<-t(sapply(1:N_REP, function(r){ d<-sim1(30,3000,moi,4,seed=round(1000*moi)+100*r+round(100*q)); analyze(d$counts,d$truth,q) }))
      m<-colMeans(rows,na.rm=TRUE)
      cat(sprintf("  q=%.2f realFDR=%.3f throwout=%.1f%% (FP-share=%.0f%%, collision=%.2f) goodCellLoss=%.1f%% falseSingleton=%.1f%%\n",
          q,m["fdr"],100*m["throwout"],100*m["frac_FP"],m["collision"],100*m["goodloss"],100*m["false_singleton"])) } }
}

## ---- dispatch ----
if (length(commandArgs(TRUE)) > 0) {
  what <- commandArgs(TRUE)[1]
  if (what %in% c("all","bh_vs_gs"))    { cat("\n== [3] BH vs GS vs BY realized FDR ==\n");   exp_bh_vs_gs() }
  if (what %in% c("all","poscontrol"))  { cat("\n== [3] on-target power (DC-TAP) ==\n");       exp_poscontrol() }
  if (what %in% c("all","perguide"))    { cat("\n== [4] error unit: per-guide FDR ==\n");      exp_perguide() }
  if (what %in% c("all","calibration")) { cat("\n== [6] gRNA-level Type-I calibration ==\n");   exp_calibration() }
  if (what %in% c("all","singleton_qc")){ cat("\n== [7] singleton QC vs FDR ==\n");            exp_singleton_qc() }
} else cat("defines exp_*() functions. Run: Rscript scripts/error_control_experiments.R {all|bh_vs_gs|poscontrol|perguide|calibration|singleton_qc}\n")
