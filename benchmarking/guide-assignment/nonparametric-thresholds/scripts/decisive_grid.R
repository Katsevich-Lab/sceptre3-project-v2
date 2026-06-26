#!/usr/bin/env Rscript
# Decisive regime map: cross the two failure-mode axes (Simpson via endo_shape_flat,
# overdispersion via Phi_noise) plus a hard low-signal regime, and measure realized FDR
# and recall for the plain ambient test, our NB candidate, and fishash. Authors' simulator.
suppressPackageStartupMessages({library(fishash); library(Matrix); library(SummarizedExperiment)})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))),error=function(e) ".")
HERE <- dirname(HERE); source(file.path(HERE,"scripts","ambient_nb_method.R"))
RES <- file.path(HERE,"results","benchmark_update"); dir.create(RES,showWarnings=FALSE,recursive=TRUE)

score <- function(A, truth){A<-as(A,"lgCMatrix");truth<-as(truth,"lgCMatrix")
  na<-sum(A);nt<-sum(truth);tp<-sum(A&truth)
  c(fdr=if(na>0)(na-tp)/na else 0, recall=if(nt>0)tp/nt else NA_real_)}

regimes <- list(
  clean        = list(flat=0,   phi=0, snr=4, cpc=100),  # no Simpson, equidispersed, strong signal
  simpson      = list(flat=1,   phi=0, snr=4, cpc=100),  # Simpson only
  overdisp     = list(flat=0,   phi=1, snr=4, cpc=100),  # overdispersion only
  both         = list(flat=1,   phi=1, snr=4, cpc=100),  # Simpson + overdispersion
  hard_lowsig  = list(flat=0.5, phi=0, snr=1, cpc=20))   # low-signal hard regime (fishash's edge)
NREP <- 5; Q <- 0.05
rows <- list()
for (rg in names(regimes)) {
  p <- regimes[[rg]]
  for (rep in seq_len(NREP)) {
    set.seed(100*which(names(regimes)==rg)+rep)
    sim <- simulate_guidebender2(n_guides=200, n_cells=4000, moi=0.3, hurdle_prob=0.1,
             snr=p$snr, count_per_cell=p$cpc, frac_noise_endo=0.75,
             endo_shape_flat=p$flat, Phi_noise=p$phi, return_sparse_only=TRUE)
    N <- assay(sim,"counts"); truth <- assay(sim,"ground_truth")
    amb <- score(ambient_nb_assign(N,q=Q,n_iter=1,overdispersion=FALSE)$assigned, truth)
    nb  <- score(ambient_nb_assign(N,q=Q,n_iter=3,overdispersion=TRUE)$assigned,  truth)
    fsh <- score(assay(fishash(N,padj_cutoff=Q),"assigned"), truth)
    rows[[length(rows)+1]] <- data.frame(regime=rg, rep=rep,
      amb_fdr=amb[1],amb_rec=amb[2], nb_fdr=nb[1],nb_rec=nb[2], fsh_fdr=fsh[1],fsh_rec=fsh[2])
    cat(sprintf("  %-12s rep%d | ambient FDR %.3f rec %.2f | NB FDR %.3f rec %.2f | fishash FDR %.3f rec %.2f\n",
        rg,rep, amb[1],amb[2], nb[1],nb[2], fsh[1],fsh[2]))
  }
}
df <- do.call(rbind, rows); write.csv(df, file.path(RES,"decisive_grid.csv"), row.names=FALSE)
cat("\n================ SUMMARY (mean over",NREP,"reps; #fail = reps with FDR>",Q,") ================\n")
cat(sprintf("%-12s | %-22s | %-22s | %-22s\n","regime","ambient (Poisson)","NB candidate (ours)","fishash"))
cat(sprintf("%-12s | %-22s | %-22s | %-22s\n","","FDR   recall  #fail","FDR   recall  #fail","FDR   recall  #fail"))
for (rg in names(regimes)) { s<-df[df$regime==rg,]
  cat(sprintf("%-12s | %.3f  %.2f    %d/%d | %.3f  %.2f    %d/%d | %.3f  %.2f    %d/%d\n", rg,
    mean(s$amb_fdr),mean(s$amb_rec),sum(s$amb_fdr>Q),NREP,
    mean(s$nb_fdr), mean(s$nb_rec), sum(s$nb_fdr>Q), NREP,
    mean(s$fsh_fdr),mean(s$fsh_rec),sum(s$fsh_fdr>Q),NREP)) }
cat("\nwrote", file.path(RES,"decisive_grid.csv"),"\n")
