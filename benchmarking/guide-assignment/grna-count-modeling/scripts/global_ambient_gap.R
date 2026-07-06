## ============================================================================
## GLOBAL clean-gap ambient-Poisson analysis.
## Generalizes the RAMAC case study (replogle_ambient_poisson.qmd) from ONE guide
## to MANY guides across MANY datasets. For each guide with a clean empty gap
## between its ambient and signal modes, isolate the below-gap ambient sample and
## ask: is it Poisson? Single Poisson vs fishash-style rank-1 DENOISED-depth-mixed
## Poisson. The headline per-guide metric mirrors RAMAC's count-2 excess:
##   c2_oe = observed(count=2) / expected(count=2)  under each null.
## RAMAC: single ~53/6=8.8x, denoised ~53/9=5.9x. Does this generalize?
##
## Usage:  Rscript scripts/global_ambient_gap.R <dataset_name> <matrix_rds_path>
## Writes: results/global_ambient_poisson/perguide_<dataset_name>.csv
## Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix)})

args <- commandArgs(trailingOnly=TRUE)
NAME <- args[1]; PATH <- path.expand(args[2])
OUT  <- "results/global_ambient_poisson"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

## ---- rank-1 Poisson denoised ambient depth + clean-gap detector ----
## detect_gap: first sufficiently-wide empty band after the ambient mode (ratio>=3, width>=3,
## >=50 ambient & >=10 signal cells, ambient top <=40). Both from the shared lib:
source("scripts/ambient_lib.R")   # fit_rank1_denoised(), detect_gap()

cat(sprintf("[%s] loading %s\n", NAME, PATH))
mc <- as(readRDS(PATH), "CsparseMatrix")
if(nrow(mc) > ncol(mc)) mc <- as(t(mc), "CsparseMatrix")   # ensure guides x cells
G <- nrow(mc); C <- ncol(mc)
cat(sprintf("[%s] %d guides x %d cells, %d nnz\n", NAME, G, C, length(mc@x)))

fit <- fit_rank1_denoised(mc)
d <- fit$d; dsum <- sum(d)
## depth histogram for the denoised-depth-mixed Poisson expected PMF (shared)
dpos <- d[d>0]; ld <- log(dpos)
hd <- hist(ld, breaks=250, plot=FALSE); dbin <- exp(hd$mids); wbin <- hd$counts
Wtot <- sum(wbin)

gv <- fit$gv; xv <- fit$xv
pos_by_g <- split(seq_along(gv), gv)        # nonzero positions per guide

rows <- vector("list", G); kk <- 0L
for(g in seq_len(G)){
  pos <- pos_by_g[[as.character(g)]]
  if(is.null(pos)) next
  vg <- xv[pos]                              # this guide's positive counts
  gap <- detect_gap(vg)
  if(is.null(gap)) next
  lo <- gap$lo; hi <- gap$hi
  is_amb <- vg <= lo                         # positive ambient counts (1..lo)
  is_sig <- vg >= hi
  n_sig <- sum(is_sig)
  n_amb_tot <- C - n_sig                     # ambient cells incl. zeros (gap is empty)
  ks <- 0:lo
  ## observed ambient histogram (0..lo); obs[0] = ambient cells minus positive-ambient
  tab_pos <- tabulate(vg[is_amb], nbins=lo)
  obs <- c(n_amb_tot - sum(tab_pos), tab_pos)
  ## a_g from ambient cells only (exclude the guide's signal cells' depths & counts)
  d_sig <- sum(d[fit$cv[pos[is_sig]]])
  amb_count_sum <- sum(vg[is_amb])
  a_g <- amb_count_sum / max(dsum - d_sig, 1e-9)
  ## expected PMFs over the n_amb_tot ambient cells
  lam1 <- amb_count_sum / n_amb_tot          # single Poisson
  exp_single <- dpois(ks, lam1) * n_amb_tot
  ## denoised-depth-mixed Poisson via shared depth histogram
  exp_den <- sapply(ks, function(k) sum(wbin * dpois(k, a_g*dbin)) / Wtot) * n_amb_tot
  ## metrics (guard tiny expecteds)
  c2_oe_single <- if(lo>=2) obs[3]/max(exp_single[3],1e-9) else NA
  c2_oe_den    <- if(lo>=2) obs[3]/max(exp_den[3],1e-9)   else NA
  tail_obs <- sum(obs[ks>=2]); tail_exp_den <- sum(exp_den[ks>=2]); tail_exp_sing <- sum(exp_single[ks>=2])
  tail_oe_den  <- if(lo>=2) tail_obs/max(tail_exp_den,1e-9) else NA
  tail_oe_sing <- if(lo>=2) tail_obs/max(tail_exp_sing,1e-9) else NA
  maxrate <- a_g * max(d)
  kk <- kk + 1L
  rows[[kk]] <- data.frame(dataset=NAME, guide=rownames(mc)[g],
    gap_lo=lo, gap_hi=hi, gap_ratio=round(hi/lo,2),
    n_amb=sum(is_amb), n_sig=n_sig, a_g=a_g, max_amb_rate=round(maxrate,4),
    obs2=if(lo>=2) obs[3] else NA, exp2_single=if(lo>=2) round(exp_single[3],2) else NA,
    exp2_den=if(lo>=2) round(exp_den[3],2) else NA,
    c2_oe_single=round(c2_oe_single,2), c2_oe_den=round(c2_oe_den,2),
    tail_oe_single=round(tail_oe_sing,2), tail_oe_den=round(tail_oe_den,2),
    stringsAsFactors=FALSE)
}
res <- if(kk>0) do.call(rbind, rows[seq_len(kk)]) else data.frame()
cat(sprintf("[%s] clean-gap guides: %d of %d (%.1f%%)\n", NAME, kk, G, 100*kk/G))
if(kk>0){
  cat(sprintf("[%s] median c2_oe: single=%.2f  denoised=%.2f | median tail_oe denoised=%.2f\n",
      NAME, median(res$c2_oe_single,na.rm=TRUE), median(res$c2_oe_den,na.rm=TRUE),
      median(res$tail_oe_den,na.rm=TRUE)))
  write.csv(res, file.path(OUT, sprintf("perguide_%s.csv", NAME)), row.names=FALSE)
  cat(sprintf("[%s] wrote %s\n", NAME, file.path(OUT, sprintf("perguide_%s.csv", NAME))))
}
