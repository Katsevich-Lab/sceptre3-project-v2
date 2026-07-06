# ============================================================================
# scripts/ambient_lib.R -- shared ambient-model helpers (source-safe: functions only).
#
# fit_rank1_denoised(counts): rank-1 Poisson ambient field a_g * d_c fit by masked EM
# (mask out entries whose upper-tail Poisson p < mask_p each outer iter, refit). Returns
# list(a, d, mask, gv, cv, xv, colN). Promoted from ambient_validation.R; also used by
# global_ambient_gap.R and weak_integration_prevalence.R (which ignore $mask).
# ============================================================================
fit_rank1_denoised <- function(counts, outer=6, inner=12, mask_p=1e-6){
  mc <- as(counts, "CsparseMatrix"); G <- nrow(mc); C <- ncol(mc)
  gv <- mc@i + 1L; cv <- rep.int(seq_len(C), diff(mc@p)); xv <- mc@x
  rowN <- as.numeric(Matrix::rowSums(mc)); colN <- as.numeric(Matrix::colSums(mc))
  fastsum <- function(val, idx, n){ o<-numeric(n); s<-rowsum(val, idx); o[as.integer(rownames(s))]<-s[,1]; o }
  d <- colN; d[d==0] <- 1e-6; a <- rowN/sum(rowN); mask <- logical(length(xv))
  for (o in seq_len(outer)){
    for (it in seq_len(inner)){
      if(any(mask)){ mrN<-fastsum(xv[mask],gv[mask],G); mrD<-fastsum(d[cv[mask]],gv[mask],G) } else { mrN<-numeric(G); mrD<-numeric(G) }
      a <- (rowN-mrN)/(sum(d)-mrD); a[!is.finite(a)|a<0]<-0; a<-a/sum(a)
      if(any(mask)){ mcN<-fastsum(xv[mask],cv[mask],C); mcA<-fastsum(a[gv[mask]],cv[mask],C) } else { mcN<-numeric(C); mcA<-numeric(C) }
      d <- (colN-mcN)/(1-mcA); d[!is.finite(d)|d<0]<-0
    }
    mu <- a[gv]*d[cv]; mask <- ppois(xv-1, mu, lower.tail=FALSE) < mask_p
  }
  list(a=a, d=d, mask=mask, gv=gv, cv=cv, xv=xv, colN=colN)
}

# detect_gap(vg): find a low ambient mode / high signal mode separated by an empty gap in a guide's
# positive counts vg. Params default to the values previously hard-coded as globals in
# global_ambient_gap.R / weak_integration_prevalence.R (GAP_RATIO=3, GAP_ABS=3, MIN_AMB=50,
# MIN_SIG=10, MAX_LO=40). Returns list(lo, hi) or NULL.
detect_gap <- function(vg, min_amb=50, min_sig=10, max_lo=40, gap_ratio=3, gap_abs=3){
  if(length(vg) < min_amb + min_sig) return(NULL)
  occ <- sort(unique(vg))
  if(length(occ) < 2) return(NULL)
  for(i in seq_len(length(occ)-1)){
    lo <- occ[i]; hi <- occ[i+1]
    if(lo > max_lo) return(NULL)            # ambient mode must be a low mode
    if(hi/lo >= gap_ratio && hi-lo >= gap_abs){
      n_amb <- sum(vg <= lo); n_sig <- sum(vg >= hi)
      if(n_amb >= min_amb && n_sig >= min_sig) return(list(lo=lo, hi=hi))
    }
  }
  NULL
}
