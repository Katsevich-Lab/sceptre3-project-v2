## ============================================================================
## Collaborator writeup §3: reproduce the RAMAC two-figure diagnostic in the two
## OTHER cell types (EndoC + CD4 T-cell). For each dataset we pick one clean-gap
## guide, run the rank-1 DENOISED depth-mixed Poisson fit, and produce:
##   (a) gap_hist_endoc_cd4.png  — the full count distribution (low ambient mode,
##       empty gap, high signal mode), mirroring gap_hist_RAMAC.png.
##   (b) ambient_fit_endoc_cd4.png — the below-gap ambient counts vs the depth-mixed
##       Poisson a_g*d_c (orange line/points, log1p y-axis), mirroring
##       ramac_ambient_fit.png; count-2 observed-vs-predicted annotated per panel.
##
## The fit + count-2 metrics mirror scripts/global_ambient_gap.R EXACTLY (same
## fit_rank1_denoised + detect_gap, same a_g-from-ambient + shared-depth-histogram
## expected PMF), so obs2/exp2_single/exp2_den reproduce the perguide CSVs.
##
## Run from nonparametric-thresholds/.  Packages: Matrix, ggplot2.
## ============================================================================
suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
OUT <- "results/collaborator_writeup"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
MODEL_COL <- "#D55E00"

## ---- COPIED verbatim from scripts/global_ambient_gap.R (do NOT source it) ----
fit_rank1_denoised <- function(mc, outer=6, inner=12, mask_p=1e-6){
  G <- nrow(mc); C <- ncol(mc)
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
  list(a=a, d=d, gv=gv, cv=cv, xv=xv, colN=colN)
}
GAP_RATIO <- 3; GAP_ABS <- 3; MIN_AMB <- 50; MIN_SIG <- 10; MAX_LO <- 40
detect_gap <- function(vg){
  if(length(vg) < MIN_AMB + MIN_SIG) return(NULL)
  occ <- sort(unique(vg))
  if(length(occ) < 2) return(NULL)
  for(i in seq_len(length(occ)-1)){
    lo <- occ[i]; hi <- occ[i+1]
    if(lo > MAX_LO) return(NULL)
    if(hi/lo >= GAP_RATIO && hi-lo >= GAP_ABS){
      n_amb <- sum(vg <= lo); n_sig <- sum(vg >= hi)
      if(n_amb >= MIN_AMB && n_sig >= MIN_SIG) return(list(lo=lo, hi=hi))
    }
  }
  NULL
}

## ---- datasets + hand-picked clean-gap guides (verified vs perguide CSVs) ------
## Guides chosen for a WIDE, CONTIGUOUS low mode (occupied at consecutive counts 1..gap_lo
## with no internal empty stretch) that actually hosts weak integrations, so the below-gap
## excess is visible across the whole low mode (not crammed into count 2):
##   EndoC KCNJ11-2 (gap 4, contiguous [1,2,3,4], modest excess — dense soup);
##   CD4 PIP-2 (gap 5, large excess — thin soup).
## (WFS1-3, the prior EndoC pick, was rejected: its low mode [1,2,3,8] is NON-contiguous —
##  a stray cell at 8 above an empty 4-7 stretch, so "gap 8" is misleading.)
DATASETS <- list(
  list(key="endoc", label="EndoC (T2D perturb-seq)", guide="KCNJ11-2",
       path="~/data/external/perturbseq-survey/endoc_t2d_perturbseq_GSE273677/grna_matrix.rds"),
  list(key="cd4",   label="CD4 T-cell perturb-seq",   guide="PIP-2",
       path="~/data/external/perturbseq-survey/tcell_cd4_perturbseq_GSE314342/grna_matrix.rds")
)

## For one dataset: load, fit, isolate the chosen guide's ambient sample, and
## return everything both figures need + the count-2 verification numbers.
analyze_one <- function(ds){
  cat(sprintf("\n[%s] loading %s\n", ds$key, ds$path))
  mc <- as(readRDS(path.expand(ds$path)), "CsparseMatrix")
  if(nrow(mc) > ncol(mc)) mc <- as(t(mc), "CsparseMatrix")   # ensure guides x cells
  G <- nrow(mc); C <- ncol(mc)
  cat(sprintf("[%s] %d guides x %d cells\n", ds$key, G, C))

  fit <- fit_rank1_denoised(mc)
  d <- fit$d; dsum <- sum(d)
  ## shared depth histogram for the denoised-depth-mixed Poisson expected PMF
  dpos <- d[d>0]; ld <- log(dpos)
  hd <- hist(ld, breaks=250, plot=FALSE); dbin <- exp(hd$mids); wbin <- hd$counts; Wtot <- sum(wbin)

  g <- match(ds$guide, rownames(mc))
  if(is.na(g)) stop(sprintf("[%s] guide %s not found", ds$key, ds$guide))
  pos <- which(fit$gv == g)                 # nonzero positions for this guide
  vg  <- fit$xv[pos]                          # this guide's positive counts
  gap <- detect_gap(vg)
  if(is.null(gap)) stop(sprintf("[%s] guide %s has NO clean gap under detect_gap()", ds$key, ds$guide))
  lo <- gap$lo; hi <- gap$hi
  is_amb <- vg <= lo; is_sig <- vg >= hi
  n_sig <- sum(is_sig); n_amb_tot <- C - n_sig   # ambient cells incl. zeros (gap is empty)
  ks <- 0:lo

  ## observed ambient histogram over 0..lo (obs[0] = ambient cells minus positive-ambient)
  tab_pos <- tabulate(vg[is_amb], nbins=lo); obs <- c(n_amb_tot - sum(tab_pos), tab_pos)
  ## a_g from ambient cells only (exclude the signal cells' depths & counts)
  d_sig <- sum(d[fit$cv[pos[is_sig]]]); amb_count_sum <- sum(vg[is_amb])
  a_g <- amb_count_sum / max(dsum - d_sig, 1e-9)
  ## expected PMFs over the n_amb_tot ambient cells
  lam1 <- amb_count_sum / n_amb_tot
  exp_single <- dpois(ks, lam1) * n_amb_tot
  exp_den <- sapply(ks, function(k) sum(wbin * dpois(k, a_g*dbin)) / Wtot) * n_amb_tot

  ## full count distribution for the gap histogram (all C cells)
  vfull <- numeric(C); vfull[fit$cv[pos]] <- vg     # 0 elsewhere
  list(ds=ds, lo=lo, hi=hi, n_amb=n_amb_tot, n_amb_pos=sum(is_amb), n_sig=n_sig,
       a_g=a_g, ks=ks, obs=obs, exp_single=exp_single, exp_den=exp_den,
       vg=vg, vfull=vfull, maxk=max(vg))
}

res <- lapply(DATASETS, analyze_one)

## ---- HONEST VERIFICATION: count-2 obs vs pred (single + denoised) vs CSV ------
csv_endoc <- read.csv("results/global_ambient_poisson/perguide_endoc_t2d.csv", stringsAsFactors=FALSE)
csv_cd4   <- read.csv("results/global_ambient_poisson/perguide_tcell_cd4.csv", stringsAsFactors=FALSE)
csv_lookup <- function(key, guide){
  df <- if(key=="endoc") csv_endoc else csv_cd4
  df[df$guide==guide, c("obs2","exp2_single","exp2_den","c2_oe_single","c2_oe_den"), drop=FALSE]
}
cat("\n================ COUNT-2 VERIFICATION (computed vs perguide CSV) ================\n")
for(r in res){
  ci <- match(2, r$ks)                       # index of count==2
  csv <- csv_lookup(r$ds$key, r$ds$guide)
  cat(sprintf("\n[%s] guide %s | gap %d -> %d | n_amb=%d (pos=%d) n_sig=%d\n",
      r$ds$key, r$ds$guide, r$lo, r$hi, r$n_amb, r$n_amb_pos, r$n_sig))
  cat(sprintf("   count-2 computed: obs=%d  exp_single=%.2f (%.2fx)  exp_den=%.2f (%.2fx)\n",
      r$obs[ci], r$exp_single[ci], r$obs[ci]/r$exp_single[ci], r$exp_den[ci], r$obs[ci]/r$exp_den[ci]))
  cat(sprintf("   count-2 CSV     : obs=%d  exp_single=%.2f (%.2fx)  exp_den=%.2f (%.2fx)\n",
      csv$obs2, csv$exp2_single, csv$c2_oe_single, csv$exp2_den, csv$c2_oe_den))
  ok <- (r$obs[ci]==csv$obs2) &&
        isTRUE(all.equal(round(r$exp_single[ci],2), csv$exp2_single, tolerance=0.02)) &&
        isTRUE(all.equal(round(r$exp_den[ci],2),    csv$exp2_den,    tolerance=0.02))
  cat(sprintf("   MATCH: %s\n", if(ok) "YES" else "*** NO — investigate ***"))
}

## ============================================================================
## FIGURE (a): gap histogram — two panels, one guide each (full distribution)
## Log-log density mirroring gap_hist_RAMAC.png, but as a compact ggplot facet.
## ============================================================================
FLOOR <- 3e-8
be <- function(maxc, thr=12, nlog=32){
  b <- c(0:thr, round(exp(seq(log(thr+1), log(maxc+1), length.out=nlog))))
  b <- sort(unique(b[b>=0 & b<=maxc])); c(b, maxc+1)
}
gap_df <- do.call(rbind, lapply(res, function(r){
  v <- r$vfull; C <- length(v); mxk <- r$maxk
  tabk <- tabulate(v+1, nbins=mxk+1); E <- be(mxk); mb <- length(E)-1
  lo <- E[1:mb]; hi <- E[2:(mb+1)]-1; wid <- hi-lo+1
  dens <- sapply(1:mb, function(i) sum(tabk[(lo[i]+1):(hi[i]+1)]))/C/wid
  data.frame(panel=sprintf("%s\nguide %s (gap %d–%d)", r$ds$label, r$ds$guide, r$lo, r$hi),
             xl=lo+0.5+1, xr=hi+1.5, ymin=FLOOR, ymax=pmax(dens, FLOOR),
             stringsAsFactors=FALSE)
}))
## empty-gap shading rectangles per panel
gap_rects <- do.call(rbind, lapply(res, function(r){
  data.frame(panel=sprintf("%s\nguide %s (gap %d–%d)", r$ds$label, r$ds$guide, r$lo, r$hi),
             xmin=r$lo+1.5, xmax=r$hi+0.5, ymin=FLOOR, ymax=2,
             lab=sprintf("empty gap\n(no cells\ncount %d–%d)", r$lo+1, r$hi-1),
             lx=sqrt((r$lo+1.5)*(r$hi+0.5)), ly=0.18, stringsAsFactors=FALSE)
}))
panel_levels <- unique(gap_df$panel)
gap_df$panel   <- factor(gap_df$panel,   levels=panel_levels)
gap_rects$panel <- factor(gap_rects$panel, levels=panel_levels)

pa <- ggplot() +
  geom_rect(data=gap_rects, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill="#D55E00", alpha=0.10) +
  geom_rect(data=gap_df, aes(xmin=xl, xmax=xr, ymin=ymin, ymax=ymax),
            fill="grey82", colour="grey55", linewidth=0.25) +
  geom_text(data=gap_rects, aes(x=lx, y=ly, label=lab), colour="darkorange3", size=4.3) +
  scale_x_log10() +
  scale_y_log10(labels=scales::label_number(drop0trailing=TRUE)) +
  facet_wrap(~panel, nrow=1, scales="free_x") +
  labs(x="gRNA UMI count + 1", y="density   P(count) / count") +
  theme_bw(base_size=16) +
  theme(strip.text=element_text(face="bold", size=14))
ggsave(file.path(OUT, "gap_hist_endoc_cd4.png"), pa, width=10, height=4.6, dpi=130)
cat("\nwrote gap_hist_endoc_cd4.png\n")

## ============================================================================
## FIGURE (b): ambient-mode fit — two panels, mirroring ramac_ambient_fit.png.
## Below-gap ambient counts (grey bars) vs depth-mixed Poisson (orange), log1p y.
## ============================================================================
## show the FULL low mode (counts 0..gap_lo) so the weak-integration excess is
## visible across the whole mode, not pooled into a tail bin.
LMAX <- max(vapply(res, function(r) r$lo, numeric(1)))
fit_df <- do.call(rbind, lapply(res, function(r){
  ks <- 0:r$lo
  data.frame(panel=sprintf("%s — guide %s", r$ds$label, r$ds$guide),
             klab=as.character(ks),
             observed=r$obs[match(ks, r$ks)], model=r$exp_den[match(ks, r$ks)],
             stringsAsFactors=FALSE)
}))
fit_df$klab <- factor(fit_df$klab, levels=as.character(0:LMAX))
## weak-integration annotation per panel: TOTAL excess over counts >= 2 (the whole
## low-mode region above pure ambient), placed top-right where the bars are short.
anno_df <- do.call(rbind, lapply(res, function(r){
  reg <- r$ks >= 2
  obs_r <- sum(r$obs[reg]); pred_r <- sum(r$exp_den[reg])
  data.frame(panel=sprintf("%s — guide %s", r$ds$label, r$ds$guide),
             klab=factor(as.character(r$lo), levels=as.character(0:LMAX)),
             ty=max(r$obs)*0.7,
             lab=sprintf("counts 2–%d: %d obs vs ~%.0f pred (%.1f×)", r$lo, obs_r, pred_r, obs_r/pred_r),
             stringsAsFactors=FALSE)
}))
plvl <- unique(fit_df$panel)
fit_df$panel  <- factor(fit_df$panel,  levels=plvl)
anno_df$panel <- factor(anno_df$panel, levels=plvl)

pb <- ggplot(fit_df, aes(klab, observed)) +
  geom_col(fill="grey80", width=0.66) +
  geom_line(aes(y=model, group=1), colour=MODEL_COL, linewidth=0.9) +
  geom_point(aes(y=model), colour=MODEL_COL, size=2.6) +
  geom_text(data=anno_df, aes(x=klab, y=ty, label=lab), inherit.aes=FALSE,
            hjust=1, vjust=1, size=4.3, colour="grey20") +
  facet_wrap(~panel, nrow=1, scales="free") +
  scale_y_continuous(trans=scales::trans_new("log1p", log1p, expm1),
                     breaks=c(0,1,10,100,1000,10000),
                     labels=scales::label_comma(), expand=expansion(mult=c(0,0.10))) +
  labs(x="gRNA UMI count", y="number of cells (log1p — axis floor at 0)") +
  theme_bw(base_size=16) +
  theme(strip.text=element_text(face="bold", size=14))
ggsave(file.path(OUT, "ambient_fit_endoc_cd4.png"), pb, width=10, height=4.6, dpi=130)
cat("wrote ambient_fit_endoc_cd4.png\n")
