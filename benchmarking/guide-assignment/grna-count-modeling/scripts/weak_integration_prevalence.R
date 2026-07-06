## ============================================================================
## WEAK-INTEGRATION PREVALENCE across the three clean-gap datasets.
##
## Question: among WIDE-GAP guides (low mode wide enough to host weak
## integrations), how many actually SUFFER from the weak-integration issue --
## i.e. a genuine, statistically-significant EXCESS of below-gap cells above the
## ambient depth-mixed Poisson in the counts-2..gap_lo region (real weak
## integrations hiding below the gap), versus a low mode fully explained by
## ambient soup? This measures how prevalent weak integration is, and should
## track soup density (dense soup -> few guides affected; thin soup -> many).
##
## fit_rank1_denoised() and detect_gap() now come from the shared scripts/ambient_lib.R
## (do NOT source global_ambient_gap.R -- it runs top-level CLI code).
##
## Read-only. Writes ONLY new files under results/collaborator_writeup/.
##   Run from grna-count-modeling/.
## ============================================================================
suppressMessages({library(Matrix)})

OUT <- "results/collaborator_writeup"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

## ---- rank-1 Poisson denoised ambient depth + clean-gap detector (shared lib) ----
source("scripts/ambient_lib.R")   # fit_rank1_denoised(), detect_gap()  (defaults match the old globals)
source("scripts/datasets.R")      # dataset_paths()

## ============================================================================
## Per-dataset analysis: clean-gap detection + weak-integration excess test.
## For each clean-gap guide we compute, over the AMBIENT (below-gap) cells and
## the counts-2..gap_lo weak-integration region:
##   - observed cell count in region
##   - expected cell count under the DENOISED depth-mixed Poisson (a_g*d_c)
##   - p = upper-tail Poisson prob of the excess: ppois(obs-1, exp, lower=FALSE)
## Effect size = obs/exp in region. BH within dataset; suffers = p_adj < 0.05.
## We also record obs2/exp2_den for cross-check against the pipeline perguide CSV.
## ============================================================================
analyze_dataset <- function(NAME, PATH){
  cat(sprintf("\n[%s] loading %s\n", NAME, PATH))
  mc <- tryCatch(as(readRDS(path.expand(PATH)), "CsparseMatrix"),
                 error=function(e){ cat(sprintf("[%s] LOAD FAILED: %s\n", NAME, conditionMessage(e))); NULL })
  if(is.null(mc)) return(NULL)
  if(nrow(mc) > ncol(mc)) mc <- as(t(mc), "CsparseMatrix")   # ensure guides x cells
  G <- nrow(mc); C <- ncol(mc)
  cat(sprintf("[%s] %d guides x %d cells, %d nnz\n", NAME, G, C, length(mc@x)))

  fit <- tryCatch(fit_rank1_denoised(mc),
                  error=function(e){ cat(sprintf("[%s] FIT FAILED: %s\n", NAME, conditionMessage(e))); NULL })
  if(is.null(fit)) return(NULL)
  d <- fit$d; dsum <- sum(d)
  if(!all(is.finite(d)) || dsum <= 0){ cat(sprintf("[%s] FIT non-finite depths\n", NAME)); return(NULL) }

  gv <- fit$gv; xv <- fit$xv; cv <- fit$cv
  pos_by_g <- split(seq_along(gv), gv)

  rows <- vector("list", G); kk <- 0L
  for(g in seq_len(G)){
    pos <- pos_by_g[[as.character(g)]]
    if(is.null(pos)) next
    vg <- xv[pos]                              # this guide's positive counts
    gap <- detect_gap(vg)
    if(is.null(gap)) next
    lo <- gap$lo; hi <- gap$hi
    is_amb <- vg <= lo
    is_sig <- vg >= hi
    n_sig <- sum(is_sig)
    n_amb_tot <- C - n_sig                     # ambient cells incl. zeros (gap is empty)

    ## a_g from ambient cells only (matches global_ambient_gap.R exactly)
    d_sig <- sum(d[cv[pos[is_sig]]])
    amb_count_sum <- sum(vg[is_amb])
    a_g <- amb_count_sum / max(dsum - d_sig, 1e-9)

    ## expected ambient histogram (0..lo) under denoised depth-mixed Poisson.
    ## Expectation over the ACTUAL ambient cells' depths (per-cell dpois), which
    ## for the region-sum is exact. Ambient cells = all cells except the guide's
    ## signal cells. Build per-count expected via depth-histogram (fast, matches
    ## the pipeline's shared-histogram approach).
    ks <- 0:lo
    ## depth vector over ambient cells: all cells minus signal cells of this guide
    d_all <- d
    ## expected using per-cell depths of ambient cells only:
    sig_cells <- cv[pos[is_sig]]
    is_amb_cell <- rep(TRUE, C); is_amb_cell[sig_cells] <- FALSE
    d_amb <- d_all[is_amb_cell]
    d_amb <- d_amb[d_amb > 0]
    ## exact per-cell expected counts in region 2..lo (denoised):
    exp_region_den <- 0
    obs_region <- sum(vg[is_amb] >= 2 & vg[is_amb] <= lo)
    if(lo >= 2){
      for(k in 2:lo) exp_region_den <- exp_region_den + sum(dpois(k, a_g * d_amb))
      ## exp at exactly 2 for cross-check
      exp2_den <- sum(dpois(2, a_g * d_amb))
      obs2 <- sum(vg[is_amb] == 2)
    } else {
      exp2_den <- NA; obs2 <- NA
    }

    oe_region <- if(lo>=2 && exp_region_den>0) obs_region/exp_region_den else NA
    ## upper-tail Poisson p-value for the excess of cells in region
    p <- if(lo>=2) ppois(obs_region - 1, exp_region_den, lower.tail=FALSE) else NA

    kk <- kk + 1L
    rows[[kk]] <- data.frame(dataset=NAME, guide=rownames(mc)[g],
      gap_lo=lo, gap_hi=hi, n_amb=sum(is_amb), a_g=a_g,
      obs_region=obs_region, exp_region_den=round(exp_region_den,3),
      oe_region=round(oe_region,3), p=p,
      obs2=obs2, exp2_den=round(exp2_den,3),
      stringsAsFactors=FALSE)
  }
  res <- if(kk>0) do.call(rbind, rows[seq_len(kk)]) else data.frame()
  cat(sprintf("[%s] clean-gap guides: %d of %d\n", NAME, kk, G))
  res
}

## ---- run all three datasets (paths from the shared registry; local display names kept) ----
.reg <- dataset_paths()
DATASETS <- list(
  c("replogle_rd7", .reg[["replogle-rd7"]]),
  c("endoc_t2d",    .reg[["endoc_t2d"]]),
  c("tcell_cd4",    .reg[["cd4_tcell"]])
)

all_res <- list()
for(ds in DATASETS){
  r <- analyze_dataset(ds[1], ds[2])
  if(!is.null(r) && nrow(r)>0) all_res[[ds[1]]] <- r
}
res <- do.call(rbind, all_res)

## ---- BH within each dataset over ALL clean-gap guides with a valid p ----
res$p_adj <- NA_real_
for(ds in unique(res$dataset)){
  idx <- which(res$dataset==ds & is.finite(res$p))
  res$p_adj[idx] <- p.adjust(res$p[idx], method="BH")
}
res$suffers <- is.finite(res$p_adj) & res$p_adj < 0.05

## ============================================================================
## Cross-check: reproduce pipeline obs2 / exp2_den for the SAME guides.
## ============================================================================
cat("\n===== CROSS-CHECK vs pipeline perguide CSVs (obs2, exp2_den) =====\n")
xcheck_files <- c(replogle_rd7="results/global_ambient_poisson/perguide_replogle_rd7.csv",
                  endoc_t2d   ="results/global_ambient_poisson/perguide_endoc_t2d.csv",
                  tcell_cd4   ="results/global_ambient_poisson/perguide_tcell_cd4.csv")
for(ds in names(xcheck_files)){
  f <- xcheck_files[[ds]]
  if(!file.exists(f)){ cat(sprintf("[%s] cross-check CSV missing: %s\n", ds, f)); next }
  ref <- read.csv(f, stringsAsFactors=FALSE)
  mine <- res[res$dataset==ds, ]
  m <- merge(mine[,c("guide","obs2","exp2_den")], ref[,c("guide","obs2","exp2_den")],
             by="guide", suffixes=c("_mine","_ref"))
  if(nrow(m)==0){ cat(sprintf("[%s] no overlapping guides for cross-check\n", ds)); next }
  obs2_match <- all(m$obs2_mine == m$obs2_ref, na.rm=TRUE)
  ## exp2_den rounded to 1 dp in ref; compare within tolerance
  exp2_err <- max(abs(m$exp2_den_mine - m$exp2_den_ref), na.rm=TRUE)
  cat(sprintf("[%s] %d guides matched | obs2 identical: %s | max |exp2_den diff|: %.4f\n",
      ds, nrow(m), obs2_match, exp2_err))
}

## ============================================================================
## Summary: at primary threshold gap_lo>=5 & n_amb>=50, + sensitivity ladder.
## ============================================================================
summarize_at <- function(res, gaplo_thresh, n_amb_thresh=50){
  out <- list(); i <- 0L
  for(ds in unique(res$dataset)){
    sub <- res[res$dataset==ds & res$gap_lo >= gaplo_thresh & res$n_amb >= n_amb_thresh & is.finite(res$p), ]
    n_wide <- nrow(sub)
    n_suf  <- sum(sub$suffers)
    pct    <- if(n_wide>0) 100*n_suf/n_wide else NA
    med_oe <- if(n_suf>0) median(sub$oe_region[sub$suffers], na.rm=TRUE) else NA
    i <- i+1L
    out[[i]] <- data.frame(dataset=ds, gap_lo_thresh=gaplo_thresh,
      n_wide_gap=n_wide, n_suffers=n_suf, pct_suffers=round(pct,1),
      median_oe_among_suffers=round(med_oe,2), stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}

ladder <- c(4,5,7,10)
summ <- do.call(rbind, lapply(ladder, function(t) summarize_at(res, t)))
summ <- summ[order(match(summ$dataset, c("replogle_rd7","endoc_t2d","tcell_cd4")), summ$gap_lo_thresh), ]

cat("\n===== SUMMARY (n_amb>=50; sensitivity ladder over gap_lo) =====\n")
print(summ, row.names=FALSE)

cat("\n===== PRIMARY (gap_lo>=5) =====\n")
print(summ[summ$gap_lo_thresh==5, ], row.names=FALSE)

## ---- write outputs ----
per_guide_out <- res[order(res$dataset, -res$suffers, res$p),
  c("dataset","guide","gap_lo","gap_hi","n_amb","obs_region","exp_region_den",
    "oe_region","p","p_adj","suffers")]
write.csv(per_guide_out, file.path(OUT, "weak_integration_prevalence.csv"), row.names=FALSE)
write.csv(summ, file.path(OUT, "weak_integration_prevalence_summary.csv"), row.names=FALSE)
cat(sprintf("\nwrote %s\n", file.path(OUT, "weak_integration_prevalence.csv")))
cat(sprintf("wrote %s\n", file.path(OUT, "weak_integration_prevalence_summary.csv")))

## ---- summary bar figure: percent suffering per dataset at primary threshold ----
prim <- summ[summ$gap_lo_thresh==5, ]
prim <- prim[order(match(prim$dataset, c("endoc_t2d","replogle_rd7","tcell_cd4"))), ]
png(file.path(OUT, "weak_integration_prevalence.png"), width=1600, height=1100, res=200)
op <- par(mar=c(5,5,4,2))
bp <- barplot(prim$pct_suffers, names.arg=prim$dataset,
   ylim=c(0, max(100, max(prim$pct_suffers,na.rm=TRUE)*1.15)),
   col=c("#4C72B0","#DD8452","#55A868"), border=NA,
   ylab="% of wide-gap guides suffering weak integration",
   main="Weak-integration prevalence (gap_lo>=5, n_amb>=50)\nExpected: EndoC (dense soup) low -> CD4 (thin soup) high")
text(bp, prim$pct_suffers,
     labels=sprintf("%.0f%%\n(%d/%d)", prim$pct_suffers, prim$n_suffers, prim$n_wide_gap),
     pos=3, cex=0.9, xpd=NA)
par(op); dev.off()
cat(sprintf("wrote %s\n", file.path(OUT, "weak_integration_prevalence.png")))

cat("\nDONE.\n")
