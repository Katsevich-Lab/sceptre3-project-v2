#!/usr/bin/env Rscript
# ============================================================================
# DC-TAP K562 high-MOI positive-control comparison: fishash vs fishash+ vs CLEANSER.
# 12 positive-control guides target GATA1 / HDAC6 / MYC (TSS + enhancer, 4 each), all in the
# 93-gene DC-TAP GEX panel. For each method we assign guides to cells, then ask: do the assigned
# cells actually show knockdown of the target gene? A better assignment recovers MORE truly-
# perturbed cells (recall) WITHOUT diluting the knockdown (precision).
#
# fishash  = contingency_assign(cell_margin="observed")   (reproduces fishash)
# fishash+ = contingency_assign(cell_margin="ambient")    (depth_fix)
# CLEANSER = run_external_assign.py cleanser --dc  -> calls.csv  (guide,cell pairs)
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2) })
source("scripts/contingency_method.R"); source("scripts/datasets.R")
DS  <- dctap_sceptre_dir("highmoi")
OUT <- "results/ambient_ceiling"

gm  <- load_grna_matrix(file.path(DS,"grna_matrix_aligned.rds"))
rm_ <- as(readRDS(file.path(DS,"response_matrix.rds")),"CsparseMatrix")
tdf <- read.csv(file.path(DS,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
stopifnot(identical(colnames(gm), colnames(rm_)))
gene_of <- setNames(tdf$grna_target, tdf$grna_id)                 # ENSG target per guide
sym_of  <- setNames(tdf$grna_target_symbol, tdf$grna_id)
pos_guides <- tdf$grna_id[tdf$grna_target_symbol != "non-targeting"]
nt_guides  <- tdf$grna_id[tdf$grna_target_symbol == "non-targeting"]

## ---- assignments -----------------------------------------------------------
cat("running fishash (observed) + fishash+ (ambient)...\n")
A_fish  <- as(fishash_assign(gm)$assigned, "lgCMatrix")
A_fishp <- as(fishashplus_assign(gm)$assigned, "lgCMatrix")
A_clns  <- as(load_cleanser_assignment(gm, "dctap_k562_highmoi"), "lgCMatrix")
methods <- list(fishash=A_fish, `fishash+`=A_fishp, CLEANSER=A_clns)
cat(sprintf("total assigned: fishash=%d  fishash+=%d  CLEANSER=%d\n",
            sum(A_fish), sum(A_fishp), sum(A_clns)))

## ---- normalized target expression -----------------------------------------
lib <- Matrix::colSums(rm_); sf <- lib / median(lib)             # per-cell size factor
norm_expr <- function(gene) { x <- as.numeric(rm_[gene, ]); log1p(x / sf * 1) }  # log1p CP-median

## ---- per-guide knockdown per method ---------------------------------------
# METHOD-INDEPENDENT per-gene negative control: cells with NO real integration of the target gene's
# guides (max raw count over that gene's 4 guides <= 2). Identical baseline for all three methods, so
# only the treated (assigned) set varies -> a fair cross-method comparison. NOTE: at MOI ~9 a strict
# "zero targeting guides at all" control is only ~540 cells and, if defined per-method, shrinks for
# high-recall methods (fishash+ 1472 vs fishash 2905) and is biased by their own assignment; the
# per-gene count<=2 control keeps ~1.6-2k cells and is method-independent.
genes <- unique(unname(gene_of[pos_guides]))
ctrl_fixed <- lapply(setNames(genes, genes), function(G){
  gg <- names(gene_of)[gene_of == G]
  as.numeric(sparseMatrixStats::colMaxs(gm[gg, , drop=FALSE])) <= 2
})
res <- list()
for (mn in names(methods)) {
  A <- methods[[mn]]
  for (g in pos_guides) {
    G <- gene_of[g]; e <- norm_expr(G)
    treated <- as.logical(A[g, ])
    control <- ctrl_fixed[[G]] & !treated                        # no-target-integration baseline
    nt <- sum(treated)
    if (nt >= 3 && sum(control) >= 20) {
      mt <- mean(expm1(e[treated])); mc <- mean(expm1(e[control]))
      lfc <- log2((mt + 1e-6) / (mc + 1e-6))
      tt  <- tryCatch(t.test(e[treated], e[control]), error=function(x) NULL)
      tstat <- if(!is.null(tt)) unname(tt$statistic) else NA
      p     <- if(!is.null(tt)) tt$p.value else NA
    } else { lfc<-NA; tstat<-NA; p<-NA }
    res[[length(res)+1]] <- data.frame(method=mn, guide=g, gene=sym_of[g],
      type=ifelse(grepl("TSS",g),"TSS","enhancer"),
      n_assigned=nt, log2fc=lfc, tstat=tstat, neglog10p=-log10(p))
  }
}
R <- do.call(rbind, res)
write.csv(R, file.path(OUT,"dctap_poscontrol_compare.csv"), row.names=FALSE)
cat("\n=== per-guide (target-gene knockdown; negative log2fc = knockdown) ===\n")
print(R[order(R$gene,R$guide,R$method), c("gene","guide","type","method","n_assigned","log2fc","neglog10p")],
      row.names=FALSE, digits=3)

R$method <- factor(R$method, levels=c("fishash","fishash+","CLEANSER"))
cols <- c(fishash="#7f8c8d", `fishash+`="#c0392b", CLEANSER="#2980b9")  # grey / red / blue (distinct)

# cap -log10p: p can underflow to 0 (Inf) for the strongest knockdowns -> keep those bars, capped+starred
fin <- R$neglog10p[is.finite(R$neglog10p)]
CAP <- ceiling(max(fin, na.rm=TRUE)/50)*50            # e.g. 350
R$power <- pmin(R$neglog10p, CAP); R$power[is.infinite(R$neglog10p)] <- CAP
R$underflow <- is.infinite(R$neglog10p)               # p below double precision (~1e-308)

## (A) recall vs precision: n cells assigned (x) vs knockdown log2FC (y), per guide
pA <- ggplot(R[is.finite(R$log2fc),], aes(n_assigned, log2fc, colour=method, shape=type)) +
  geom_hline(yintercept=0, colour="grey60") +
  geom_line(aes(group=guide), colour="grey70", linewidth=0.3) +
  geom_point(size=2.8) +
  facet_wrap(~gene, scales="free_x") +
  scale_colour_manual(values=cols) +
  labs(title="DC-TAP high-MOI positive controls: recall vs knockdown, by method",
       subtitle="Each point = one guide under one method. Lower y = stronger target knockdown; higher x = more cells assigned. Lines join the same guide across methods.",
       x="cells assigned to the guide", y="target-gene log2FC (assigned vs no-target-integration control)") +
  theme_bw(base_size=12)
ggsave(file.path(OUT,"dctap_poscontrol_recall_vs_knockdown.png"), pA, width=12, height=5, dpi=140)

## (B) detection power: -log10 p of target knockdown per guide (Inf capped at CAP, starred)
pB <- ggplot(R[!is.na(R$power),], aes(guide, power, fill=method)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  geom_text(aes(label=ifelse(underflow,"*","")), position=position_dodge(0.8),
            vjust=-0.1, size=5, show.legend=FALSE) +
  facet_wrap(~gene, scales="free_x") +
  scale_fill_manual(values=cols) +
  labs(title="DC-TAP high-MOI positive controls: knockdown detection power",
       subtitle=sprintf("-log10 p (Welch t), assigned vs no-target-integration control. Higher = stronger detection. * = p underflowed to 0 (capped at %g).", CAP),
       x=NULL, y="-log10 p (target knockdown)") +
  theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(OUT,"dctap_poscontrol_power.png"), pB, width=12, height=5.5, dpi=140)

## (C) estimated knockdown per guide: PERCENTAGE change in target expression, vertical bars
##     (parallels the power plot: guide on x, facet by gene, dodge by method)
R$pct <- (2^R$log2fc - 1) * 100                        # % change in target expr (negative = knockdown)
pC <- ggplot(R[is.finite(R$pct),], aes(guide, pct, fill=method)) +
  geom_col(position=position_dodge(0.8), width=0.75) +
  geom_hline(yintercept=0, colour="grey40") +
  geom_text(aes(label=sprintf("%.0f%%", pct)), position=position_dodge(0.8),
            vjust=1.25, size=2.7) +
  facet_wrap(~gene, scales="free_x") +
  scale_fill_manual(values=cols) +
  labs(title="DC-TAP high-MOI positive controls: estimated target-gene knockdown",
       subtitle="Percentage change in target-gene expression, assigned vs no-target-integration control cells (more negative = stronger knockdown).",
       x=NULL, y="target-gene expression change (%)") +
  theme_bw(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(OUT,"dctap_poscontrol_knockdown.png"), pC, width=12, height=5.5, dpi=140)

## aggregate summary per method
agg <- do.call(rbind, lapply(split(R, R$method), function(d){
  d<-d[is.finite(d$log2fc),]
  data.frame(method=d$method[1], guides_scored=nrow(d),
    total_assigned_cells=sum(d$n_assigned), median_log2fc=median(d$log2fc),
    median_neglog10p=median(d$neglog10p[is.finite(d$neglog10p)]),
    n_knockdown_detected=sum(d$neglog10p>2 & d$log2fc<0, na.rm=TRUE)) }))
write.csv(agg, file.path(OUT,"dctap_poscontrol_summary.csv"), row.names=FALSE)
cat("\n=== aggregate per method ===\n"); print(agg, row.names=FALSE, digits=3)
cat("\nwrote dctap_poscontrol_* to", OUT, "\n")
