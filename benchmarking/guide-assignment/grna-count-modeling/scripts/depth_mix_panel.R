#!/usr/bin/env Rscript
# Cheap per-guide MARGINAL DEPTH-MIXED POISSON panel from the cached fishash+ fit.
# No refit: reads results/ambient_ceiling/fit_cache/<ds>_fit.rds (Rg, Cc, Tn, guide_mean) and
# one row of the original count matrix. lambda_{gc} = (Rg[g]/Tn) * Cc[c]; the fitted marginal
# predicted #cells with count k = sum_c dpois(k, lambda_{gc}).
#
# Usage:  source("scripts/depth_mix_panel.R")
#         depth_mix_panel("replogle-rd7", "8832_TFAM_P1P2_ENSG00000108064")   # guide id or index
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(tidyr) })

.DS_PATHS <- local({
  IN  <- "/Users/ekatsevi/data/projects/sceptre3/benchmarking/guide_assignment/input_data"
  SUR <- "/Users/ekatsevi/data/external/perturbseq-survey"
  list(barnyard_lrb100_0hr=file.path(IN,"barnyard_lrb100_0hr/sceptre/grna_matrix.rds"),
       barnyard_lrb100_72hr=file.path(IN,"barnyard_lrb100_72hr/sceptre/grna_matrix.rds"),
       barnyard_mch2_0hr=file.path(IN,"barnyard_mch2_0hr/sceptre/grna_matrix.rds"),
       barnyard_mch2_72hr=file.path(IN,"barnyard_mch2_72hr/sceptre/grna_matrix.rds"),
       gasperini=file.path(IN,"gasperini/sceptre/grna_matrix.rds"),
       `replogle-rd7`=file.path(IN,"replogle-rd7/sceptre/grna_matrix.rds"),
       mccutcheon="/Users/ekatsevi/data/external/mccutcheon-2023-GSE218988/grna_matrix.rds",
       a549=file.path(SUR,"a549_crispri_perturbseq_GSE236304/grna_matrix.rds"),
       cd8_tcell=file.path(SUR,"cd8_tcell_perturbcite_GSE279498/grna_matrix.rds"),
       dctap_k562_highmoi=file.path(SUR,"dctap_k562_highmoi_GSE303901/grna_matrix.rds"),
       dctap_k562_lowmoi=file.path(SUR,"dctap_k562_lowmoi_GSE303901/grna_matrix.rds"),
       endoc_t2d=file.path(SUR,"endoc_t2d_perturbseq_GSE273677/grna_matrix.rds"),
       gastric_organoid=file.path(SUR,"gastric_organoid_cropseq_GSE280506/grna_matrix.rds"),
       invivo_cortex=file.path(SUR,"invivo_cortex_perturbseq_GSE249416/grna_matrix.rds"),
       ipsc=file.path(SUR,"ipsc_crispri_hipsci_figshare27989294/grna_matrix.rds"),
       erythroid_multiome=file.path(SUR,"perturb_multiome_erythroid_GSE274113/grna_matrix.rds"),
       cd4_tcell=file.path(SUR,"tcell_cd4_perturbseq_GSE314342/grna_matrix.rds"))
})
.CACHE <- "results/ambient_ceiling/fit_cache"

depth_mix_panel <- function(dataset, guide, KMAX = 25L, save = TRUE, outdir = "results/ambient_ceiling") {
  fit <- readRDS(file.path(.CACHE, paste0(dataset, "_fit.rds")))
  g   <- if (is.numeric(guide)) as.integer(guide) else match(guide, fit$guides)
  if (is.na(g)) stop("guide not found in ", dataset)
  gid <- fit$guides[g]
  # observed counts: read just this guide's row from the original matrix
  counts <- readRDS(.DS_PATHS[[dataset]])
  y   <- as.numeric(counts[g, ]); rm(counts)
  lam <- (fit$Rg[g] / fit$Tn) * fit$Cc                     # lambda_{gc} = a_g d_c, all cells
  mu  <- fit$guide_mean[g]; n <- length(y)

  ks  <- 0:KMAX
  obs <- vapply(ks, function(k) if (k<KMAX) sum(y==k) else sum(y>=k), numeric(1))
  dm  <- vapply(ks, function(k) if (k<KMAX) sum(dpois(k,lam)) else sum(1-ppois(k-1,lam)), numeric(1))
  sp  <- vapply(ks, function(k) if (k<KMAX) n*dpois(k,mu)    else n*(1-ppois(k-1,mu)),   numeric(1))
  klab<- ifelse(ks<KMAX, as.character(ks), paste0(KMAX,"+"))
  # ambient ceiling straight from the cached per-guide table for consistency
  pgt <- read.csv(file.path(outdir, "per_guide_ambient_ceiling.csv"))
  ceil<- pgt$ambient_ceiling[pgt$dataset==dataset & pgt$guide==gid]
  ceil<- if (length(ceil)) ceil[1] else max(which(obs>0))-1

  D  <- data.frame(k=ks, klab=factor(klab, levels=klab), observed=obs, depthmix=dm, single=sp)
  Dl <- pivot_longer(D, c(depthmix, single), names_to="model", values_to="pred")
  Dl$model <- factor(Dl$model, c("depthmix","single"),
    labels=c("fitted marginal depth-mixed Poisson (a_g d_c)","single Poisson (guide mean)"))
  p <- ggplot(D, aes(klab, observed)) +
    geom_col(fill="grey82", width=0.72) +
    geom_vline(xintercept=ceil+1, linetype="dashed", colour="grey40") +
    annotate("text", x=ceil+1, y=max(obs), hjust=-0.08, vjust=1,
             label=sprintf("ambient ceiling = %d", ceil), size=4, colour="grey30") +
    geom_line(data=Dl, aes(klab, pred, colour=model, group=model), linewidth=0.8) +
    geom_point(data=Dl, aes(klab, pred, colour=model), size=2.2) +
    scale_y_continuous(trans=scales::trans_new("log1p",log1p,expm1),
                       breaks=c(0,1,10,100,1000,10000,1e5,1e6), labels=scales::label_comma(),
                       expand=expansion(mult=c(0,0.05))) +
    scale_colour_manual(values=c("fitted marginal depth-mixed Poisson (a_g d_c)"="#D55E00",
                                 "single Poisson (guide mean)"="#0072B2")) +
    labs(title=sprintf("%s  %s", dataset, gid),
         subtitle=sprintf("var/mean=%.1f  ceiling=%d  (depth-mixed Poisson = ambient model)",
                          var(y)/mean(y), ceil),
         x="gRNA UMI count", y="number of cells (log1p; floor 0)", colour=NULL) +
    theme_bw(base_size=13) + theme(legend.position="top", axis.text.x=element_text(size=9))
  if (save) {
    f <- file.path(outdir, sprintf("panel_%s_%s.png", dataset, gsub("[^A-Za-z0-9]+","_", gid)))
    ggsave(f, p, width=11, height=6, dpi=140); cat("wrote", f, "\n")
  }
  invisible(list(plot=p, table=D, ceiling=ceil, vmr=var(y)/mean(y)))
}

# CLI: Rscript scripts/depth_mix_panel.R <dataset> <guide>
if (sys.nframe() == 0) {
  a <- commandArgs(trailingOnly=TRUE)
  if (length(a) >= 2) depth_mix_panel(a[1], a[2])
}
