#!/usr/bin/env Rscript
# ============================================================================
# Sensitivity of fishash+ to the INITIALIZATION of the assigned set A.
#   default : A_0 = {}          -> background_1 = raw counts (initial depth = raw library size)
#   clns    : A_0 = {N_gc > 2}  -> background_1 masks >2 (initial depth = CLEANSER-style <=2 depth)
# Both then run the identical fishash+ iteration to a fixed point. We compare the FINAL assigned
# sets: if they coincide (Jaccard ~ 1) the fixed point is initialization-insensitive.
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(fishash); library(extraDistr)
  library(sparseMatrixStats); library(ggplot2) })
source("scripts/contingency_method.R")
OUT <- "results/ambient_ceiling"

source("scripts/datasets.R")
paths <- dataset_paths()

ceil_per_guide <- function(counts, A) {
  ct <- as(counts,"TsparseMatrix"); i<-ct@i+1L; j<-ct@j+1L; y<-as.numeric(ct@x)
  a <- as.logical(A[cbind(i,j)]); a[is.na(a)]<-FALSE; amb <- !a
  as.numeric(tapply(y[amb], factor(i[amb], levels=seq_len(nrow(counts))), max))
}

summ <- list(); ceilrows <- list()
for (nm in names(paths)) {
  t0 <- proc.time()["elapsed"]
  counts <- load_grna_matrix(paths[[nm]])
  cat(sprintf("[%s] %s  %d x %d ...\n", format(Sys.time(),"%H:%M:%S"), nm, nrow(counts), ncol(counts))); flush.console()
  Adef <- fishashplus_assign(counts)$assigned                          # A0 = {} (default)
  Acln <- fishashplus_assign(counts, init_assigned = counts > 2)$assigned   # A0 = {N>2} (CLEANSER-depth seed)
  Adef <- as(Adef,"lgCMatrix"); Acln <- as(Acln,"lgCMatrix")
  inter <- sum(Adef & Acln); uni <- sum(Adef | Acln)
  jacc  <- if (uni>0) inter/uni else 1
  cd <- ceil_per_guide(counts, Adef); cc <- ceil_per_guide(counts, Acln)
  ok <- is.finite(cd) & is.finite(cc)
  summ[[nm]] <- data.frame(dataset=nm, n_default=sum(Adef), n_clns_init=sum(Acln),
    intersection=inter, union=uni, jaccard=jacc, n_flipped=uni-inter,
    ceil_max_absdiff = if(any(ok)) max(abs(cd[ok]-cc[ok])) else NA,
    ceil_mean_absdiff = if(any(ok)) mean(abs(cd[ok]-cc[ok])) else NA)
  ceilrows[[nm]] <- data.frame(dataset=nm, guide=seq_along(cd), ceil_default=cd, ceil_clns=cc)
  cat(sprintf("   default=%d  clns-init=%d  jaccard=%.4f  flipped=%d  (%.1fs)\n",
      sum(Adef), sum(Acln), jacc, uni-inter, proc.time()["elapsed"]-t0)); flush.console()
  saveRDS(list(summ=summ, ceilrows=ceilrows), file.path(OUT,"_init_progress.rds"))
  rm(counts, Adef, Acln); gc()
}
S <- do.call(rbind, summ); write.csv(S, file.path(OUT,"init_sensitivity_summary.csv"), row.names=FALSE)
CE <- do.call(rbind, ceilrows)

cat("\n==== INIT SENSITIVITY (default A={} vs CLEANSER-depth A={N>2}) ====\n"); print(S, row.names=FALSE)

# figure: Jaccard of final assigned sets per dataset
S$dataset <- factor(S$dataset, levels=S$dataset[order(S$jaccard)])
p1 <- ggplot(S, aes(jaccard, dataset)) +
  geom_col(fill="#3b6ea5") + geom_vline(xintercept=1, linetype="dashed", colour="grey50") +
  geom_text(aes(label=sprintf("%.3f", jaccard)), hjust=-0.1, size=3) +
  coord_cartesian(xlim=c(min(0.9,min(S$jaccard)*0.99),1.02)) +
  labs(title="fishash+ initialization sensitivity: final-assignment agreement",
       subtitle="Jaccard between fixed points from A0={} (default) vs A0={N>2} (CLEANSER-depth seed). 1 = identical.",
       x="Jaccard( default , CLEANSER-init )  [note x-axis starts near 0.9]", y=NULL) +
  theme_bw(base_size=12)
ggsave(file.path(OUT,"init_sensitivity_jaccard.png"), p1, width=9, height=6, dpi=150)

# per-guide ceiling agreement (do the seeds change the plausibly-ambient ceilings?)
CEp <- CE[is.finite(CE$ceil_default) & is.finite(CE$ceil_clns), ]
p2 <- ggplot(CEp, aes(ceil_default, ceil_clns)) +
  geom_abline(slope=1, intercept=0, colour="grey55", linetype="dashed") +
  geom_jitter(width=0.15, height=0.15, alpha=0.25, size=0.5, colour="#3b6ea5") +
  facet_wrap(~dataset, scales="free") +
  labs(title="Per-guide ambient ceiling: default init vs CLEANSER-depth init",
       x="ceiling (default A0={})", y="ceiling (CLEANSER-depth A0={N>2})") +
  theme_bw(base_size=10)
ggsave(file.path(OUT,"init_sensitivity_ceilings.png"), p2, width=14, height=10, dpi=150)
cat("\nwrote init_sensitivity_summary.csv + init_sensitivity_jaccard.png + init_sensitivity_ceilings.png\n")
