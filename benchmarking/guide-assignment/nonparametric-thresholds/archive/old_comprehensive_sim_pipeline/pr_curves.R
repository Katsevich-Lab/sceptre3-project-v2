#!/usr/bin/env Rscript
# Precision/recall curves on the comprehensive-sim benchmark subset (full ground
# truth). Sweep the FDR knob for ours (ambient, BH q) and fishash (padj_cutoff);
# Otsu/valley are single operating points. Split by difficulty regime so we can
# read off whether, at matched precision, ours equals fishash. (geomux curve is
# added by combine step using its Python sweep.)
suppressPackageStartupMessages({library(Matrix); library(fishash); library(SummarizedExperiment)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
EXP  <- file.path(HERE, "results", "comprehensive_bench")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

gm <- readRDS(file.path(EXP,"gm.rds")); truth <- readRDS(file.path(EXP,"truth.rds"))
guides <- read.csv(file.path(EXP,"guides.csv"))
Ncell <- ncol(gm); Ng <- nrow(gm)
mu <- as.integer(sub("mu","", vapply(strsplit(guides$group,"_"), `[`, character(1), 2)))
regime <- ifelse(mu <= 15, "hard (mu<=15)", ifelse(mu >= 150, "easy (mu>=150)", "moderate (mu=50)"))

# truth pairs (perturbed & observed) as keys, and per-guide regime
tt <- as(truth, "TsparseMatrix"); gg <- as(gm, "TsparseMatrix")
obs_key <- (gg@i)*Ncell + gg@j           # observed nonzero pairs
tru_full <- (tt@i)*Ncell + tt@j
truth_key <- intersect(tru_full, obs_key)  # perturbed AND observed
guide_of <- function(key) key %/% Ncell + 1L
truth_g  <- guide_of(truth_key)

pr_one <- function(assigned_key, assigned_g) {   # micro-averaged P/R per regime
  do.call(rbind, lapply(unique(regime), function(rg) {
    gid <- which(regime == rg)
    A <- assigned_key[assigned_g %in% gid]; T <- truth_key[truth_g %in% gid]
    tp <- length(intersect(A, T))
    data.frame(regime = rg, precision = if(length(A)) tp/length(A) else NA_real_,
               recall = if(length(T)) tp/length(T) else NA_real_, n_assigned = length(A))
  }))
}

rows <- list()
# ours: sweep BH q
amb <- ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=1)
ai <- amb$i; aj <- amb$j; ap <- amb$pval; akey <- (ai-1L)*Ncell + (aj-1L); ag <- ai  # 0-based key
padj <- p.adjust(ap, "BH")
for (q in c(.001,.005,.01,.02,.05,.1,.2,.35,.5,.7)) {
  sel <- padj < q
  pr <- pr_one(akey[sel], ag[sel]); pr$method <- "ours (ambient)"; pr$cutoff <- q; rows[[length(rows)+1]] <- pr
}
# fishash: sweep padj_cutoff (default method GS = block-dependent correction)
for (c0 in c(.01,.05,.1,.2,.35,.5,.7)) {
  fr <- fishash(gm, padj_cutoff=c0); cd <- as.data.frame(colData(fr)); cd$bc <- rownames(cd)
  asn <- strsplit(as.character(cd$assignment), ","); ci <- match(rep(cd$bc, lengths(asn)), colnames(gm))
  gi <- match(unlist(asn), rownames(gm)); ok <- !is.na(gi)&!is.na(ci)
  key <- (gi[ok]-1L)*Ncell + (ci[ok]-1L)   # 0-based key
  pr <- pr_one(key, gi[ok]); pr$method <- "fishash"; pr$cutoff <- c0; rows[[length(rows)+1]] <- pr
  cat("fishash padj", c0, "done\n")
}
# otsu / valley single points
for (mn in c("otsu","valley")) {
  A <- as(readRDS(file.path(EXP, paste0("assign_", mn, ".rds"))), "TsparseMatrix")
  key <- (A@i)*Ncell + A@j; g <- A@i + 1L
  pr <- pr_one(key, g); pr$method <- mn; pr$cutoff <- NA; rows[[length(rows)+1]] <- pr
}
df <- do.call(rbind, rows)
write.csv(df, file.path(HERE,"results","pr_curves.csv"), row.names=FALSE)
cat("\nwrote results/pr_curves.csv\n"); print(head(df[df$method=="ours (ambient)",],6), row.names=FALSE)
