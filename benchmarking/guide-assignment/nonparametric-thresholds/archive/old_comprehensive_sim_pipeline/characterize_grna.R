#!/usr/bin/env Rscript
# Characterize gRNA count distributions across diverse REAL datasets to derive
# sensible simulation parameter ranges for the sum-process simulator
#   y = Pois(lib*mu_exo) + [endo] NB(theta_endo, lib*mu_endo) + [pert] NB(theta_pert, lib*mu_pert)
# Per guide we split cells into perturbed/background with the ambient test, then
# estimate: prob_pert, library-normalized mu_pert, NB theta_pert, background rate,
# and mode separation. Aggregate per dataset -> realistic ranges.

suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
SURV <- path.expand("~/data/external/perturbseq-survey")
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||is.na(a)) b else a

# discover survey datasets + add real Gasperini/Replogle as references
ds <- list()
for (d in list.dirs(SURV, recursive=FALSE)) { f <- file.path(d,"grna_matrix.rds")
  if (file.exists(f)) ds[[basename(d)]] <- f }
ds[["gasperini (ref)"]]  <- file.path(DATA,"gasperini/sceptre/grna_matrix.rds")
ds[["replogle-med (ref)"]] <- file.path(DATA,"replogle-rd7_medium/sceptre/grna_matrix.rds")

characterize <- function(gm, max_cells=40000, max_guides=600, seed=1) {
  set.seed(seed)
  if (ncol(gm) > max_cells) gm <- gm[, sort(sample(ncol(gm), max_cells))]
  gm <- as(gm, "CsparseMatrix"); storage.mode(gm@x) <- "double"
  lib_nnz <- Matrix::colSums(gm > 0); lib_umi <- Matrix::colSums(gm)
  lib_scaled <- lib_nnz / mean(lib_nnz)
  A  <- ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=1)$assignment_matrix
  Ar <- as(A, "RsparseMatrix"); gmr <- as(gm, "RsparseMatrix")
  ng <- nrow(gm); gsel <- if (ng > max_guides) sort(sample(ng, max_guides)) else seq_len(ng)
  rows <- lapply(gsel, function(g) {
    counts <- as.numeric(gmr[g,]); asn <- as.numeric(Ar[g,]) > 0
    np <- sum(asn); if (np < 5) return(NULL)
    pc <- counts[asn]; bc <- counts[!asn]; lsp <- lib_scaled[asn]
    m <- mean(pc); v <- var(pc)
    data.frame(
      prob_pert  = np/length(counts),
      mu_pert    = mean(pc/lsp),                          # library-normalized signal level
      mu_pert_raw= m,
      theta_pert = if (v>m) min(m^2/(v-m), 1000) else 1000,  # NB dispersion (method of moments)
      bg_rate    = mean(bc/lib_scaled[!asn]),             # library-normalized background
      sep_gap    = log1p(median(pc)) - log1p(as.numeric(quantile(bc,.95))))  # signal above bg bulk
  })
  pg <- do.call(rbind, rows)
  list(per_guide=pg, n_guides=ng, n_cells=ncol(gm),
       moi_assigned=median(Matrix::colSums(A)), lib_nnz_med=median(lib_nnz), lib_umi_med=median(lib_umi))
}

summ <- list(); allpg <- list()
for (nm in names(ds)) {
  gm <- readRDS(ds[[nm]]); r <- characterize(gm)
  pg <- r$per_guide
  q2 <- function(x) round(median(x, na.rm=TRUE),2)
  summ[[nm]] <- data.frame(dataset=nm, n_guides=r$n_guides, n_cells=r$n_cells,
    moi=r$moi_assigned, lib_nnz=r$lib_nnz_med,
    prob_pert=signif(median(pg$prob_pert,na.rm=TRUE),2),
    mu_pert=q2(pg$mu_pert), theta_pert=q2(pg$theta_pert),
    bg_rate=signif(median(pg$bg_rate,na.rm=TRUE),2), sep_gap=q2(pg$sep_gap),
    n_guides_used=nrow(pg))
  allpg[[nm]] <- pg %>% mutate(dataset=nm)
  cat("done", nm, "\n")
}
S <- do.call(rbind, summ)
cat("\n===== per-dataset gRNA parameter estimates (medians) =====\n"); print(S, row.names=FALSE)
write.csv(S, file.path(HERE,"results","grna_param_summary.csv"), row.names=FALSE)

PG <- bind_rows(allpg)
saveRDS(PG, file.path(HERE,"results","grna_param_perguide.rds"))
cat("\n===== realistic ranges (10th-90th pct across all guides, all datasets) =====\n")
for (v in c("prob_pert","mu_pert","theta_pert","bg_rate","sep_gap")) {
  qs <- quantile(PG[[v]], c(.1,.5,.9), na.rm=TRUE)
  cat(sprintf("  %-11s 10%%=%-8.3g median=%-8.3g 90%%=%-8.3g\n", v, qs[1], qs[2], qs[3]))
}

# plots: distribution of each parameter by dataset
pl <- PG %>% tidyr::pivot_longer(c(mu_pert,theta_pert,prob_pert,sep_gap), names_to="param", values_to="val") %>%
  mutate(val = ifelse(param %in% c("mu_pert","theta_pert") & val>0, val, val))
p <- ggplot(pl, aes(dataset, val, fill=grepl("ref", dataset))) +
  geom_boxplot(outlier.size=0.2) + facet_wrap(~param, scales="free_y") +
  scale_fill_manual(values=c(`FALSE`="#2c7fb8",`TRUE`="grey70"), guide="none") +
  scale_y_continuous(trans=scales::pseudo_log_trans(base=10)) +
  labs(title="gRNA distribution parameters across diverse real datasets (per guide)",
       subtitle="blue = 2024-26 survey; grey = Gasperini/Replogle reference; y on pseudo-log",
       x=NULL, y="value (pseudo-log10)") +
  theme_bw(base_size=8) + theme(axis.text.x=element_text(angle=40,hjust=1))
ggsave(file.path(HERE,"results","grna_param_distributions.png"), p, width=11, height=7, dpi=120)
cat("\nWrote results/grna_param_summary.csv + grna_param_distributions.png\n")
