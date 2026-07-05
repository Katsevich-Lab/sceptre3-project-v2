#!/usr/bin/env Rscript
# Verify the comprehensive sim spans the real separation range, and show the
# method difficulty gradient across its parameter grid.
suppressPackageStartupMessages({library(Matrix); library(ggplot2); library(dplyr); library(tidyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
GA   <- normalizePath(file.path(HERE, ".."))
DATA <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data")
source(file.path(GA, "guide-assignment-pipeline", "bin", "script", "lib", "threshold_methods.R"))

gm <- as(readRDS(file.path(DATA,"sims_comprehensive/sceptre/grna_matrix.rds")), "RsparseMatrix")
truth <- as(readRDS(file.path(DATA,"sims_comprehensive/true_pert_matrix.rds")), "RsparseMatrix")
gn <- rownames(gm)
moi <- vapply(strsplit(gn,"_"), `[`, character(1), 2)
mu  <- as.integer(sub("mu","",vapply(strsplit(gn,"_"), `[`, character(1), 3)))
th  <- vapply(strsplit(gn,"_"), `[`, character(1), 4)

# assignments for the 3 methods (whole-matrix for ambient)
Aamb <- as(ambient_test_assign(gm, q=0.05)$assignment_matrix, "RsparseMatrix")
jac1 <- function(pred, tr){u<-sum(pred|tr); if(u==0) NA else sum(pred&tr)/u}
rows <- lapply(seq_along(gn), function(g){
  counts <- as.numeric(gm[g,]); tr <- as.numeric(truth[g,])>0 & counts>0
  v <- smoothed_valley_threshold(counts)
  to <- otsu_threshold_log1p(counts)$t; tv <- v$t
  data.frame(moi=moi[g], mu=mu[g], th=th[g],
    sep_gap = if(isTRUE(v$ok)) log1p(v$mode2)-log1p(v$mode1) else NA_real_,
    otsu   = jac1(if(is.finite(to)) counts>=to else logical(length(counts)), tr),
    valley = jac1(if(is.finite(tv)) counts>=tv else logical(length(counts)), tr),
    ambient= jac1(as.numeric(Aamb[g,])>0, tr))
})
df <- bind_rows(rows)
agg <- df %>% group_by(moi,mu,th) %>% summarize(across(c(sep_gap,otsu,valley,ambient), ~round(mean(.,na.rm=TRUE),3)), .groups="drop")
cat("===== comprehensive sim: separation + method Jaccard by param group =====\n")
print(as.data.frame(agg), row.names=FALSE)
write.csv(agg, file.path(HERE,"results","comprehensive_sim_grid.csv"), row.names=FALSE)

cat(sprintf("\nsim sep_gap range: %.2f - %.2f (real data span ~0.7 - 5.5)\n",
            min(df$sep_gap,na.rm=TRUE), max(df$sep_gap,na.rm=TRUE)))

# plot: ambient/otsu/valley Jaccard vs mu_pert, faceted by MOI x theta
pl <- agg %>% pivot_longer(c(otsu,valley,ambient), names_to="method", values_to="jaccard")
p <- ggplot(pl, aes(mu, jaccard, color=method, shape=method)) +
  geom_line() + geom_point(size=2) + facet_grid(moi~paste0("theta=",th)) +
  scale_x_continuous(trans="log2", breaks=c(15,50,150,500)) +
  labs(title="Comprehensive sim: method Jaccard across the realistic gRNA range",
       subtitle="mu_pert spans real signal levels; rows=MOI regime (real lib sizes), cols=dispersion",
       x=expression(mu[pert]), y="mean per-guide Jaccard") +
  theme_bw(base_size=9) + theme(legend.position="top")
ggsave(file.path(HERE,"results","Comprehensive_Sim_Grid.png"), p, width=9, height=5, dpi=120)
cat("Wrote results/Comprehensive_Sim_Grid.png\n")
