#!/usr/bin/env Rscript
# Headline benchmark: full method panel across the comprehensive sim grid.
# Methods: ours (ambient), Otsu, valley, real geomux (default & nolor), real
# fishash. Per-guide Jaccard/precision/recall vs ground truth, aggregated by the
# (MOI, mu_pert, theta_pert) parameter group.
suppressPackageStartupMessages({library(Matrix); library(fishash); library(SummarizedExperiment)
  library(ggplot2); library(dplyr); library(tidyr)})
HERE <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
HERE <- dirname(HERE)  # scripts/ -> folder root (added by reorg)
EXP  <- file.path(HERE, "results", "comprehensive_bench")

gm <- readRDS(file.path(EXP,"gm.rds")); truth <- readRDS(file.path(EXP,"truth.rds"))
guides <- read.csv(file.path(EXP,"guides.csv")); cells <- read.csv(file.path(EXP,"cells.csv"))
gnames <- rownames(gm); gmr <- as(gm,"RsparseMatrix"); tr_r <- as(truth,"RsparseMatrix")
truth_sets <- lapply(seq_along(gnames), function(g) which(as.numeric(tr_r[g,])>0 & as.numeric(gmr[g,])>0))

# per-guide predicted sets, from a guides x cells matrix
sets_from_matrix <- function(A){A<-as(A,"RsparseMatrix"); lapply(seq_along(gnames), function(g) which(as.numeric(A[g,])>0))}
# per-guide predicted sets, from a (cell_barcode, guide) call table
sets_from_calls <- function(fp){
  d <- read.csv(fp); ci <- match(as.character(d$cell_barcode), cells$barcode); gi <- match(as.character(d$guide), gnames)
  ok <- !is.na(ci)&!is.na(gi); s <- split(ci[ok], factor(gi[ok], levels=seq_along(gnames)))
  lapply(s, function(x) if(is.null(x)) integer(0) else x)
}
metrics <- function(psets){
  do.call(rbind, lapply(seq_along(gnames), function(g){
    p<-psets[[g]]; t<-truth_sets[[g]]; i<-length(intersect(p,t)); u<-length(union(p,t))
    data.frame(jac=if(u==0)NA else i/u, prec=if(length(p))i/length(p) else NA, rec=if(length(t))i/length(t) else NA)
  }))
}

# fishash (run here)
cat("running fishash...\n")
fr <- fishash(gm); cd <- as.data.frame(colData(fr)); cd$bc <- rownames(cd)
asn <- strsplit(as.character(cd$assignment), ","); fc <- data.frame(cell_barcode=rep(cd$bc, lengths(asn)), guide=unlist(asn))
fc <- fc[!is.na(fc$guide) & fc$guide!="",]; write.csv(fc, file.path(EXP,"fishash_calls.csv"), row.names=FALSE)

panel <- list(
  ours    = sets_from_matrix(readRDS(file.path(EXP,"assign_ambient.rds"))),
  otsu    = sets_from_matrix(readRDS(file.path(EXP,"assign_otsu.rds"))),
  valley  = sets_from_matrix(readRDS(file.path(EXP,"assign_valley.rds"))),
  geomux  = sets_from_calls(file.path(EXP,"geomux_default_calls.csv")),
  `geomux_nolor` = sets_from_calls(file.path(EXP,"geomux_nolor_calls.csv")),
  fishash = sets_from_calls(file.path(EXP,"fishash_calls.csv")))

res <- bind_rows(lapply(names(panel), function(mn){ m<-metrics(panel[[mn]]); m$method<-mn; m$guide<-gnames; m }))
res <- res %>% left_join(guides, by=c("guide")) %>%
  separate(group, into=c("moi","mu","th"), sep="_") %>% mutate(mu=as.integer(sub("mu","",mu)))
agg <- res %>% group_by(method,moi,mu,th) %>%
  summarize(jaccard=mean(jac,na.rm=TRUE), precision=mean(prec,na.rm=TRUE), recall=mean(rec,na.rm=TRUE), .groups="drop")
overall <- res %>% group_by(method) %>% summarize(jaccard=round(mean(jac,na.rm=TRUE),3),
  precision=round(mean(prec,na.rm=TRUE),3), recall=round(mean(rec,na.rm=TRUE),3), .groups="drop") %>% arrange(-jaccard)
cat("\n===== overall mean across the comprehensive grid =====\n"); print(as.data.frame(overall), row.names=FALSE)
write.csv(agg, file.path(HERE,"results","comprehensive_panel_grid.csv"), row.names=FALSE)
write.csv(overall, file.path(HERE,"results","comprehensive_panel_overall.csv"), row.names=FALSE)

p <- ggplot(agg, aes(mu, jaccard, color=method, shape=method)) +
  geom_line() + geom_point(size=1.7) + facet_grid(moi~paste0("theta=",th)) +
  scale_x_continuous(trans="log2", breaks=c(5,15,50,150,500)) +
  labs(title="Full method panel across the comprehensive sim (real gRNA parameter range)",
       subtitle="per-guide Jaccard vs ground truth; rows=MOI regime, cols=dispersion; x=signal level mu_pert",
       x=expression(mu[pert]), y="mean Jaccard") + theme_bw(base_size=9) + theme(legend.position="top")
ggsave(file.path(HERE,"results","Comprehensive_Panel.png"), p, width=10, height=5.5, dpi=120)
cat("\nWrote results/Comprehensive_Panel.png + grid/overall csv\n")
