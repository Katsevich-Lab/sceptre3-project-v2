#!/usr/bin/env Rscript
# Test the recommended method (fishash + ambient-depth cell margin) on the realistic comprehensive
# sim grid (1200 guides x 12k cells; lowMOI/highMOI x mu_pert x theta), vs fishash and the bare
# ambient test. Per-guide Jaccard/precision/recall vs ground truth (perturbed AND observed),
# aggregated by parameter group. -> results/benchmark_update/comprehensive_depthfix.{png,csv}
suppressPackageStartupMessages({library(Matrix); library(fishash); library(sparseMatrixStats)
  library(ggplot2); library(dplyr); library(tidyr)})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))),error=function(e) ".")
HERE <- dirname(HERE); GA <- normalizePath(file.path(HERE,".."))
source(file.path(HERE,"scripts","contingency_method.R"))
source(file.path(GA,"guide-assignment-pipeline","bin","script","lib","threshold_methods.R"))
EXP <- file.path(HERE,"results","comprehensive_bench"); RES <- file.path(HERE,"results","benchmark_update")

gm <- readRDS(file.path(EXP,"gm.rds")); truth <- readRDS(file.path(EXP,"truth.rds"))
guides <- read.csv(file.path(EXP,"guides.csv")); gnames <- rownames(gm)
gmr <- as(gm,"RsparseMatrix"); tr_r <- as(truth,"RsparseMatrix")
truth_sets <- lapply(seq_along(gnames), function(g) which(as.numeric(tr_r[g,])>0 & as.numeric(gmr[g,])>0))
sets_from_matrix <- function(A){A<-as(A,"RsparseMatrix"); lapply(seq_along(gnames), function(g) which(as.numeric(A[g,])>0))}
metrics <- function(psets) do.call(rbind, lapply(seq_along(gnames), function(g){
  p<-psets[[g]]; t<-truth_sets[[g]]; i<-length(intersect(p,t)); u<-length(union(p,t))
  data.frame(jac=if(u==0)NA else i/u, prec=if(length(p))i/length(p) else NA, rec=if(length(t))i/length(t) else NA)}))

cat("running ambient test...\n");  amb <- ambient_test_assign(gm, q=0.05, model="hypergeometric", n_iter=1)$assignment_matrix
cat("running fishash (contingency, observed margin)...\n"); fsh <- contingency_assign(gm, cell_margin="observed", tail="hyper", fdr="GS")$assigned
cat("running proposed (contingency, ambient-depth margin)...\n"); prop<- contingency_assign(gm, cell_margin="ambient",  tail="hyper", fdr="GS")$assigned

panel <- list(`ambient test`=sets_from_matrix(amb), `fishash`=sets_from_matrix(fsh), `proposed (depth fix)`=sets_from_matrix(prop))
res <- bind_rows(lapply(names(panel), function(mn){m<-metrics(panel[[mn]]); m$method<-mn; m$guide<-gnames; m})) %>%
  left_join(guides, by="guide") %>% separate(group, into=c("moi","mu","th"), sep="_") %>% mutate(mu=as.integer(sub("mu","",mu)))
agg <- res %>% group_by(method,moi,mu,th) %>% summarize(jaccard=mean(jac,na.rm=TRUE),precision=mean(prec,na.rm=TRUE),recall=mean(rec,na.rm=TRUE),.groups="drop")
overall <- res %>% group_by(method) %>% summarize(jaccard=round(mean(jac,na.rm=TRUE),3),precision=round(mean(prec,na.rm=TRUE),3),recall=round(mean(rec,na.rm=TRUE),3),.groups="drop") %>% arrange(-jaccard)
byMOI <- res %>% group_by(method,moi) %>% summarize(jaccard=round(mean(jac,na.rm=TRUE),3),recall=round(mean(rec,na.rm=TRUE),3),.groups="drop")
write.csv(agg, file.path(RES,"comprehensive_depthfix.csv"), row.names=FALSE)
cat("\n=== overall (mean over the grid) ===\n"); print(as.data.frame(overall),row.names=FALSE)
cat("\n=== by MOI regime (where the depth fix should matter) ===\n"); print(as.data.frame(byMOI),row.names=FALSE)

agg$method<-factor(agg$method,levels=c("ambient test","fishash","proposed (depth fix)"))
cols<-c(`ambient test`="#888780", fishash="#1D9E75", `proposed (depth fix)`="#185FA5")
p<-ggplot(agg, aes(mu, jaccard, colour=method)) + geom_line(linewidth=0.7)+geom_point(size=1.8)+
  facet_grid(paste0("MOI: ",moi)~paste0("theta=",th)) + scale_x_continuous(trans="log2",breaks=c(5,15,50,150,500))+
  scale_colour_manual(values=cols)+
  labs(title="Recommended method vs fishash vs ambient test — comprehensive sim (realistic gRNA range)",
       subtitle="per-guide Jaccard vs ground truth; rows = MOI regime, cols = dispersion; x = signal level mu_pert",
       x=expression(mu[pert]), y="mean Jaccard")+theme_bw(base_size=10)+theme(legend.position="top",legend.title=element_blank())
ggsave(file.path(RES,"comprehensive_depthfix.png"), p, width=11, height=6, dpi=140)
cat("\nwrote comprehensive_depthfix.{png,csv}\n")
