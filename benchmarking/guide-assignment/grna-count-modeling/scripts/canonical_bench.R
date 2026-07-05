#!/usr/bin/env Rscript
# Ablate the two canonical fixes vs fishash across MOI (Opening 1: depth/cell-margin) and
# overdispersion (Opening 2: beta-binomial tail). Methods:
#   fishash        = contingency(observed margin, hyper)   [reproduces fishash]
#   + depth        = contingency(AMBIENT margin,  hyper)
#   + depth + NB    = contingency(AMBIENT margin,  betabinom)   [the canonical method]
suppressPackageStartupMessages({library(fishash); library(Matrix); library(SummarizedExperiment)
  library(sparseMatrixStats); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))),error=function(e) ".")
HERE <- dirname(HERE); source(file.path(HERE,"scripts","contingency_method.R"))
RES <- file.path(HERE,"results","benchmark_update")
score <- function(A,truth){A<-as(A,"lgCMatrix");truth<-as(truth,"lgCMatrix");na<-sum(A);nt<-sum(truth);tp<-sum(A&truth)
  prec<-if(na>0)tp/na else NA; rec<-if(nt>0)tp/nt else NA
  c(fdr=if(na>0)(na-tp)/na else 0, recall=rec, f1=if(!is.na(prec)&!is.na(rec)&(prec+rec)>0)2*prec*rec/(prec+rec) else NA)}
regimes <- list(
  `MOI 0.3`      =list(moi=0.3,phi=0), `MOI 1`=list(moi=1,phi=0),
  `MOI 3`        =list(moi=3,  phi=0), `MOI 5`=list(moi=5,phi=0),
  `MOI 0.3 +OD`  =list(moi=0.3,phi=1), `MOI 1 +OD`=list(moi=1,phi=1))
NREP<-4; Q<-0.05; rows<-list()
for (rg in names(regimes)){p<-regimes[[rg]]
  for (rep in seq_len(NREP)){ set.seed(200*which(names(regimes)==rg)+rep)
    sim<-simulate_guidebender2(n_guides=200,n_cells=4000,moi=p$moi,hurdle_prob=0.1,snr=4,
          count_per_cell=100,frac_noise_endo=0.75,endo_shape_flat=0.5,Phi_noise=p$phi,return_sparse_only=TRUE)
    N<-assay(sim,"counts"); tr<-assay(sim,"ground_truth")
    m<-list(
      `fishash`        =score(contingency_assign(N,cell_margin="observed",tail="hyper",fdr="GS")$assigned,tr),
      `+ depth`        =score(contingency_assign(N,cell_margin="ambient", tail="hyper",fdr="GS")$assigned,tr),
      `+ depth + NB`   =score(contingency_assign(N,cell_margin="ambient", tail="nb",fdr="GS")$assigned,tr))
    for (mn in names(m)) rows[[length(rows)+1]]<-data.frame(regime=rg,rep=rep,method=mn,
      fdr=m[[mn]][1],recall=m[[mn]][2],f1=m[[mn]][3])
  }
  s<-bind_rows(rows)%>%filter(regime==rg)%>%group_by(method)%>%summarise(fdr=mean(fdr),rec=mean(recall),f1=mean(f1),.groups="drop")
  cat(sprintf("%-12s | %s\n", rg, paste(sprintf("%s F%.2f r%.2f fdr%.3f",s$method,s$f1,s$rec,s$fdr),collapse=" | ")))
}
df<-bind_rows(rows); write.csv(df,file.path(RES,"canonical_bench.csv"),row.names=FALSE)
agg<-df%>%group_by(regime,method)%>%summarise(fdr=mean(fdr),recall=mean(recall),f1=mean(f1),.groups="drop")
agg$regime<-factor(agg$regime,levels=names(regimes)); agg$method<-factor(agg$method,levels=c("fishash","+ depth","+ depth + NB"))
cols<-c(fishash="#1D9E75",`+ depth`="#185FA5",`+ depth + NB`="#D4537E")
th<-theme_bw(base_size=10)+theme(legend.position="top",legend.title=element_blank(),plot.title=element_text(face="bold",size=11),axis.text.x=element_text(angle=20,hjust=1))
p1<-ggplot(agg,aes(regime,recall,colour=method,group=method))+geom_line()+geom_point(size=2.4)+scale_colour_manual(values=cols)+labs(title="Recall (higher better)",x=NULL,y="recall")+th
p2<-ggplot(agg,aes(regime,fdr,colour=method,group=method))+geom_line()+geom_point(size=2.4)+geom_hline(yintercept=Q,linetype="dashed",colour="grey50")+scale_colour_manual(values=cols)+labs(title="Realized FDR (dashed=5%)",x=NULL,y="FDR")+th
p3<-ggplot(agg,aes(regime,f1,colour=method,group=method))+geom_line()+geom_point(size=2.4)+scale_colour_manual(values=cols)+labs(title="F1 (higher better)",x=NULL,y="F1")+th
ggsave(file.path(RES,"canonical_bench.png"),(p1|p3|p2)+plot_layout(guides="collect")&theme(legend.position="top"),width=14,height=4.6,dpi=150)
cat("\nwrote canonical_bench.{csv,png}\n")
