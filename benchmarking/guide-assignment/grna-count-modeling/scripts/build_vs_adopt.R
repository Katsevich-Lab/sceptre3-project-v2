#!/usr/bin/env Rscript
# Build-vs-adopt evidence: across Simpson x overdispersion regimes, compare the plain ambient
# (Poisson) test, two variants of our decontaminated NB candidate (global vs robust phi), and
# fishash. Shows the candidate is bracketed by two failure modes (FDR vs power) while fishash
# attains both. Produces results/benchmark_update/build_vs_adopt.{csv,png}.
suppressPackageStartupMessages({library(fishash); library(Matrix); library(SummarizedExperiment)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})
HERE <- tryCatch(dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)))),error=function(e) ".")
HERE <- dirname(HERE); source(file.path(HERE,"scripts","ambient_nb_method.R"))
RES <- file.path(HERE,"results","benchmark_update"); dir.create(RES,showWarnings=FALSE,recursive=TRUE)
score <- function(A,truth){A<-as(A,"lgCMatrix");truth<-as(truth,"lgCMatrix");na<-sum(A);nt<-sum(truth);tp<-sum(A&truth)
  c(fdr=if(na>0)(na-tp)/na else 0, recall=if(nt>0)tp/nt else NA_real_)}
regimes <- list(clean=c(0,0,4,100), simpson=c(1,0,4,100), overdisp=c(0,1,4,100),
                both=c(1,1,4,100), hard=c(0.5,0,1,20))
NREP<-4; Q<-0.05; rows<-list()
for (rg in names(regimes)){p<-regimes[[rg]]
  for (rep in seq_len(NREP)){ set.seed(100*which(names(regimes)==rg)+rep)
    sim<-simulate_guidebender2(n_guides=200,n_cells=4000,moi=0.3,hurdle_prob=0.1,snr=p[3],
          count_per_cell=p[4],frac_noise_endo=0.75,endo_shape_flat=p[1],Phi_noise=p[2],return_sparse_only=TRUE)
    N<-assay(sim,"counts"); tr<-assay(sim,"ground_truth")
    m<-list(
      `ambient (Poisson)`        = score(ambient_nb_assign(N,q=Q,n_iter=1,overdispersion=FALSE)$assigned,tr),
      `NB candidate (phi global)`= score(ambient_nb_assign(N,q=Q,n_iter=3,overdispersion=TRUE,phi_method="global")$assigned,tr),
      `NB candidate (phi robust)`= score(ambient_nb_assign(N,q=Q,n_iter=3,overdispersion=TRUE,phi_method="robust")$assigned,tr),
      `fishash`                  = score(assay(fishash(N,padj_cutoff=Q),"assigned"),tr))
    for (mn in names(m)) rows[[length(rows)+1]]<-data.frame(regime=rg,rep=rep,method=mn,fdr=m[[mn]][1],recall=m[[mn]][2])
  }
  cat("done", rg, "\n")
}
df<-do.call(rbind,rows); write.csv(df,file.path(RES,"build_vs_adopt.csv"),row.names=FALSE)
agg<-df%>%group_by(regime,method)%>%summarise(fdr=mean(fdr),recall=mean(recall),.groups="drop")
agg$regime<-factor(agg$regime,levels=names(regimes))
agg$method<-factor(agg$method,levels=c("ambient (Poisson)","NB candidate (phi global)","NB candidate (phi robust)","fishash"))
cols<-c(`ambient (Poisson)`="#888780",`NB candidate (phi global)`="#D4537E",`NB candidate (phi robust)`="#D85A30",`fishash`="#1D9E75")
th<-theme_bw(base_size=11)+theme(legend.position="top",legend.title=element_blank(),plot.title=element_text(face="bold",size=11))
p1<-ggplot(agg,aes(regime,fdr,fill=method))+geom_col(position=position_dodge(0.8),width=0.75)+
  geom_hline(yintercept=Q,linetype="dashed",colour="grey40")+scale_fill_manual(values=cols)+
  labs(title="Realized FDR by regime (dashed = nominal 5%)",y="FDR",x=NULL)+th
p2<-ggplot(agg,aes(regime,recall,fill=method))+geom_col(position=position_dodge(0.8),width=0.75)+
  scale_fill_manual(values=cols)+labs(title="Recall by regime",y="recall",x=NULL)+th
ggsave(file.path(RES,"build_vs_adopt.png"),p1/p2+plot_layout(guides="collect")&theme(legend.position="top"),
       width=10,height=7,dpi=150)
cat("\nwrote build_vs_adopt.{csv,png}\n")
agg%>%mutate(fdr=round(fdr,3),recall=round(recall,2))%>%pivot_wider(names_from=method,values_from=c(fdr,recall))%>%as.data.frame()%>%print()
