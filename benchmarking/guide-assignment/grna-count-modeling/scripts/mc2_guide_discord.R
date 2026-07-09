#!/usr/bin/env Rscript
# Generalized: for dataset <ds>, show most-discordant (top) vs most-concordant (bottom) guides
# by per-guide CLEANSER-vs-fishash+ Jaccard, with per-guide count histograms + ambient calls.
suppressMessages({library(Matrix); library(ggplot2); library(dplyr)})
source("scripts/datasets.R")
ds <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(ds)) stop("usage: guide_discord.R <dataset>")
ADIR <- "results/ambient_ceiling/writeup/assign"
gm <- load_grna_matrix(ds)
Af <- as(readRDS(file.path(ADIR,paste0(ds,"_fishash.rds"))),    "lgCMatrix")
Ap <- as(readRDS(file.path(ADIR,paste0(ds,"_fishashplus.rds"))),"lgCMatrix")
cf <- c(file.path("results/ambient_ceiling/cleanser_calls_cs",paste0(ds,".csv")),   # prefer chemistry-correct
        file.path("results/ambient_ceiling/cleanser_calls",   paste0(ds,".csv")))
cf <- cf[file.exists(cf)][1]; chem <- if (grepl("cleanser_calls_cs", cf)) "--cs" else "--dc"
cl <- read.csv(cf); gg<-match(cl$guide,rownames(gm)); cc<-match(cl$cell,colnames(gm)); ok<-!is.na(gg)&!is.na(cc)
Ac <- as(sparseMatrix(i=gg[ok], j=cc[ok], x=TRUE, dims=dim(gm), dimnames=dimnames(gm)), "lgCMatrix")

G<-nrow(gm); ct<-as(gm,"TsparseMatrix"); gi<-ct@i+1L; cj<-ct@j+1L; xv<-ct@x; idx<-cbind(gi,cj)
nz<-tabulate(gi,nbins=G); inter<-rowSums(Ac&Ap); uni<-rowSums(Ac|Ap); jcp<-inter/uni
nAc<-rowSums(Ac); nAp<-rowSums(Ap)
for(thr in c(50,30,20,10,1)){elig<-which(nz>=thr & uni>=5); if(length(elig)>=6)break}
ord<-elig[order(jcp[elig])]; sel<-c(head(ord,3),tail(ord,3)); grp<-c(rep("discordant",3),rep("concordant",3))
cat(sprintf("%s (CLEANSER %s): %d guides x %d cells; median CLEANSER-f+ Jaccard=%.2f (over %d guides)\n",
    ds, chem, G, ncol(gm), median(jcp[elig]), length(elig)))
print(data.frame(guide=rownames(gm)[sel],grp,J=round(jcp[sel],2),n_nonzero=nz[sel],
      n_assigned_cleanser=nAc[sel],n_assigned_fishashplus=nAp[sel]),row.names=FALSE)

amb<-list(cleanser=!Ac[idx],fishash=!Af[idx],`fishash+`=!Ap[idx]); bars<-list(); lines<-list()
for(r in seq_along(sel)){g<-sel[r];s<-gi==g;vals<-xv[s];ambg<-lapply(amb,function(a)a[s])
  maxv<-max(vals);nb<-max(6,min(22,ceiling(12*log10(maxv+1))))
  edges<-seq(log10(0.9),log10(maxv*1.12),length.out=nb+1);bin<-cut(log10(vals),edges,labels=FALSE,include.lowest=TRUE)
  tot<-tabulate(bin,nbins=nb);ctr<-10^((edges[1:nb]+edges[2:(nb+1)])/2)
  lab<-sprintf("%s: %s\n%s  (J=%.2f, n=%d)",ds,grp[r],rownames(gm)[g],jcp[g],nz[g])
  bars[[length(bars)+1]]<-data.frame(panel=lab,ord=r,xmin=10^edges[1:nb],xmax=10^edges[2:(nb+1)],total=tot)
  for(m in names(ambg)){na<-tabulate(bin[ambg[[m]]],nbins=nb);keep<-tot>0
    lines[[length(lines)+1]]<-data.frame(panel=lab,ord=r,x=ctr[keep],method=m,n=na[keep])}}
bars_df<-bind_rows(bars);lines_df<-bind_rows(lines)
lev<-bars_df%>%distinct(panel,ord)%>%arrange(ord)%>%pull(panel)
bars_df$panel<-factor(bars_df$panel,levels=lev);lines_df$panel<-factor(lines_df$panel,levels=lev)
lines_df$method<-factor(lines_df$method,levels=c("cleanser","fishash","fishash+"))
pal<-c("cleanser"="#e34948","fishash"="#eda100","fishash+"="#2a78d6");SHP<-c("cleanser"=17,"fishash"=15,"fishash+"=16)
DODGE<-c("cleanser"=1/1.07,"fishash"=1,"fishash+"=1.07);FLOOR<-0.5;bars_df<-subset(bars_df,total>0)
lines_df$n_plot<-pmax(lines_df$n,FLOOR);lines_df$xd<-lines_df$x*DODGE[as.character(lines_df$method)]
p<-ggplot()+
  geom_rect(data=bars_df,aes(xmin=xmin,xmax=xmax,ymin=FLOOR,ymax=total),fill="grey88",color="grey74",linewidth=0.12)+
  geom_line(data=lines_df,aes(xd,n_plot,color=method),linewidth=0.5,alpha=0.9)+
  geom_point(data=lines_df,aes(xd,n_plot,color=method,shape=method),size=1.4,alpha=0.95)+
  facet_wrap(~panel,scales="free",ncol=3)+
  scale_x_log10()+scale_y_log10()+scale_color_manual(values=pal)+scale_shape_manual(values=SHP)+
  labs(x="guide UMI count across cells (log scale; lines slightly x-offset)",y="number of cells",
       color="ambient calls",shape="ambient calls",
       title=sprintf("%s (CLEANSER %s): CLEANSER vs fishash+ discordant (top) and concordant (bottom) guides",ds,chem))+
  theme_bw(base_size=10)+theme(legend.position="bottom",strip.text=element_text(size=8),plot.title=element_text(size=11))
ggsave(sprintf("results/ambient_ceiling/%s_discordant_hist.png",ds),p,width=11,height=6.6,dpi=150)
cat(sprintf("\nwrote results/ambient_ceiling/%s_discordant_hist.png\n",ds))
