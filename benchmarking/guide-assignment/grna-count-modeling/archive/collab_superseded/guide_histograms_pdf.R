## Scrollable multi-page PDF of the count distribution of each of the 109 power-positive clean-gap
## Replogle guides (the aggregate cohort), so we can eyeball low-mode structure. Per guide: number of
## cells at each positive gRNA UMI count, on log-log axes; points coloured by mode (low = count<=gap_lo,
## high = count>=gap_hi); the empty gap shaded; title = guide + [gap_lo -> gap_hi].
suppressMessages({library(Matrix); library(ondisc); library(sceptre); library(ggplot2)})
OUT_SRC<-"results/global_ambient_poisson"; OUT<-"results/collaborator_writeup"
Dm<-path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp<-path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")
mc<-as(readRDS(Dm),"CsparseMatrix"); resp<-initialize_odm_from_backing_file(file.path(Dp,"response.odm"))
so<-readRDS(file.path(Dp,"sceptre_object.rds")); lib<-exp(so@covariate_matrix[,"log(response_n_umis)"])
tdf<-so@grna_target_data_frame; gids<-rownames(resp)
ntg<-tdf$grna_id[tdf$grna_target=="non-targeting"]; ntmask<-rep(FALSE,ncol(mc)); for(g in ntg) ntmask<-ntmask|(as.numeric(mc[g,])>=30)
pg<-read.csv(file.path(OUT_SRC,"perguide_replogle_rd7.csv"),stringsAsFactors=FALSE); pg<-pg[pg$gap_lo>=2,]
pg$target<-tdf$grna_target[match(pg$guide,tdf$grna_id)]; pg<-pg[!is.na(pg$target)&pg$target!="non-targeting"&pg$target%in%gids,]

pts<-list(); rects<-list(); k<-0
for(i in seq_len(nrow(pg))){ r<-match(pg$guide[i],rownames(mc)); v<-as.numeric(mc[r,]); cp<-as.numeric(resp[pg$target[i],])/lib*1e4
  base<-mean(cp[ntmask&v==0]); gh<-pg$gap_hi[i]; lo<-pg$gap_lo[i]; st<-which(v>=gh)
  if(base<0.5||!length(st)||mean(cp[st])/base>0.7) next
  k<-k+1
  tab<-table(v[v>0]); cnt<-as.integer(names(tab)); nc<-as.integer(tab)
  lab<-sprintf("%s  [%d→%d]", sub("_P1P2.*|_P1_.*","",pg$guide[i]), lo, gh)
  mode<-ifelse(cnt<=lo,"low mode (≤ gap_lo)", ifelse(cnt>=gh,"high mode (≥ gap_hi)","gap"))
  pts[[k]]<-data.frame(guide=lab, gap_lo=lo, count=cnt, ncells=nc, mode=mode, stringsAsFactors=FALSE)
  rects[[k]]<-data.frame(guide=lab, gap_lo=lo, xmin=lo, xmax=gh, stringsAsFactors=FALSE)
}
P<-do.call(rbind,pts); R<-do.call(rbind,rects)
ord<-unique(P[order(P$gap_lo,P$guide),"guide"]); P$guide<-factor(P$guide,levels=ord); R$guide<-factor(R$guide,levels=ord)
cat(sprintf("guides: %d\n", length(ord)))

PER<-20; pages<-split(ord, ceiling(seq_along(ord)/PER))
pdf(file.path(OUT,"guide_histograms_109.pdf"), width=14, height=9)
for(pg_lv in pages){
  Pi<-P[P$guide %in% pg_lv,]; Ri<-R[R$guide %in% pg_lv,]
  Pi$guide<-factor(Pi$guide,levels=pg_lv); Ri$guide<-factor(Ri$guide,levels=pg_lv)
  print(ggplot(Pi, aes(count, ncells)) +
    geom_rect(data=Ri, aes(xmin=xmin, xmax=xmax, ymin=0.5, ymax=Inf), inherit.aes=FALSE, fill="grey88") +
    geom_segment(aes(xend=count, yend=0.5, colour=mode), linewidth=0.35) +
    geom_point(aes(colour=mode), size=0.9) +
    facet_wrap(~guide, ncol=5, scales="free") +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values=c("low mode (≤ gap_lo)"="#0072B2","high mode (≥ gap_hi)"="#D55E00","gap"="grey60")) +
    labs(x="gRNA UMI count (log)", y="number of cells (log)", colour=NULL,
         title="Replogle power-positive clean-gap guides: per-guide count distribution (grey band = empty gap)") +
    theme_bw(base_size=9) + theme(legend.position="top", strip.text=element_text(size=7.5)))
}
dev.off()
cat("wrote", file.path(OUT,"guide_histograms_109.pdf"), "\n")
