#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Matrix);library(ggplot2);library(grid)})
self <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE))
source(file.path(dirname(normalizePath(self)), 'datasets.R'))
input <- dctap_sceptre_dir("highmoi")
out <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"),
              "association/pos-control/dctap-highmoi-mixture-fishashplus-20261008")
x<-as(readRDS(file.path(input,'grna_matrix_aligned.rds')),'CsparseMatrix')
m<-readRDS(file.path(out,'assignments_mixture.rds'));f<-readRDS(file.path(out,'assignments_fishashplus.rds'))
tdf<-read.csv(file.path(input,'grna_target_data_frame.csv'))
symbols<-setNames(tdf$grna_target_symbol,tdf$grna_target)
joint<-list()
for(scope in c('promoter_union','individual_guides')){
  a<-readRDS(file.path(out,paste0(scope,'_mixture_results.rds')))
  b<-readRDS(file.path(out,paste0(scope,'_fishashplus_results.rds')))
  key<-intersect(c('grna_id','grna_target','response_id'),names(a))
  d<-merge(a,b,by=key,suffixes=c('_mixture','_fishashplus'))
  d$gene<-unname(symbols[d$grna_target])
  d$type<-if(scope=='promoter_union')'Promoter union' else ifelse(grepl('TSS',d$grna_id),'Promoter','Enhancer')
  for(arm in c('mixture','fishashplus')){
    assn<-if(arm=='mixture')m else f
    d[[paste0('n_treated_',arm)]]<-vapply(seq_len(nrow(d)),function(k){
      ids<-if(scope=='individual_guides')d$grna_id[k] else tdf$grna_id[tdf$grna_target==d$grna_target[k]&grepl('TSS',tdf$grna_id)]
      sum(colSums(assn[ids,,drop=FALSE])>0)
    },numeric(1))
  }
  d$delta_log2fc<-d$log_2_fold_change_fishashplus-d$log_2_fold_change_mixture
  write.csv(d,file.path(out,paste0(scope,'_paired.csv')),row.names=FALSE)
  joint[[scope]]<-d
}
capture.output(lapply(joint,function(d)d[,c('gene','type','n_treated_mixture','n_treated_fishashplus',
  'log_2_fold_change_mixture','log_2_fold_change_fishashplus','p_value_mixture','p_value_fishashplus',
  'stage_mixture','stage_fishashplus')]),file=file.path(out,'paired_summary.txt'))
draw_scatter<-function(){
  par(mfrow=c(2,2),mar=c(4.5,4.5,3,1),oma=c(2,0,4,0),family='sans')
  cols<-c(GATA1='#43865E',HDAC6='#387CA1',MYC='#C34E3D')
  for(scope in names(joint)){
    d<-joint[[scope]];c<-cols[d$gene];shape<-ifelse(d$type=='Enhancer',1,16)
    title<-if(scope=='promoter_union')'3 genes: union of promoter guides' else '12 individual guides'
    a<-d$log_2_fold_change_mixture;b<-d$log_2_fold_change_fishashplus;lim<-range(c(a,b));span<-diff(lim)
    lim<-lim+c(-1,1)*span*.1
    plot(a,b,xlim=lim,ylim=lim,pch=shape,col=c,cex=1.3,
      xlab='Mixture: log2 fold change',ylab='Fishash+: log2 fold change',main=paste(title,'— effects'))
    abline(0,1,col='grey50',lty=2)
    if(scope=='promoter_union')text(a,b,d$gene,pos=3,cex=.85,col=c)
    else legend('topleft',legend=c('Promoter','Enhancer'),pch=c(16,1),bty='n',cex=.8)
    a<--log10(d$p_value_mixture);b<--log10(d$p_value_fishashplus);lim<-range(c(0,a,b))
    empirical<-d$stage_fishashplus==3|d$stage_mixture==3
    plot(a,b,xlim=lim+c(0,10),ylim=lim+c(0,10),pch=ifelse(empirical,8,shape),col=c,cex=1.3,
      xlab='Mixture: -log10(p)',ylab='Fishash+: -log10(p)',main=paste(title,'— p-values'))
    abline(0,1,col='grey50',lty=2)
    if(scope=='promoter_union')text(a,b,d$gene,pos=3,cex=.85,col=c)
    if(any(empirical))text(a[empirical],b[empirical],'GATA1-TSS-2: empirical fallback',pos=2,cex=.7)
    if(scope=='individual_guides')legend('topleft',legend=names(cols),col=cols,pch=16,bty='n',cex=.8)
  }
  mtext('High-MOI DC-TAP K562: mixture versus Poisson Fishash+',outer=TRUE,side=3,line=2,cex=1.35,font=2)
  mtext('Same 10,932 cells and fixed expression / guide-count covariates; SCEPTRE CRT',outer=TRUE,side=3,line=.5,cex=.95)
  mtext('Dashed lines = equality. Top-row GATA1 p-values hit 1e-250. Star = empirical fallback (p = 0.0004).',outer=TRUE,side=1,line=.3,cex=.9)
}
png(file.path(out,'positive-control-scatter.png'),width=2000,height=1500,res=160);draw_scatter();dev.off()
pdf(file.path(out,'positive-control-scatter.pdf'),width=12.5,height=9.375);draw_scatter();dev.off()
# Distribution plot uses the 12 individual guide tests, split by guide type.
d<-joint$individual_guides
long<-rbind(data.frame(method='Mixture',type=d$type,log2fc=d$log_2_fold_change_mixture,logp=-log10(d$p_value_mixture)),
            data.frame(method='Fishash+',type=d$type,log2fc=d$log_2_fold_change_fishashplus,logp=-log10(d$p_value_fishashplus)))
long$method<-factor(long$method,c('Mixture','Fishash+'))
violin<-function(var,title,ylab){ggplot(long,aes(method,.data[[var]],fill=method))+
  geom_violin(trim=TRUE,alpha=.4,color=NA)+geom_boxplot(width=.12,outlier.shape=NA,fill='white')+
  geom_point(position=position_jitter(width=.03,height=0,seed=20261008),size=2)+
  facet_wrap(~type)+scale_fill_manual(values=c(Mixture='#43865E','Fishash+'='#387CA1'),guide='none')+
  labs(title=title,x=NULL,y=ylab,subtitle='Six guides per type; two per gene (not independent biological replicates)')+
  theme_bw(base_size=13)+theme(panel.grid.minor=element_blank())}
p1<-violin('log2fc','Individual-guide effect distributions','Log2 fold change')
p2<-violin('logp','Individual-guide p-value distributions','-log10(p)')+
  labs(caption='One Fishash+ GATA1 promoter guide uses empirical fallback p=0.0004; treat its p-value rank cautiously.')
draw_v<-function(){grid.newpage();pushViewport(viewport(layout=grid.layout(1,2)))
  print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1));print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2));popViewport()}
png(file.path(out,'positive-control-violins.png'),width=2200,height=1000,res=160);draw_v();dev.off()
pdf(file.path(out,'positive-control-violins.pdf'),width=13.75,height=6.25);draw_v();dev.off()
# Area-proportional overlap diagram.
nf<-sum(f);nm<-sum(m);sh<-sum(f*m);r<-sqrt(nm/nf)
lens<-function(d){if(d<=1-r)return(pi*r*r);if(d>=1+r)return(0)
  acos((d*d+1-r*r)/(2*d))+r*r*acos((d*d+r*r-1)/(2*d*r))-
    .5*sqrt((-d+1+r)*(d+1-r)*(d-1+r)*(d+1+r))}
dist<-uniroot(function(d)lens(d)-pi*sh/nf,c(1-r+1e-9,1+r-1e-9))$root
draw_venn<-function(){
 par(mar=c(1,1,4,1),family='sans');plot.new();plot.window(xlim=c(-1.9,1.9),ylim=c(-1.3,1.25),asp=1)
 th<-seq(0,2*pi,length.out=1000)
 polygon(cos(th),sin(th),col='#DCECF5',border='#387CA1',lwd=2)
 polygon(-dist+r*cos(th),r*sin(th),col=adjustcolor('#43865E',.35),border='#43865E',lwd=2)
 text(-.12,.16,'Shared',font=2,cex=1.2);text(-.12,-.05,format(sh,big.mark=','),font=2,cex=1.8)
 text(-.12,-.27,sprintf('%.2f%% of mixture retained',100*sh/nm),cex=1)
 segments(-dist-r,0,-1.3,0,col='#43865E');text(-1.34,.1,'Mixture only',adj=1,font=2);text(-1.34,-.12,format(nm-sh,big.mark=','),adj=1,cex=1.3)
 segments(.995,0,1.3,0,col='#387CA1');text(1.34,.1,'Fishash+ only',adj=0,font=2);text(1.34,-.12,format(nf-sh,big.mark=','),adj=0,cex=1.3)
 text(0,1.13,paste0('Mixture: ',format(nm,big.mark=','),'    Fishash+: ',format(nf,big.mark=',')),cex=1.15)
 title('High-MOI DC-TAP K562: call overlap',cex.main=1.4)
 mtext('110 guides × 10,932 cells · each call = one guide–cell pair',side=3,line=.3,cex=.95)
 text(0,-1.16,'Circle areas proportional to calls',cex=.9)
}
png(file.path(out,'call-overlap-venn.png'),width=1700,height=1100,res=160);draw_venn();dev.off()
pdf(file.path(out,'call-overlap-venn.pdf'),width=10.625,height=6.875);draw_venn();dev.off()
# Nonzero per-pair UMI distributions and per-cell called/uncalled guide counts.
tr<-as(x,'TsparseMatrix');v<-tr@x
plots<-list()
for(arm in c('Mixture','Fishash+')){
 a<-if(arm=='Mixture')m else f;called<-as.logical(a[cbind(tr@i+1L,tr@j+1L)])
 bins<-cut(v,c(.5,1.5,2.5,4.5,9.5,19.5,49.5,99.5,199.5,499.5,Inf),
           labels=c('1','2','3-4','5-9','10-19','20-49','50-99','100-199','200-499','500+'))
 tab<-as.data.frame(table(bins,called));names(tab)<-c('UMIs','called','n')
 tab$group<-ifelse(tab$called=='TRUE','Called','Not called');tab$method<-arm
 tab$percent<-ave(tab$n,tab$group,FUN=function(z)100*z/sum(z))
 plots[[arm]]<-tab
}
counts<-do.call(rbind,plots);write.csv(counts,file.path(out,'nonzero_umi_distributions.csv'),row.names=FALSE)
p<-ggplot(counts,aes(UMIs,percent,fill=group))+geom_col(position='dodge')+facet_wrap(~method)+
 scale_fill_manual(values=c(Called='#387CA1','Not called'='#C34E3D'))+
 labs(title='High-MOI DC-TAP: UMI counts in called versus uncalled pairs',
      subtitle='Nonzero guide-cell pairs only; percentages normalized separately within each group',
      y='Percent within group',x='UMIs per guide-cell pair')+
 theme_bw(base_size=13)+theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='bottom')
ggsave(file.path(out,'count-distributions.png'),p,width=13,height=6,dpi=160)
ggsave(file.path(out,'count-distributions.pdf'),p,width=13,height=6)
cat('Plotting complete\n');print(joint$promoter_union[,c('gene','n_treated_mixture','n_treated_fishashplus','delta_log2fc')])
