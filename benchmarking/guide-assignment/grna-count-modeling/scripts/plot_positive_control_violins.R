#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(grid)})
out <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"),
              "association/pos-control/fishashplus-comparison-20261008")
d <- readRDS(file.path(out,'paired_positive_controls.rds'))
df <- data.frame(method=factor(rep(c('Count >= 5','Fishash+'),each=nrow(d)),
                              levels=c('Count >= 5','Fishash+')),
                 log2fc=c(d$log_2_fold_change_threshold5,d$log_2_fold_change_fishashplus),
                 logp=c(-log10(d$p_value_threshold5),-log10(d$p_value_fishashplus)))
include_mixture <- 'mixture' %in% commandArgs(trailingOnly=TRUE)
suffix <- if(include_mixture)'-three-methods' else ''
if(include_mixture){
  m<-readRDS(file.path(out,'full_mixture_results.rds'))
  df$method<-as.character(df$method)
  df<-rbind(df,data.frame(method='Mixture',log2fc=m$log_2_fold_change,
                          logp=-log10(m$p_value)))
  df$method<-factor(df$method,levels=c('Count >= 5','Mixture','Fishash+'))
}
stopifnot(all(is.finite(df$log2fc)),all(is.finite(df$logp)))
col <- c('Count >= 5'='#C34E3D','Fishash+'='#387CA1','Mixture'='#43865E')
base <- function(value,title,ylab,subtitle){
  ggplot(df,aes(x=method,y=.data[[value]],fill=method))+
    geom_violin(trim=TRUE,scale='width',alpha=.45,color=NA,width=.85)+
    geom_boxplot(width=.14,fill='white',color='#22313B',outlier.shape=NA,linewidth=.55)+
    geom_point(position=position_jitter(width=.035,height=0,seed=20261008),
               size=.7,alpha=.27,color='#22313B')+
    stat_summary(fun=median,geom='point',shape=21,fill='white',color='#22313B',size=2.5)+
    scale_fill_manual(values=col,guide='none')+
    labs(title=title,subtitle=subtitle,x=NULL,y=ylab)+
    theme_bw(base_size=15)+theme(panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),
      plot.subtitle=element_text(size=11),plot.margin=margin(10,14,10,14))
}
fcmed <- tapply(df$log2fc,df$method,median)
pmed <- tapply(df$logp,df$method,median)
fc <- base('log2fc','Log2 fold changes','SCEPTRE log2 fold change',
  paste('Medians:',paste(sprintf('%s %.3f',names(fcmed),fcmed),collapse='; ')))+
  geom_hline(yintercept=0,color='grey60',linetype=2)
pv <- base('logp','P-values on a log scale','-log10(SCEPTRE CRT p-value)',
  paste('Medians:',paste(sprintf('%s %.1f',names(pmed),pmed),collapse='; ')))+
  geom_hline(yintercept=250,color='grey50',linetype=2)+
  labs(caption='p = 1e-250 floor: 71 threshold controls; 62 Fishash+ controls')+
  theme(plot.caption=element_text(size=10,hjust=0))
if(include_mixture)pv<-pv+labs(caption=paste0('p = 1e-250 floor: threshold 71; mixture ',
                                      sum(m$p_value<=1.01e-250),'; Fishash+ 62'))
draw <- function(){
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,2,heights=unit(c(.10,.9),'npc'))))
  grid.text('Gasperini positive controls: assignment distributions',
            x=.5,y=.68,gp=gpar(fontsize=20,fontface='bold'),
            vp=viewport(layout.pos.row=1,layout.pos.col=1:2))
  grid.text('Same 377 target-gene pairs and 207,324 cells in each method; boxes show median and interquartile range',
            x=.5,y=.18,gp=gpar(fontsize=12),
            vp=viewport(layout.pos.row=1,layout.pos.col=1:2))
  print(fc,vp=viewport(layout.pos.row=2,layout.pos.col=1))
  print(pv,vp=viewport(layout.pos.row=2,layout.pos.col=2))
  popViewport()
}
png(file.path(out,paste0('positive-control-violins',suffix,'.png')),width=2100,height=1200,res=160);draw();dev.off()
pdf(file.path(out,paste0('positive-control-violins',suffix,'.pdf')),width=13.125,height=7.5);draw();dev.off()
cat('Log2FC medians:',fcmed,'\n-log10p medians:',pmed,'\n')
