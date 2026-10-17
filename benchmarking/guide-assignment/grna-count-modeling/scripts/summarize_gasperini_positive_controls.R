#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Matrix))
out <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"),
              "association/pos-control/fishashplus-comparison-20261008")
t <- as.data.frame(readRDS(file.path(out,'full_threshold5_results.rds')))
f <- as.data.frame(readRDS(file.path(out,'full_fishashplus_results.rds')))
keys<-c('response_id','grna_target')
stopifnot(nrow(t)==377L,nrow(f)==377L,!anyDuplicated(t[keys]),!anyDuplicated(f[keys]))
t$bh_adjusted_p<-p.adjust(t$p_value,'BH');f$bh_adjusted_p<-p.adjust(f$p_value,'BH')
d<-merge(t,f,by=keys,suffixes=c('_threshold5','_fishashplus'),sort=FALSE)
stopifnot(nrow(d)==377L)
dat<-readRDS(file.path(out,'paired_inputs.rds'))
for (arm in c('threshold5','fishashplus')) {
  treated <- vapply(d$grna_target,function(g){
    idx<-which(dat$targets$grna_target==g)
    sum(colSums(dat[[arm]][idx,,drop=FALSE])>0)
  },numeric(1))
  d[[paste0('n_treated_',arm)]]<-treated
}
d$additional_treated_cells<-d$n_treated_fishashplus-d$n_treated_threshold5
d$fraction_added_treated<-d$additional_treated_cells/d$n_treated_threshold5
d$delta_log2_fold_change<-d$log_2_fold_change_fishashplus-d$log_2_fold_change_threshold5
d$delta_fold_change<-d$fold_change_fishashplus-d$fold_change_threshold5
d$stronger_knockdown<-d$log_2_fold_change_fishashplus<d$log_2_fold_change_threshold5
d$smaller_p_value<-d$p_value_fishashplus<d$p_value_threshold5
d$detected_threshold5<-d$bh_adjusted_p_threshold5<=.1
d$detected_fishashplus<-d$bh_adjusted_p_fishashplus<=.1
d$p_at_floor_threshold5<-d$p_value_threshold5<=1.01e-250
d$p_at_floor_fishashplus<-d$p_value_fishashplus<=1.01e-250
write.csv(d,file.path(out,'paired_positive_controls.csv'),row.names=FALSE)
saveRDS(d,file.path(out,'paired_positive_controls.rds'))
valid<-is.finite(d$log_2_fold_change_threshold5)&is.finite(d$log_2_fold_change_fishashplus)&
  is.finite(d$p_value_threshold5)&is.finite(d$p_value_fishashplus)
summary<-list(n_pairs=nrow(d),valid_pairs=sum(valid),
  fold_change_spearman=cor(d$log_2_fold_change_threshold5,d$log_2_fold_change_fishashplus,method='spearman',use='complete.obs'),
  median_log2fc_threshold5=median(d$log_2_fold_change_threshold5,na.rm=TRUE),
  median_log2fc_fishashplus=median(d$log_2_fold_change_fishashplus,na.rm=TRUE),
  median_delta_log2fc=median(d$delta_log2_fold_change,na.rm=TRUE),
  stronger_knockdowns=sum(d$stronger_knockdown,na.rm=TRUE),
  smaller_p_values=sum(d$smaller_p_value,na.rm=TRUE),
  bh10_detected_threshold5=sum(d$detected_threshold5,na.rm=TRUE),
  bh10_detected_fishashplus=sum(d$detected_fishashplus,na.rm=TRUE),
  gained_bh10=sum(d$detected_fishashplus & !d$detected_threshold5,na.rm=TRUE),
  lost_bh10=sum(!d$detected_fishashplus & d$detected_threshold5,na.rm=TRUE),
  p_zero_threshold5=sum(d$p_value_threshold5==0,na.rm=TRUE),
  p_zero_fishashplus=sum(d$p_value_fishashplus==0,na.rm=TRUE),
  p_at_floor_threshold5=sum(d$p_at_floor_threshold5),
  p_at_floor_fishashplus=sum(d$p_at_floor_fishashplus),
  p_at_floor_both=sum(d$p_at_floor_threshold5 & d$p_at_floor_fishashplus),
  median_treated_threshold5=median(d$n_treated_threshold5),
  median_treated_fishashplus=median(d$n_treated_fishashplus),
  median_percent_added_treated=100*median(d$fraction_added_treated))
print(summary);saveRDS(summary,file.path(out,'comparison_summary.rds'))
capture.output(print(summary),file=file.path(out,'comparison_summary.txt'))
floorp<-.Machine$double.xmin
xp<--log10(pmax(d$p_value_threshold5,floorp));yp<--log10(pmax(d$p_value_fishashplus,floorp))
draw<-function(){
  par(mfrow=c(1,3),mar=c(5,5,3,1),oma=c(1,0,4,0),family='sans',bg='white')
  lim<-range(c(d$log_2_fold_change_threshold5,d$log_2_fold_change_fishashplus),finite=TRUE)
  plot(d$log_2_fold_change_threshold5,d$log_2_fold_change_fishashplus,
    pch=16,cex=.65,col=adjustcolor('#387CA1',.7),xlim=lim,ylim=lim,
    xlab='Count >= 5: log2 fold change',ylab='Fishash+: log2 fold change',main='Effect estimates')
  abline(0,1,col='grey50',lty=2)
  abline(h=0,v=0,col='grey85')
  text(lim[1],lim[2],paste0('Spearman r = ',round(summary$fold_change_spearman,3)),adj=c(0,1),cex=.85)
  lim<-range(c(xp,yp),finite=TRUE)
  plot(xp,yp,pch=ifelse(d$p_at_floor_threshold5|d$p_at_floor_fishashplus,17,16),
    cex=.65,col=adjustcolor('#387CA1',.7),xlim=lim,ylim=lim,
    xlab='Count >= 5: -log10(p)',ylab='Fishash+: -log10(p)',main='SCEPTRE CRT p-values')
  abline(0,1,col='grey50',lty=2)
  plot(100*d$fraction_added_treated,d$delta_log2_fold_change,pch=16,cex=.65,
    col=adjustcolor('#387CA1',.7),xlab='Additional treated cells (%)',
    ylab='Fishash+ minus threshold: log2 FC',main='Change versus added calls')
  abline(h=0,col='grey50',lty=2)
  mtext('Gasperini positive controls: paired assignment comparison',outer=TRUE,side=3,line=2,cex=1.4,font=2)
  mtext('377 target-gene pairs · union of two guides per target · same 207,324 cells and full-library covariates',
        outer=TRUE,side=3,line=.4,cex=.94)
  if(summary$p_at_floor_threshold5+summary$p_at_floor_fishashplus>0)
    mtext('Triangles indicate at least one p-value at the SCEPTRE numerical floor of 1e-250.',outer=TRUE,side=1,line=0,cex=.85)
}
png(file.path(out,'positive-control-comparison.png'),width=2400,height=1000,res=160);draw();dev.off()
pdf(file.path(out,'positive-control-comparison.pdf'),width=15,height=6.25);draw();dev.off()
