#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(Matrix))
base <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
out <- file.path(base,'outputs/gasperini/fishashplus-poisson-20261008')
mix_fp <- file.path(base,'outputs/gasperini/sceptre-pipeline/grna_assignment_matrix.rds')
raw <- readRDS(mix_fp)
x <- as(readRDS(file.path(base,'input_data/gasperini/sceptre/grna_matrix.rds')),'CsparseMatrix')
f <- readRDS(file.path(out,'grna_assignment_matrix.rds'))
stopifnot(!anyDuplicated(rownames(raw)),setequal(rownames(raw),rownames(x)),
          ncol(raw)==ncol(x),identical(dimnames(x),dimnames(f)))
row_map<-match(rownames(x),rownames(raw))
m<-as(raw[row_map,],'CsparseMatrix')
# Pipeline get_grna_assignments exports original cell indices without names.
# The saved Nextflow command points to this full Gasperini ODM and stops before QC.
colnames(m)<-colnames(x)
stopifnot(identical(dimnames(m),dimnames(x)),sum(m*(x>0))==sum(m))
saveRDS(m,file.path(out,'mixture_assignment_matrix_aligned.rds'))
write.csv(data.frame(current_row=seq_along(row_map),saved_mixture_row=row_map,
                     guide=rownames(x)),file.path(out,'mixture_guide_row_alignment.csv'),row.names=FALSE)
shared<-sum(m*f);nm<-sum(m);nf<-sum(f)
stats<-list(shared=shared,mixture_only=nm-shared,fishashplus_only=nf-shared,
            mixture_total=nm,fishashplus_total=nf,
            fraction_mixture_retained=shared/nm,
            fraction_fishashplus_shared=shared/nf,
            source=mix_fp,source_md5=tools::md5sum(mix_fp),
            cell_alignment='Original ODM column indices retained by Nextflow; stop before QC; every mixture call verified against nonzero current input',
            rows_reordered=sum(row_map!=seq_along(row_map)))
saveRDS(stats,file.path(out,'mixture_overlap_metadata.rds'))
capture.output(print(stats),file=file.path(out,'mixture_overlap_summary.txt'))
# Circle areas and intersection area proportional to numbers of calls.
rf<-1;rm<-sqrt(nm/nf)
lens<-function(d){
  if(d<=abs(rf-rm))return(pi*min(rf,rm)^2)
  if(d>=rf+rm)return(0)
  rf^2*acos((d^2+rf^2-rm^2)/(2*d*rf))+
    rm^2*acos((d^2+rm^2-rf^2)/(2*d*rm))-
    .5*sqrt((-d+rf+rm)*(d+rf-rm)*(d-rf+rm)*(d+rf+rm))
}
distance<-uniroot(function(d)lens(d)-pi*shared/nf,
                  c(abs(rf-rm)+1e-9,rf+rm-1e-9))$root
draw<-function(){
  par(mar=c(1,1,5,1),family='sans',bg='white')
  plot.new();plot.window(xlim=c(-1.9,1.9),ylim=c(-1.35,1.25),asp=1)
  th<-seq(0,2*pi,length.out=1000)
  polygon(rf*cos(th),rf*sin(th),col='#DCECF5',border='#387CA1',lwd=2.5)
  polygon(-distance+rm*cos(th),rm*sin(th),
          col=adjustcolor('#C34E3D',.28),border='#C34E3D',lwd=2.5)
  text(0,.18,'Shared',cex=1.2,font=2,col='#22313B')
  text(0,-.06,format(shared,big.mark=','),cex=1.65,font=2,col='#22313B')
  text(0,-.29,sprintf('%.2f%% of mixture calls retained',100*shared/nm),cex=1.1,col='#22313B')
  segments(-distance-rm+.002,0,-1.28,0,col='#C34E3D',lwd=1.5)
  text(-1.32,.11,'Mixture only',adj=1,cex=1.15,font=2,col='#22313B')
  text(-1.32,-.14,format(nm-shared,big.mark=','),adj=1,cex=1.4,col='#22313B')
  segments(.995,0,1.28,0,col='#387CA1',lwd=1.5)
  text(1.32,.11,'Fishash+ only',adj=0,cex=1.15,font=2,col='#22313B')
  text(1.32,-.14,format(nf-shared,big.mark=','),adj=0,cex=1.4,col='#22313B')
  text(0,1.15,paste0('Mixture: ',format(nm,big.mark=','),'    Fishash+: ',format(nf,big.mark=',')),
       cex=1.1,col='#22313B')
  text(0,-1.19,'Circle areas proportional to calls; mixture-only region is a very thin sliver.',cex=.92,col='#4C5B66')
  title(main='Full Gasperini: mixture versus Fishash+',cex.main=1.45,col.main='#22313B')
  mtext('13,077 guides × 207,324 cells · rows aligned by guide identity; original cell indices preserved',
        side=3,line=.5,cex=.94,col='#4C5B66')
}
png(file.path(out,'fishashplus-vs-mixture-venn.png'),width=1700,height=1100,res=160);draw();dev.off()
pdf(file.path(out,'fishashplus-vs-mixture-venn.pdf'),width=10.625,height=6.875);draw();dev.off()
print(stats)
