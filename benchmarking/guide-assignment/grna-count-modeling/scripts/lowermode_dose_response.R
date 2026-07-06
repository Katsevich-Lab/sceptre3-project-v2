## ============================================================================
## What are the below-signal ("lower-mode") cells, really? Doublet, artifact, or
## real low-dose integration? And does the answer generalize?
##  (A) RNA-content doublet check: 2-strong-guide barcodes are ~1x single-cell RNA
##      (1.05x), NOT ~2x -> they are single cells (co-integration/MOI), not doublets.
##  (B) Dose-response (Replogle): target knockdown deepens smoothly with the
##      below-gap UMI count -> real low-dose integrations (count ~ dose). The count
##      does NOT track total gRNA (refutes the chimera/artifact hypothesis).
##  (C) within-guide paired test (controls for composition) + per-guide panels.
##  (D) cross-dataset: the dose-response holds on Replogle, EndoC, CD4 T-cell.
## Together with the barnyard Poisson result (ambient_null_poisson_vs_nb.R), this
## says the count-2 "excess" is weak-integration SIGNAL, not overdispersed soup.
## Run from grna-count-modeling/.  (cross-dataset tcell step is slow, ~7 min.)
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats); library(ondisc); library(sceptre)})
source("scripts/datasets.R")   # load_replogle_rd7_de()
OUT<-"results/global_ambient_poisson"
rd<-load_replogle_rd7_de(); mc<-rd$mc; totgrna<-as.numeric(colSums(mc))
resp<-rd$resp; so<-rd$so
lib<-exp(so@covariate_matrix[,"log(response_n_umis)"]); tdf<-so@grna_target_data_frame; gids<-rownames(resp)
cdf<-so@covariate_data_frame; umi<-cdf$response_n_umis
ntg<-tdf$grna_id[tdf$grna_target=="non-targeting"]; ntmask<-rep(FALSE,ncol(mc)); for(g in ntg) ntmask<-ntmask|(as.numeric(mc[g,])>=30)
nstrong<-as.numeric(colSums(mc>=30))

## (A) RNA-content doublet check ----------------------------------------------
r1<-median(umi[nstrong==1])
cat("(A) RNA content by # strong guides (doublet => ~2x):\n")
for(k in 1:4){ sel<-if(k<4) nstrong==k else nstrong>=4; cat(sprintf("    %d strong: n=%d med_UMIs=%.0f (%.2fx single)\n",k,sum(sel),median(umi[sel]),median(umi[sel])/r1)) }

pg<-read.csv(file.path(OUT,"perguide_replogle_rd7.csv"),stringsAsFactors=FALSE); pg<-pg[pg$gap_lo>=2,]
pg$target<-tdf$grna_target[match(pg$guide,tdf$grna_id)]; pg<-pg[!is.na(pg$target)&pg$target!="non-targeting"&pg$target%in%gids,]
gene<-function(g) sub("^[0-9]+_([^_]+)_.*","\\1",g)

## (B)+(C) dose-response + within-guide, per power-positive guide --------------
cnt<-c();rat<-c();tg<-c();sigr<-c(); lo<-c();hi<-c(); info<-list()
for(i in seq_len(nrow(pg))){ r<-match(pg$guide[i],rownames(mc)); v<-as.numeric(mc[r,]); cp<-as.numeric(resp[pg$target[i],])/lib*1e4
  base<-mean(cp[ntmask&v==0]); gh<-pg$gap_hi[i]; st<-which(v>=gh); if(base<0.5||!length(st)||mean(cp[st])/base>0.7) next
  bg<-which(v>=1&v<=pg$gap_lo[i]); cnt<-c(cnt,v[bg]);rat<-c(rat,cp[bg]/base);tg<-c(tg,totgrna[bg]);sigr<-c(sigr,cp[st]/base)
  L<-which(v>=2&v<=3); H<-which(v>=5&v<=pg$gap_lo[i]); if(length(L)>=15&&length(H)>=15){lo<-c(lo,mean(cp[L])/base);hi<-c(hi,mean(cp[H])/base)}
  if(length(bg)>=60) info[[pg$guide[i]]]<-list(v=v[bg],r=cp[bg]/base,sigr=mean(cp[st])/base,glo=pg$gap_lo[i],ghi=gh,nbg=length(bg)) }
brk<-list(`1`=1,`2`=2,`3`=3,`4`=4,`5`=5,`6-7`=6:7,`8-11`=8:11,`12-19`=12:19,`20-40`=20:40); xs<-c();ys<-c()
cat("\n(B) Replogle dose-response (target expr/baseline by below-gap count):\n")
for(nm in names(brk)){s<-cnt%in%brk[[nm]];if(sum(s)<20)next;xs<-c(xs,mean(cnt[s]));ys<-c(ys,mean(rat[s]));cat(sprintf("    count %-6s n=%6d  expr/base=%.3f\n",nm,sum(s),mean(rat[s])))}
cat(sprintf("    signal mode: %.3f\n",mean(sigr)))
cat(sprintf("    cor(below-gap count>=2, total gRNA) = %+.3f (artifact would be POSITIVE; it isn't)\n",cor(cnt[cnt>=2],tg[cnt>=2])))
cat(sprintf("(C) within-guide paired (count 2-3 vs >=5), %d guides: hi<lo in %d, median drop %.3f, t p=%.3f\n",
    length(lo),sum(hi<lo),median(lo-hi),t.test(lo,hi,paired=TRUE)$p.value))

png(file.path(OUT,"dose_response.png"),width=1050,height=680,res=140)
par(mar=c(4.6,4.8,3.6,1),mgp=c(2.9,0.8,0),cex.axis=1.0,cex.lab=1.1)
plot(xs,ys,log="x",type="b",pch=16,col="firebrick",cex=1.4,ylim=c(0,1.08),xlim=c(1,max(xs)*60),
     xlab="guide UMI count (below-gap lower mode)",ylab="target expression / NT baseline",
     main="Dose-response: knockdown deepens with the below-gap count\n(=> real low-dose integrations, not artifact/doublet)")
abline(h=1,col="grey45",lty=2); text(xs,ys,sprintf("%.2f",ys),pos=3,cex=0.7,col="firebrick")
points(max(xs)*30,mean(sigr),pch=17,col="steelblue4",cex=1.7); text(max(xs)*30,mean(sigr),sprintf("signal %.2f",mean(sigr)),pos=3,cex=0.8,col="steelblue4")
dev.off()
## per-guide panels (top 9 by below-gap n)
sel<-names(info)[order(-sapply(info,function(z)z$nbg))][1:9]
png(file.path(OUT,"dose_response_perguide.png"),width=1500,height=1300,res=140); par(mfrow=c(3,3),mar=c(4,4.2,3,1),mgp=c(2.3,0.7,0),cex.axis=0.9)
for(g in sel){z<-info[[g]];pb<-list(`2`=2,`3`=3,`4-5`=4:5,`6-9`=6:9,`10-19`=10:19,`20+`=20:9999);gx<-c();gy<-c()
  for(nm in names(pb)){s<-z$v%in%pb[[nm]];if(sum(s)<8)next;gx<-c(gx,mean(z$v[s]));gy<-c(gy,mean(z$r[s]))}
  plot(gx,gy,log="x",type="b",pch=16,col="firebrick",cex=1.3,ylim=c(0,1.15),xlim=c(1.7,z$ghi*1.5),xlab="guide UMI count",ylab="expr/baseline",main=sprintf("%s (gap %d-%d, n=%d)",gene(g),z$glo,z$ghi,z$nbg))
  abline(h=1,col="grey50",lty=2);points(z$ghi*1.1,z$sigr,pch=17,col="steelblue4",cex=1.5)}
dev.off(); cat("wrote dose_response.png, dose_response_perguide.png\n")

## (D) cross-dataset dose-response (survey triples; NT baseline or complement) --
dose_survey<-function(dir,has_nt){
  g<-as(readRDS(file.path(dir,"grna_matrix_aligned.rds")),"CsparseMatrix"); r<-as(readRDS(file.path(dir,"response_matrix.rds")),"CsparseMatrix")
  td<-read.csv(file.path(dir,"grna_target_data_frame.csv"),stringsAsFactors=FALSE); libx<-as.numeric(Matrix::colSums(r)); tmap<-setNames(td$grna_target,td$grna_id); rn<-rownames(r)
  nt<-td$grna_id[td$grna_target=="non-targeting"]; ntr<-match(nt,rownames(g)); ntr<-ntr[!is.na(ntr)]
  nm<-if(has_nt&&length(ntr)) Matrix::colSums(g[ntr,,drop=FALSE]>=30)>=1 else rep(TRUE,ncol(g))
  gvv<-g@i+1L; cvv<-rep.int(seq_len(ncol(g)),diff(g@p)); xvv<-g@x; ps<-split(seq_along(gvv),gvv)
  cand<-list(); for(gi in seq_len(nrow(g))){tg<-tmap[rownames(g)[gi]];if(is.na(tg)||tg=="non-targeting"||!(tg%in%rn))next;p<-ps[[as.character(gi)]];if(is.null(p))next;vv<-xvv[p];if(sum(vv>=30)>=5&&sum(vv>=2&vv<30)>=5)cand[[length(cand)+1]]<-list(gi=gi,tg=tg,p=p,vv=vv)}
  utg<-unique(vapply(cand,function(z)z$tg,"")); Est<-t(r[utg,,drop=FALSE]); colnames(Est)<-utg
  cc<-c();rr<-c();sg<-c(); for(z in cand){v<-integer(ncol(g));v[cvv[z$p]]<-z$vv;st<-which(v>=30);cp<-as.numeric(Est[,z$tg])/libx*1e4;base<-mean(cp[nm&v==0]);if(!is.finite(base)||base<0.5||mean(cp[st])/base>0.7)next;bg<-which(v>=1&v<30);cc<-c(cc,v[bg]);rr<-c(rr,cp[bg]/base);sg<-c(sg,cp[st]/base)}
  b<-list(`1`=1,`2`=2,`3`=3,`4-5`=4:5,`6-9`=6:9,`10-19`=10:19,`20-29`=20:29);ox<-c();oy<-c();for(nm2 in names(b)){s<-cc%in%b[[nm2]];if(sum(s)<20)next;ox<-c(ox,mean(cc[s]));oy<-c(oy,mean(rr[s]))}
  list(xs=ox,ys=oy,sig=mean(sg)) }
SV<-paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
en<-dose_survey(file.path(SV,"endoc_t2d_perturbseq_GSE273677/sceptre"),TRUE)
tc<-dose_survey(file.path(SV,"tcell_cd4_perturbseq_GSE314342/sceptre"),FALSE)   # slow (~7min)
rp<-list(xs=xs,ys=ys,sig=mean(sigr)); saveRDS(list(replogle=rp,endoc=en,tcell=tc),file.path(OUT,"dose_response_xdataset.rds"))
png(file.path(OUT,"dose_response_xdataset.png"),width=1080,height=700,res=140); par(mar=c(4.6,4.8,3.6,1),mgp=c(2.9,0.8,0),cex.lab=1.1)
plot(NA,log="x",xlim=c(1,700),ylim=c(0,1.12),xlab="guide UMI count (below-signal)",ylab="target expression / baseline",main="Dose-response across datasets with clean modes + expression")
abline(h=1,col="grey50",lty=2); al<-list(rp,en,tc); cols<-c("firebrick","seagreen4","purple3"); nms<-c("Replogle","EndoC","CD4 T-cell")
for(j in 1:3){z<-al[[j]];lines(z$xs,z$ys,type="b",pch=16,col=cols[j],lwd=1.6,cex=1.2);points(300,z$sig,pch=17,col=cols[j],cex=1.5)}
legend("bottomleft",bty="n",cex=0.9,pch=16,col=cols,legend=sprintf("%s (strong: %.2f)",nms,sapply(al,function(z)z$sig)))
dev.off(); cat("wrote dose_response_xdataset.png\n")
