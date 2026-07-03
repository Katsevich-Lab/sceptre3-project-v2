## ============================================================================
## Is the ambient null overdispersed (needs NB), or is Poisson correct?
## The Replogle low-count rate is ~2.2x the rank-1 Poisson. This script works the
## whole arc:
##  (1) count-2 obs/exp under the rank-1 Poisson (~2.18x).
##  (2) rank-2 mean test: a richer mean does NOT absorb it (2.18 -> 2.15) => not
##      hidden guide x cell structure.
##  (3) ML negative-binomial dispersion fit on the low-count noise + GOF (rho~1.66;
##      one param fits counts 0-3 at obs/NB~1).
##  (4) depth_fix with the NB null through the sandbox: removes ~half the count-2
##      calls, barely touches count-3+ (which are doublet-context / real).
##  (5) THE DECISIVE TEST -- barnyard provable ambient (wrong-species guides) is
##      Poisson at count-2 (obs/exp ~0.90 after de-doubleting), NOT 2.18. So the
##      Replogle excess is WEAK-INTEGRATION SIGNAL, not overdispersed soup, and the
##      NB would suppress real weak integrations. => KEEP POISSON.
## Run from nonparametric-thresholds/.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats)})
OUT<-"results/global_ambient_poisson"
D<-path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
mc<-as(readRDS(D),"CsparseMatrix"); G<-nrow(mc); C<-ncol(mc)
gv<-mc@i+1L; cv<-rep.int(seq_len(C),diff(mc@p)); xv<-mc@x
rowN<-as.numeric(rowSums(mc)); colN<-as.numeric(colSums(mc))
fs<-function(val,idx,n){o<-numeric(n);s<-rowsum(val,idx);o[as.integer(rownames(s))]<-s[,1];o}
d<-colN;d[d==0]<-1e-6;a<-rowN/sum(rowN);mask<-logical(length(xv))
for(o in 1:6){for(it in 1:15){
 if(any(mask)){mrN<-fs(xv[mask],gv[mask],G);mrD<-fs(d[cv[mask]],gv[mask],G)}else{mrN<-numeric(G);mrD<-numeric(G)}
 a<-(rowN-mrN)/(sum(d)-mrD);a[!is.finite(a)|a<0]<-0;a<-a/sum(a)
 if(any(mask)){mcN<-fs(xv[mask],cv[mask],C);mcA<-fs(a[gv[mask]],cv[mask],C)}else{mcN<-numeric(C);mcA<-numeric(C)}
 d<-(colN-mcN)/(1-mcA);d[!is.finite(d)|d<0]<-0}
 mask<-ppois(xv-1,a[gv]*d[cv],lower.tail=FALSE)<1e-6}

## (1)+(2) count-2 obs/exp: rank-1 vs rank-2 (unsupervised, no batch) --------
is2<-xv==2; n2<-tabulate(cv[is2],nbins=C); set.seed(1); samp<-sample(seq_len(C),30000)
exp1<-sapply(samp,function(c){sum(dpois(2,a*d[c]))}); obs<-sum(n2[samp])
cat(sprintf("(1) rank-1 count-2 obs/exp = %.2f\n", obs/sum(exp1)))
u<-!mask; gn<-gv[u];cn<-cv[u];xn<-xv[u]; gm<-gv[mask];cm<-cv[mask]
W<-cbind(a,a*0.1+mean(a)*1e-3); H<-rbind(d,d*0.1+mean(d)*1e-3)
for(it in 1:40){ Mn<-W[gn,1]*H[1,cn]+W[gn,2]*H[2,cn]; rr<-xn/pmax(Mn,1e-12)
  for(k in 1:2){H[k,]<-H[k,]*fs(W[gn,k]*rr,cn,C)/pmax(sum(W[,k])-fs(W[gm,k],cm,C),1e-12);H[k,][!is.finite(H[k,])|H[k,]<0]<-0}
  Mn<-W[gn,1]*H[1,cn]+W[gn,2]*H[2,cn]; rr<-xn/pmax(Mn,1e-12)
  for(k in 1:2){W[,k]<-W[,k]*fs(H[k,cn]*rr,gn,G)/pmax(sum(H[k,])-fs(H[k,cm],gm,G),1e-12);W[,k][!is.finite(W[,k])|W[,k]<0]<-0} }
exp2<-sapply(samp,function(c){mu<-W[,1]*H[1,c]+W[,2]*H[2,c];sum(dpois(2,mu))})
cat(sprintf("(2) rank-2 count-2 obs/exp = %.2f  (a richer mean does NOT absorb it => not mean structure)\n", obs/sum(exp2)))

## (3) ML negative-binomial dispersion on the low-count noise + GOF -----------
la<-log(a[a>0]);ld<-log(d[d>0]);ha<-hist(la,breaks=140,plot=FALSE);hd<-hist(ld,breaks=320,plot=FALSE)
mu_b<-exp(as.vector(outer(ha$mids,hd$mids,"+")));w_b<-as.vector(outer(ha$counts,hd$counts,"*"));kp<-w_b>0;mu_b<-mu_b[kp];w_b<-w_b[kp];Np<-sum(w_b)
K<-6; vnm<-xv[!mask]; nobs<-numeric(K+1); for(k in 1:K) nobs[k+1]<-sum(vnm==k); nobs[1]<-Np-length(xv)
Nnm<-Np-sum(mask); ngt<-Nnm-sum(nobs)
pk<-function(rho,k) if(rho<=0) sum(w_b*dpois(k,mu_b))/Np else sum(w_b*dnbinom(k,mu=mu_b,size=1/rho))/Np
nLL<-function(rho){p<-sapply(0:K,function(k)pk(rho,k));pg<-max(1-sum(p),1e-12);-(sum(nobs*log(pmax(p,1e-300)))+ngt*log(pg))}
rho_ml<-optimize(nLL,c(1e-4,4))$minimum
cat(sprintf("\n(3) ML negative-binomial dispersion rho = %.3f\n    GOF (obs/pred): count  Poisson   NB\n", rho_ml))
for(k in 0:4) cat(sprintf("      %d   obs/Pois=%.2f  obs/NB=%.2f\n",k,nobs[k+1]/(pk(0,k)*Np),nobs[k+1]/(pk(rho_ml,k)*Np)))
saveRDS(list(a=a,d=d,rho_ml=rho_ml),file.path(OUT,"ambient_ml_fit.rds"))

## (4) NB sandbox: run depth_fix hyper vs nb(rho_ml), lower-mode + residual ----
source("scripts/contingency_method.R")
Ap<-as(contingency_assign(mc,q=0.05,refit=10,min_count=2,cell_margin="ambient",tail="hyper",fdr="GS")$assigned,"CsparseMatrix")
An<-as(contingency_assign(mc,q=0.05,refit=10,min_count=2,cell_margin="ambient",tail="nb",fdr="GS",rho_fixed=rho_ml)$assigned,"CsparseMatrix")
pg<-read.csv(file.path(OUT,"perguide_replogle_rd7.csv"),stringsAsFactors=FALSE);pg<-pg[pg$gap_lo>=2,];ridx<-match(pg$guide,rownames(mc))
nstrong<-as.numeric(colSums(mc>=30))
c2p<-c2n<-c3p<-c3n<-sigk<-sigt<-lown<-dbn<-0
for(i in seq_len(nrow(pg))){r<-ridx[i];v<-as.numeric(mc[r,]);lo<-pg$gap_lo[i];ap<-as.logical(Ap[r,]);an<-as.logical(An[r,])
  c2p<-c2p+sum(ap&v==2);c2n<-c2n+sum(an&v==2);c3p<-c3p+sum(ap&v>=3&v<=lo);c3n<-c3n+sum(an&v>=3&v<=lo)
  sigk<-sigk+sum(an&v>=pg$gap_hi[i]);sigt<-sigt+sum(v>=pg$gap_hi[i])
  ln<-which(an&v>=2&v<=lo);lown<-lown+length(ln);dbn<-dbn+sum(nstrong[ln]-(v[ln]>=30)>=1)}
cat(sprintf("\n(4) NB sandbox: count-2 calls %d->%d (%.0f%% removed); count-3+ calls %d->%d (%.0f%% removed); signal recall %.0f%%\n",
    c2p,c2n,100*(1-c2n/c2p),c3p,c3n,100*(1-c3n/c3p),100*sigk/sigt))
cat(sprintf("    of NB-residual lower-mode calls, %.0f%% are doublet-context (another strong guide)\n",100*dbn/lown))

## (5) DECISIVE: barnyard provable ambient (wrong-species) count-2 obs/exp -----
R<-"external/repro_work"
load_sample<-function(sample,purity=0.9){ gm<-as(readMM(file.path(R,paste0(sample,"_grna_counts.mtx"))),"CsparseMatrix")
  meta<-read.csv(file.path(R,paste0(sample,"_meta.csv")));guides<-read.csv(file.path(R,paste0(sample,"_guides.csv")))
  sg<-meta$homo_sum_gex+meta$mus_sum_gex;fh<-meta$homo_sum_gex/sg
  qc<-(meta$mito_sum/sg<.15)&(meta$features_gex<=6000)&(sg<=20000)&(meta$features_gex>=1500)&(sg>=3500)&(fh<(1-purity)|fh>purity);qc[is.na(qc)]<-FALSE
  list(counts=gm[,qc,drop=FALSE],guide_homo=guides$guide_type=="homo_guide",cell_homo=fh[qc]>purity) }
dc<-load_sample("mix0hr_DirectCapture",0.90); gsel<-!dc$guide_homo
cat("\n(5) DECISIVE -- barnyard provable ambient (mouse guides in human cells):\n")
for(cut in c(Inf,0.10,0.02)){ cn<-dc$counts; wf<-Matrix::colSums(cn[gsel,])/pmax(Matrix::colSums(cn),1)
  sub<-cn[gsel, dc$cell_homo & wf<=cut, drop=FALSE]; rs<-Matrix::rowSums(sub);cs<-Matrix::colSums(sub);tot<-sum(rs)
  obs2<-sum(sub@x==2); exp2<-sum(outer(rs,cs,function(x,y){dpois(2,x*y/tot)}))
  cat(sprintf("    de-doublet cut=%-4s : count-2 obs/exp = %.2f (maxrate=%.2f)\n",ifelse(is.finite(cut),cut,"none"),obs2/exp2,max(rs)*max(cs)/tot)) }
cat("    => de-doubleted pure soup ~0.90 (POISSON), vs Replogle 2.18 => Replogle excess is weak-integration signal. KEEP POISSON.\n")
