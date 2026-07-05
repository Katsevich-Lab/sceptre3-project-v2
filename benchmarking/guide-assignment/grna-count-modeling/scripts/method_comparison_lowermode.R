## ============================================================================
## Do competing methods share the lower-mode dilution? (fishash vs depth_fix)
## Runs the REAL fishash package and depth_fix on Replogle, and for the 5 example
## well-separated guides compares the on-target effect estimate from the gap-based
## threshold vs all of each method's calls. Both reach the count-2 floor; fishash
## (raw-library margin) dilutes less than depth_fix (denoised-depth margin), but
## both share the problem. -> results/global_ambient_poisson/method_compare_effectsize.png
## Run from nonparametric-thresholds/.
## ============================================================================
suppressMessages({library(Matrix); library(sparseMatrixStats); library(fishash); library(SummarizedExperiment); library(ondisc); library(sceptre)})
source("scripts/contingency_method.R")
OUT<-"results/global_ambient_poisson"
Dm<-path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")
Dp<-path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline")
mc<-as(readRDS(Dm),"CsparseMatrix")

cat("running REAL fishash...\n"); res<-fishash(mc,refit=10,padj_cutoff=0.05,exclude_empty=TRUE)
Afh<-as(assay(res,"assigned"),"CsparseMatrix"); if(is.null(rownames(Afh))) rownames(Afh)<-rownames(mc)
cat("running depth_fix (Poisson/hyper)...\n"); Adf<-as(contingency_assign(mc,q=0.05,refit=10,min_count=2,cell_margin="ambient",tail="hyper",fdr="GS")$assigned,"CsparseMatrix")
cat(sprintf("total assignments: fishash=%d (%.2f/cell)  depth_fix=%d (%.2f/cell)\n",sum(Afh),sum(Afh)/ncol(mc),sum(Adf),sum(Adf)/ncol(mc)))

pg<-read.csv(file.path(OUT,"perguide_replogle_rd7.csv"),stringsAsFactors=FALSE); pg<-pg[pg$gap_lo>=2,]; ridx<-match(pg$guide,rownames(mc))
lowfrac<-function(A){ nl<-nt<-0; for(i in seq_len(nrow(pg))){ r<-ridx[i]; v<-as.numeric(mc[r,]); a<-as.logical(A[r,]); nl<-nl+sum(a&v>=2&v<=pg$gap_lo[i]); nt<-nt+sum(a) }; nl/nt }
cat(sprintf("pooled lower-mode fraction on clean-gap guides: fishash=%.1f%%  depth_fix=%.1f%%\n\n",100*lowfrac(Afh),100*lowfrac(Adf)))

resp<-initialize_odm_from_backing_file(file.path(Dp,"response.odm")); so<-readRDS(file.path(Dp,"sceptre_object.rds"))
lib<-exp(so@covariate_matrix[,"log(response_n_umis)"]); tdf<-so@grna_target_data_frame
ntg<-tdf$grna_id[tdf$grna_target=="non-targeting"]; ntmask<-rep(FALSE,ncol(mc)); for(g in ntg) ntmask<-ntmask|(as.numeric(mc[g,])>=30)
picks<-c("1290_CCDC6_P1P2_ENSG00000108091","7424_RPAP3_P1P2_ENSG00000005175","968_C12orf60_P1P2_ENSG00000182993","9676_USP8_P1P2_ENSG00000138592","9175_TOP2A_P1P2_ENSG00000131747")
gene<-function(g) sub("^[0-9]+_([^_]+)_.*","\\1",g)
d<-data.frame()
for(g in picks){ r<-match(g,rownames(mc)); tgt<-sub(".*_(ENSG[0-9]+)$","\\1",g)
  v<-as.numeric(mc[r,]); cp<-as.numeric(resp[tgt,])/lib*1e4; base<-mean(cp[ntmask&v==0]); gh<-pg$gap_hi[pg$guide==g]
  d<-rbind(d,data.frame(gene=gene(g),eff_gap=mean(cp[v>=gh])/base,eff_fh=mean(cp[as.logical(Afh[r,])])/base,eff_df=mean(cp[as.logical(Adf[r,])])/base)) }
print(d,row.names=FALSE)
png(file.path(OUT,"method_compare_effectsize.png"),width=1200,height=680,res=140)
par(mar=c(4.6,4.8,3.6,1),mgp=c(2.9,0.8,0),cex.axis=0.95,cex.lab=1.05); nD<-nrow(d); x<-seq_len(nD); w<-0.26
plot(NA,xlim=c(.5,nD+.5),ylim=c(0,1.15),xaxt="n",xlab="",ylab="target expression / NT baseline",
     main="Do competing methods dilute too? gap-based vs all calls\n(1.0 = no knockdown; both reach the count-2 floor)")
abline(h=1,col="grey45",lty=2)
rect(x-1.5*w,0,x-0.5*w,d$eff_gap,col="grey55",border=NA); rect(x-0.5*w,0,x+0.5*w,d$eff_fh,col="darkorange2",border=NA); rect(x+0.5*w,0,x+1.5*w,d$eff_df,col="firebrick",border=NA)
for(k in 1:nD){ text(x[k]-w,d$eff_gap[k],sprintf("%.2f",d$eff_gap[k]),pos=3,cex=0.66,col="grey35"); text(x[k],d$eff_fh[k],sprintf("%.2f",d$eff_fh[k]),pos=3,cex=0.66,col="darkorange3"); text(x[k]+w,d$eff_df[k],sprintf("%.2f",d$eff_df[k]),pos=3,cex=0.66,col="firebrick") }
axis(1,at=x,labels=d$gene,padj=0.5,cex.axis=0.9)
legend("topleft",bty="n",cex=0.85,fill=c("grey55","darkorange2","firebrick"),legend=c("gap-based (signal mode only)","all fishash calls","all depth_fix calls"))
dev.off(); cat("wrote method_compare_effectsize.png\n")
