## ENO1 versions of the two RAMAC figures (gap histogram + ambient-mode depth-mixed-Poisson fit), so we
## can compare RAMAC vs ENO1 as the §3 worked example. Uses the cached denoised rank-1 field (a_g d_c).
suppressMessages({library(Matrix); library(ggplot2)})
OUT<-"results/collaborator_writeup"
mc<-as(readRDS(path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre/grna_matrix.rds")),"CsparseMatrix")
ad<-readRDS("results/replogle_ambient_poisson/replogle_rank1_ad.rds")
G<-"2640_ENO1_P1P2_ENSG00000074800"; idx<-match(G,rownames(mc)); LO<-17L; HI<-63L
v<-as.numeric(mc[idx,]); lam<-ad$a[idx]*ad$d          # per-cell denoised ambient rate a_ENO1 * d_c

## ---- (1) gap histogram: full count distribution, low mode / empty gap / high mode ---------------
tab<-table(v[v>0]); cnt<-as.integer(names(tab)); nc<-as.integer(tab)
mode<-ifelse(cnt<=LO,"low mode (≤ gap_lo = 17)", ifelse(cnt>=HI,"high mode (≥ gap_hi = 63)","gap"))
dh<-data.frame(cnt=cnt,nc=nc,mode=mode)
p1<-ggplot(dh, aes(cnt, nc)) +
  annotate("rect", xmin=LO+0.5, xmax=HI-0.5, ymin=0.6, ymax=Inf, fill="grey88") +
  annotate("text", x=sqrt((LO+0.5)*(HI-0.5)), y=6, label="empty gap\n(no cells\ncount 18–62)", colour="grey35", size=4) +
  geom_segment(aes(xend=cnt, yend=0.6, colour=mode), linewidth=0.5) + geom_point(aes(colour=mode), size=1.6) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c("low mode (≤ gap_lo = 17)"="#0072B2","high mode (≥ gap_hi = 63)"="#D55E00","gap"="grey60")) +
  labs(x="ENO1 gRNA UMI count (log)", y="number of cells (log)", colour=NULL,
       title="ENO1 count distribution") +
  theme_bw(base_size=15) + theme(legend.position="top")
ggsave(file.path(OUT,"eno1_gap_hist.png"), p1, width=9, height=5.2, dpi=130)

## ---- (2) ambient-mode fit: below-gap observed vs depth-mixed Poisson (a_g d_c) ------------------
amb<-which(v<=LO); la<-lam[amb]; ks<-0:LO
obs <-vapply(ks, function(k) sum(v[amb]==k), numeric(1))
pred<-vapply(ks, function(k) sum(dpois(k, la)), numeric(1))
reg<-ks>=2; obs_r<-sum(obs[reg]); pred_r<-sum(pred[reg])
cat(sprintf("ENO1 low-mode counts 2–%d: obs=%d vs pred=%.0f (%.1fx); count-2 obs=%d pred=%.0f (%.1fx)\n",
    LO, obs_r, pred_r, obs_r/pred_r, obs[3], pred[3], obs[3]/pred[3]))
df<-data.frame(k=ks, observed=obs, model=pred)
p2<-ggplot(df, aes(factor(k), observed)) +
  geom_col(fill="grey80", width=0.7) +
  geom_line(aes(y=model, group=1), colour="#D55E00", linewidth=0.9) +
  geom_point(aes(y=model), colour="#D55E00", size=2.4) +
  annotate("text", x=length(ks), y=max(obs)*0.6, hjust=1, size=4.5, colour="grey20",
           label=sprintf("counts 2–17: %d obs\nvs ~%.0f predicted (%.1f×)", obs_r, pred_r, obs_r/pred_r)) +
  scale_y_continuous(trans=scales::trans_new("log1p", log1p, expm1),
                     breaks=c(0,1,10,100,1000,10000,1e5,1e6), labels=scales::label_comma(),
                     expand=expansion(mult=c(0,0.06))) +
  labs(x="ENO1 gRNA UMI count", y="number of cells (log1p — floor at 0)",
       title="ENO1 low mode versus fit") +
  theme_bw(base_size=15)
ggsave(file.path(OUT,"eno1_ambient_fit.png"), p2, width=10, height=5.0, dpi=130)
cat("wrote eno1_gap_hist.png + eno1_ambient_fit.png\n")
