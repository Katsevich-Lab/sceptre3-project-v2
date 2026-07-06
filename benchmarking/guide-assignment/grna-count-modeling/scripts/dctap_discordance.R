#!/usr/bin/env Rscript
# ============================================================================
# Zoom into the dctap datasets, where CLEANSER's <=2 per-cell ambient depth diverges most
# from fishash+'s rank-1 depth. Two views:
#   (1) per-cell scatter (CLEANSER vs fishash+ depth) coloured by the cell's AMBIENT mass sitting
#       at counts > 2 -- the UMIs CLEANSER's <=2 cap throws away. Explains the inverted-U.
#   (2) individual-guide count distributions for the guides whose ambient (unassigned) counts most
#       often exceed 2 -- with the depth-mixed Poisson ambient fit, the CLEANSER <=2 cutoff, and the
#       fishash+ ambient ceiling. The shaded 3..ceiling band = ambient UMIs CLEANSER discards.
# Uses cached fit + assignment (no refit).
# ============================================================================
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(tidyr) })
OUT <- "results/ambient_ceiling"; CACHE <- file.path(OUT, "fit_cache")
SUR <- paste0(.get_config_path("LOCAL_EXTERNAL_DATA_DIR"), "perturbseq-survey")
paths <- list(
  dctap_k562_highmoi = file.path(SUR,"dctap_k562_highmoi_GSE303901/grna_matrix.rds"),
  dctap_k562_lowmoi  = file.path(SUR,"dctap_k562_lowmoi_GSE303901/grna_matrix.rds"))

load_ds <- function(nm) {
  fit <- readRDS(file.path(CACHE, paste0(nm,"_fit.rds")))
  asg <- as(readRDS(file.path(CACHE, paste0(nm,"_assigned.rds"))), "CsparseMatrix")
  counts <- as(readRDS(paths[[nm]]), "CsparseMatrix"); storage.mode(counts@x) <- "double"
  list(nm=nm, fit=fit, asg=asg, counts=counts)
}

## ---- (1) per-cell scatter coloured by ambient-mass-above-2 -------------------------
cell_rows <- list()
D <- lapply(names(paths), load_ds)
for (d in D) {
  counts <- d$counts; asg <- d$asg
  fh <- as.numeric(d$fit$Cc)
  sub <- counts; sub@x[sub@x > 2] <- 0
  cl <- Matrix::colSums(sub)                                   # CLEANSER depth
  # ambient UMIs at counts>2 (unassigned entries with count>2) -- what the <=2 cap discards
  amb_gt2 <- counts * (counts > 2) * (1 - asg)
  disc <- Matrix::colSums(amb_gt2)
  lib  <- Matrix::colSums(counts)
  cell_rows[[d$nm]] <- data.frame(dataset=d$nm, fish=fh, clean=cl,
                                  disc_frac = disc / pmax(lib,1))
}
CELL <- do.call(rbind, cell_rows)
CELL <- CELL[is.finite(CELL$fish) & is.finite(CELL$clean) & CELL$fish>0 & CELL$clean>0, ]
p1 <- ggplot(CELL, aes(fish, clean, colour = pmin(disc_frac, 0.8))) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(alpha=0.35, size=0.5) +
  facet_wrap(~dataset, scales="free") +
  scale_x_log10() + scale_y_log10() +
  scale_colour_viridis_c(option="magma", direction=-1, name="cell's ambient\nUMIs at count>2\n(/ library)") +
  labs(title="dctap: CLEANSER (<=2) vs fishash+ per-cell ambient depth",
       subtitle="colour = fraction of the cell's UMIs in AMBIENT counts >2 (what CLEANSER's <=2 cap discards). Dark cells fall furthest below y=x.",
       x="fishash+ ambient depth  d_c (log)", y="CLEANSER ambient depth (sum counts <=2, log)") +
  theme_bw(base_size=12)
ggsave(file.path(OUT,"dctap_depth_scatter_zoom.png"), p1, width=12, height=5.5, dpi=150)

## ---- (2) most-discordant individual guides in dctap_highmoi ------------------------
dh <- D[[1]]                                                   # highmoi
counts <- dh$counts; asg <- dh$asg; fit <- dh$fit
G <- nrow(counts)
ct <- as(counts,"TsparseMatrix"); gi<-ct@i+1L; ci<-ct@j+1L; y<-as.numeric(ct@x)
a  <- as.logical(asg[cbind(gi,ci)]); a[is.na(a)]<-FALSE
amb_gt2 <- y>2 & !a                                           # ambient entries above the <=2 cap
umi_disc <- tapply(y[amb_gt2], factor(gi[amb_gt2], levels=seq_len(G)), sum)
umi_disc[is.na(umi_disc)] <- 0
rank_df <- data.frame(guide=rownames(counts), umi_discarded=as.numeric(umi_disc),
                      n_cells_amb_gt2 = as.numeric(tapply(amb_gt2, factor(gi,levels=seq_len(G)), sum)))
rank_df <- rank_df[order(-rank_df$umi_discarded), ]
write.csv(rank_df, file.path(OUT,"dctap_highmoi_guide_discordance.csv"), row.names=FALSE)
cat("Top discordant guides (most ambient UMIs above CLEANSER's <=2 cap):\n")
print(head(rank_df, 8), row.names=FALSE)

pgt <- read.csv(file.path(OUT,"per_guide_ambient_ceiling.csv"))
ceil_of <- function(gid) { v<-pgt$ambient_ceiling[pgt$dataset=="dctap_k562_highmoi" & pgt$guide==gid]; if(length(v)) v[1] else NA }

sel <- head(rank_df$guide, 6)
KMAX <- 60L; ks <- 0:KMAX
Tn <- fit$Tn
rows <- list(); barrows <- list()
for (gid in sel) {
  g <- match(gid, fit$guides); yv <- as.numeric(counts[g, ]); lam <- (fit$Rg[g]/Tn)*fit$Cc
  a <- as.logical(asg[g, ]); a[is.na(a)] <- FALSE          # signal (assigned)
  ambc <- !a                                               # ambient (unassigned) cells
  obs_amb <- vapply(ks, function(k) if(k<KMAX) sum(yv[ambc]==k) else sum(yv[ambc]>=k), numeric(1))
  obs_sig <- vapply(ks, function(k) if(k<KMAX) sum(yv[a]==k)    else sum(yv[a]>=k),    numeric(1))
  # depth-mixed Poisson summed over AMBIENT cells ONLY -> integrates to #ambient cells (fit the noise
  # to the noise cells; the signal cells are excluded, exactly as collab_eno1_gapfit.R does).
  dm  <- vapply(ks, function(k) if(k<KMAX) sum(dpois(k,lam[ambc])) else sum(1-ppois(k-1,lam[ambc])), numeric(1))
  cl <- ceil_of(gid)
  rows[[gid]]    <- data.frame(guide=gid, k=ks, depthmix=dm, ceiling=cl)
  barrows[[gid]] <- rbind(data.frame(guide=gid, k=ks, n=obs_amb, cls="ambient (unassigned)"),
                          data.frame(guide=gid, k=ks, n=obs_sig, cls="signal (assigned)"))
}
GD <- do.call(rbind, rows); BR <- do.call(rbind, barrows)
lab <- sprintf("%s  (ceiling %d)", sel, sapply(sel, ceil_of))
GD$facet <- factor(sprintf("%s  (ceiling %d)", GD$guide, GD$ceiling), levels=lab)
BR$facet <- factor(sprintf("%s  (ceiling %d)", BR$guide, sapply(as.character(BR$guide), ceil_of)), levels=lab)
BR$cls   <- factor(BR$cls, levels=c("signal (assigned)","ambient (unassigned)"))  # ambient drawn on top
band <- do.call(rbind, lapply(sel, function(gid){
  cl <- ceil_of(gid); data.frame(facet=sprintf("%s  (ceiling %d)", gid, cl), xmin=2.5, xmax=cl+0.5)}))
band$facet <- factor(band$facet, levels=lab)

p2 <- ggplot() +
  geom_rect(data=band, aes(xmin=xmin, xmax=xmax, ymin=0.5, ymax=Inf), fill="#f0c6b0", alpha=0.45) +
  geom_col(data=BR, aes(as.integer(k), n, fill=cls), width=0.9) +
  geom_line(data=GD, aes(as.integer(k), depthmix), colour="#D55E00", linewidth=0.75) +
  geom_vline(xintercept=2.5, colour="#0072B2", linewidth=0.6) +
  facet_wrap(~facet, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("ambient (unassigned)"="grey72","signal (assigned)"="#7570b3"), name=NULL) +
  scale_y_continuous(trans=scales::trans_new("log1p",log1p,expm1),
                     breaks=c(0,1,10,100,1000), labels=scales::label_comma()) +
  coord_cartesian(xlim=c(0, KMAX)) +
  labs(title="dctap_k562_highmoi: most CLEANSER-discordant guides (ambient fit over ambient cells)",
       subtitle="bars = observed cells split ambient (grey) vs signal (purple); orange = depth-mixed Poisson summed over AMBIENT cells only (integrates to #ambient cells). Blue = CLEANSER <=2 cutoff; shaded 3..ceiling = ambient CLEANSER discards.",
       x="gRNA UMI count", y="cells (log1p; floor 0)", fill=NULL) +
  theme_bw(base_size=12) + theme(legend.position="top")
ggsave(file.path(OUT,"dctap_highmoi_discordant_guides.png"), p2, width=13, height=7.4, dpi=140)

cat("\nwrote:\n  ", file.path(OUT,"dctap_depth_scatter_zoom.png"),
    "\n  ", file.path(OUT,"dctap_highmoi_discordant_guides.png"),
    "\n  ", file.path(OUT,"dctap_highmoi_guide_discordance.csv"), "\n")
