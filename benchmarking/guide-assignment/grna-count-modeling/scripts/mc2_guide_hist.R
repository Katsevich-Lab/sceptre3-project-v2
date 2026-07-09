#!/usr/bin/env Rscript
suppressMessages({library(Matrix); library(ggplot2); library(dplyr)})
source("scripts/datasets.R")   # load_grna_matrix, load_cleanser_assignment

DSETS <- c("dctap_k562_highmoi","endoc_t2d","mccutcheon","gastric_organoid")  # depth desc, non-barnyard, clean chemistry
DISP  <- c(mccutcheon = "mccutcheon (CD8 T-cell CRISPRi)")      # more-informative display names
disp  <- function(ds) if (ds %in% names(DISP)) DISP[[ds]] else ds
ADIR  <- "results/ambient_ceiling/writeup/assign"

# per-guide triplet helper: returns list with gi/cj/xv and ambient masks per method
load_ds <- function(ds) {
  gm <- load_grna_matrix(ds)
  Af <- as(readRDS(file.path(ADIR, paste0(ds,"_fishash.rds"))),     "lgCMatrix")
  Ap <- as(readRDS(file.path(ADIR, paste0(ds,"_fishashplus.rds"))), "lgCMatrix")
  ct <- as(gm, "TsparseMatrix"); gi <- ct@i+1L; cj <- ct@j+1L; xv <- ct@x
  idx <- cbind(gi,cj)
  amb <- list()                                                 # CLEANSER optional (missing for big datasets)
  cf <- c(file.path("results/ambient_ceiling/cleanser_calls_cs", paste0(ds,".csv")),  # prefer chemistry-correct --cs
          file.path("results/ambient_ceiling/cleanser_calls",    paste0(ds,".csv")))
  cf <- cf[file.exists(cf)][1]
  if (!is.na(cf)) { cl <- read.csv(cf); gg <- match(cl$guide,rownames(gm)); cc <- match(cl$cell,colnames(gm)); ok <- !is.na(gg)&!is.na(cc)
    Ac <- sparseMatrix(i=gg[ok], j=cc[ok], x=TRUE, dims=dim(gm), dimnames=dimnames(gm))
    amb$cleanser <- !as(Ac,"lgCMatrix")[idx] }
  amb$fishash <- !Af[idx]; amb$`fishash+` <- !Ap[idx]           # ambient among nonzero
  list(gm=gm, G=nrow(gm), guides=rownames(gm), gi=gi, cj=cj, xv=xv, amb=amb)
}

pick_guides <- function(d) {
  G <- d$G
  nz  <- tabulate(d$gi, nbins=G)
  ceil_p <- rep(0, G); a <- d$amb$`fishash+`
  m <- tapply(d$xv[a], factor(d$gi[a], levels=seq_len(G)), max); ceil_p[!is.na(m)] <- m[!is.na(m)]
  # eligibility: enough nonzero cells for a real histogram (relax if too few guides qualify)
  for (thr in c(200,100,50,20,1)) { elig <- which(nz>=thr); if (length(elig)>=3) break }
  ord <- elig[order(ceil_p[elig], nz[elig])]                    # by ceiling asc, tie-break nz
  k <- length(ord); sel <- ord[unique(pmax(1, round(c(0.10,0.5,0.92)*k)))]
  sel <- sel[!duplicated(sel)]; while(length(sel)<3) sel <- c(sel, setdiff(ord, sel)[1])
  data.frame(guide_idx=sel, guide=d$guides[sel], nz=nz[sel], ceil=ceil_p[sel],
             rank=c("low","median","high")[seq_along(sel)])
}

shorten <- function(g) { g <- sub("^CROPseq_dCas9_DS_","",g); g <- sub("-P1P2$","",g)
  ifelse(nchar(g)>20, paste0(substr(g,1,9),"..",substr(g,nchar(g)-7,nchar(g))), g) }

bars <- list(); lines <- list(); picks <- list()
for (ds in DSETS) {
  d <- load_ds(ds); pg <- pick_guides(d); pg$dataset <- ds; picks[[ds]] <- pg
  for (r in seq_len(nrow(pg))) {
    g <- pg$guide_idx[r]; sel <- d$gi==g; vals <- d$xv[sel]
    amb <- lapply(d$amb, function(a) a[sel])
    maxv <- max(vals); nb <- max(6, min(22, ceiling(12*log10(maxv+1))))
    edges <- seq(log10(0.9), log10(maxv*1.12), length.out=nb+1)
    bin <- cut(log10(vals), edges, labels=FALSE, include.lowest=TRUE)
    tot <- tabulate(bin, nbins=nb)
    lab <- sprintf("%s\n%s  (ceil=%d, n=%d)", disp(ds), shorten(pg$guide[r]), pg$ceil[r], pg$nz[r])
    ord_key <- r*10 + match(ds,DSETS)                 # rank-major -> rows = spread level, cols = dataset
    bars[[length(bars)+1]] <- data.frame(panel=lab, ord=ord_key,
      xmin=10^edges[1:nb], xmax=10^edges[2:(nb+1)], total=tot)
    ctr <- 10^((edges[1:nb]+edges[2:(nb+1)])/2)
    for (mth in names(amb)) {
      na <- tabulate(bin[amb[[mth]]], nbins=nb)
      keep <- tot>0
      lines[[length(lines)+1]] <- data.frame(panel=lab, ord=ord_key, x=ctr[keep],
        method=mth, n=na[keep])
    }
  }
}
picks_df <- bind_rows(picks); print(picks_df[,c("dataset","rank","guide","nz","ceil")], row.names=FALSE)
bars_df <- bind_rows(bars); lines_df <- bind_rows(lines)
lev <- bars_df %>% distinct(panel, ord) %>% arrange(ord) %>% pull(panel)
bars_df$panel  <- factor(bars_df$panel, levels=lev)
lines_df$panel <- factor(lines_df$panel, levels=lev)
lines_df$method <- factor(lines_df$method, levels=c("cleanser","fishash","fishash+"))
pal   <- c("cleanser"="#e34948","fishash"="#eda100","fishash+"="#2a78d6")
SHP   <- c("cleanser"=17,"fishash"=15,"fishash+"=16)             # triangle / square / circle
DODGE <- c("cleanser"=1/1.07,"fishash"=1,"fishash+"=1.07)        # small log-x offset so coincident lines separate
FLOOR <- 0.5                                  # log-y floor; ambient count of 0 drawn here
bars_df  <- subset(bars_df, total>0)
lines_df$n_plot <- pmax(lines_df$n, FLOOR)
lines_df$xd <- lines_df$x * DODGE[as.character(lines_df$method)]

p <- ggplot() +
  geom_rect(data=bars_df, aes(xmin=xmin,xmax=xmax,ymin=FLOOR,ymax=total),
            fill="grey88", color="grey74", linewidth=0.12) +
  geom_line(data=lines_df, aes(xd,n_plot,color=method), linewidth=0.5, alpha=0.9) +
  geom_point(data=lines_df, aes(xd,n_plot,color=method,shape=method), size=1.15, alpha=0.95) +
  facet_wrap(~panel, scales="free", ncol=4) +
  scale_x_log10() + scale_y_log10() + scale_color_manual(values=pal) + scale_shape_manual(values=SHP) +
  labs(x="guide UMI counts across cells",
       y="number of cells", color="ambient calls", shape="ambient calls",
       title="Per-guide count histograms with ambient calls by method (grey = all cells)") +
  theme_bw(base_size=13) +
  theme(legend.position="bottom", strip.text=element_text(size=10.5),
        plot.title=element_text(size=15),
        axis.title=element_text(size=14),
        legend.text=element_text(size=13), legend.title=element_text(size=13))
ggsave("results/ambient_ceiling/guide_hist_methods.png", p, width=14, height=9, dpi=150)
cat("\nwrote results/ambient_ceiling/guide_hist_methods.png\n")
