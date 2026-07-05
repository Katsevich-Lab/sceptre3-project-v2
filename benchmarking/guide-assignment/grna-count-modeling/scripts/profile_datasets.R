## Profile all available Perturb-seq datasets: dimensions, MOI, library structure,
## gene-expression availability, negative-control guides. Anti-overfitting landscape.
suppressMessages(library(Matrix))
base <- path.expand("~/data/external/perturbseq-survey")
repl <- path.expand("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre")
OUT  <- "results/method_decision"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

find_grna <- function(dd){
  for (p in c(file.path(dd,"sceptre","grna_matrix_aligned.rds"),
              file.path(dd,"grna_matrix.rds"),
              file.path(dd,"grna_matrix.rds")))
    if (file.exists(p)) return(p)
  NA
}
find_tdf <- function(dd){
  p <- file.path(dd,"sceptre","grna_target_data_frame.csv"); if(file.exists(p)) return(p)
  p <- file.path(dd,"grna_target_data_frame.csv"); if(file.exists(p)) return(p); NA
}
find_resp <- function(dd){
  p <- file.path(dd,"sceptre","response_matrix.rds"); if(file.exists(p)) return(p); NA
}

profile <- function(name, dd){
  gp <- find_grna(dd); if(is.na(gp)) return(NULL)
  g <- readRDS(gp); if(!inherits(g,"CsparseMatrix")) g <- as(g,"CsparseMatrix")
  # orient guides x cells: assume more cells than guides
  if(nrow(g) > ncol(g)) g <- t(g)
  G <- nrow(g); C <- ncol(g)
  L <- colSums(g)
  mxv <- sparseMatrixStats::colMaxs(g)
  topfrac <- ifelse(L>0, mxv/L, NA)                      # ~1 => low MOI
  # MOI proxy (sparse-native): # guides with count >= max(3, 0.1*cellmax) per cell
  gv <- g@i + 1L; cv <- rep.int(seq_len(C), diff(g@p)); xv <- g@x
  thr <- pmax(3, 0.1*mxv)
  above <- xv >= thr[cv]
  ng_above <- tabulate(cv[above], nbins=C)
  tp <- find_tdf(dd); nNT <- NA; ntargets <- NA
  if(!is.na(tp)){ tdf <- read.csv(tp)
    tcol <- if("grna_target" %in% colnames(tdf)) "grna_target" else colnames(tdf)[2]
    nNT <- sum(grepl("non.?target|^NT|safe|control", tdf[[tcol]], ignore.case=TRUE))
    ntargets <- length(unique(tdf[[tcol]]))
  }
  rp <- find_resp(dd); hasGE <- !is.na(rp); geUMI <- NA; nGenes <- NA
  if(hasGE){ r <- readRDS(rp); nGenes <- nrow(r); geUMI <- median(colSums(r)) }
  data.frame(dataset=name, guides=G, cells=C,
             med_lib=median(L), med_cellmax=median(mxv),
             med_topfrac=round(median(topfrac,na.rm=TRUE),3),
             MOI_med=median(ng_above), MOI_mean=round(mean(ng_above),2),
             n_NT_guides=nNT, n_targets=ntargets,
             hasGE=hasGE, nGenes=nGenes, med_geUMI=geUMI)
}

rows <- list()
dsdirs <- list.dirs(base, recursive=FALSE)
for (dd in dsdirs){
  nm <- basename(dd)
  cat("profiling", nm, "...\n")
  r <- tryCatch(profile(nm, dd), error=function(e){cat("  ERR:", conditionMessage(e), "\n"); NULL})
  if(!is.null(r)) rows[[nm]] <- r
}
# replogle
cat("profiling replogle-rd7 ...\n")
r <- tryCatch(profile("replogle_rd7", dirname(repl)), error=function(e){cat("ERR",conditionMessage(e),"\n");NULL})
# replogle grna is at sceptre/grna_matrix.rds
if(is.null(r)){
  g <- as(readRDS(file.path(repl,"grna_matrix.rds")),"CsparseMatrix"); if(nrow(g)>ncol(g)) g<-t(g)
  L<-colSums(g); mxv<-sparseMatrixStats::colMaxs(g)
  tdf<-read.csv(file.path(repl,"grna_target_data_frame.csv"))
  rows[["replogle_rd7"]] <- data.frame(dataset="replogle_rd7",guides=nrow(g),cells=ncol(g),
    med_lib=median(L),med_cellmax=median(mxv),med_topfrac=round(median(mxv/pmax(L,1)),3),
    MOI_med=NA,MOI_mean=NA,n_NT_guides=sum(grepl("non.?target",tdf$grna_target,ignore.case=TRUE)),
    n_targets=length(unique(tdf$grna_target)),hasGE=FALSE,nGenes=NA,med_geUMI=NA)
} else rows[["replogle_rd7"]] <- r

tab <- do.call(rbind, rows); rownames(tab) <- NULL
write.csv(tab, file.path(OUT,"dataset_profiles.csv"), row.names=FALSE)
print(tab, row.names=FALSE)
