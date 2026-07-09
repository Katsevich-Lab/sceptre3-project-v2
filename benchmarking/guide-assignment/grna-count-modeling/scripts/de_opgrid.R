#!/usr/bin/env Rscript
## ============================================================================
## de_opgrid.R  --  DE pushthrough over ASSIGNMENT OPERATING POINTS (not methods).
## Settles the error-control questions empirically: for the SAME fishash+ p-values,
## vary the multiple-testing correction (GS / BH / per-guide-BH / BY), the FDR
## level q, and the per-cell QC (singleton filter vs keep-all), then run sceptre's
## REAL run_calibration_check (neg-control Type-I) + run_power_check (on-target power).
##
## The fishash+ GS-converged fit is computed ONCE per dataset (cached); each
## operating point is a cheap re-threshold of the shared p-values -> A -> sceptre DE.
## Downstream sceptre block is faithful to scripts/de_native.R.
##
## Usage: Rscript scripts/de_opgrid.R <dataset> <max_cells> <out_dir>
## Outputs: <out_dir>/summaries.csv (one row per operating point; checkpointed/resumable),
##          <out_dir>/fit.rds (cached fishash+ fit), plus per-op calibration/power csvs.
## ============================================================================
suppressPackageStartupMessages({
  library(Matrix); library(sceptre); library(sparseMatrixStats); library(matrixStats)
})
source(file.path(getwd(),"scripts","sim_de_lib.R"))
source(file.path(getwd(),"scripts","contingency_method.R"))
source(file.path(getwd(),"scripts","error_control_experiments.R"))   # cassign / apply_corr / apply_perguide
`%||%` <- function(a,b) if(is.null(a)) b else a

args <- commandArgs(trailingOnly=TRUE)
dataset   <- args[1]
max_cells <- as.integer(if(length(args)>=2) args[2] else "60000")
out_dir   <- args[3]; dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
ncpus     <- as.integer(Sys.getenv("NCPUS","6"))
NCAL      <- as.integer(Sys.getenv("NCAL","2000"))   # calibration pairs

## -------- operating-point grid (filtered by MOI below) --------
GRID <- rbind(
  data.frame(corr="GS",         q=0.05, qc="none"),      # fishash default
  data.frame(corr="BH",         q=0.05, qc="none"),
  data.frame(corr="perguideBH", q=0.05, qc="none"),
  data.frame(corr="BY",         q=0.05, qc="none"),
  data.frame(corr="BH",         q=0.01, qc="none"),      # q-sweep (cost-asymmetry test)
  data.frame(corr="BH",         q=0.10, qc="none"),
  data.frame(corr="BH",         q=0.20, qc="none"),
  data.frame(corr="BH",         q=0.40, qc="none"),
  data.frame(corr="perguideBH", q=0.20, qc="none"),
  data.frame(corr="BH",         q=0.10, qc="singleton")  # hard per-cell QC (low-MOI only)
)

## ============================ load dataset (once) ============================
spec <- REAL_DATASETS[[dataset]]
real <- load_real_dataset(dataset, max_cells=max_cells, seed=1)
gm   <- as(real$grna_mat, "CsparseMatrix")
GN   <- rownames(gm); CN <- colnames(gm)

## MOI from RAW data (method-independent adaptive threshold; == de_native.R)
.mx <- sparseMatrixStats::colMaxs(gm); .C <- ncol(gm)
.gv <- gm@i+1L; .cv <- rep.int(seq_len(.C), diff(gm@p)); .xv <- gm@x
.above <- .xv >= pmax(3, 0.2*.mx)[.cv]; .gpc <- tabulate(.cv[.above], nbins=.C)
moi <- if(median(.gpc) >= 2) "high" else "low"
cat(sprintf("[%s] %d guides x %d cells | MOI=%s (median %d guides/cell)\n",
            dataset, nrow(gm), ncol(gm), moi, median(.gpc)))
if(moi=="high") GRID <- GRID[GRID$qc!="singleton",]   # singleton QC meaningless at high MOI

## ============================ fishash+ fit (once, cached) ====================
fit_fp <- file.path(out_dir,"fit.rds")
if(file.exists(fit_fp)){ fit <- readRDS(fit_fp); cat("  loaded cached fit\n")
} else {
  t0 <- Sys.time()
  fit <- cassign(gm, q=0.05, cell_margin="ambient", tail="hyper", fdr="GS")
  cat(sprintf("  fishash+ GS fit: %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
  saveRDS(fit, fit_fp)
}

## ============================ response + target df + cov (once) ==============
if(real$mode=="survey"){
  tdf  <- read.csv(file.path(spec$dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)[, c("grna_id","grna_target")]
  ## MEMORY FIX: subset the response to on-target + <=3000 random genes (was: full 36k-gene
  ## matrix held in memory and forked per parallel worker -> the OOM). Mirrors the ondisc path.
  resp_full <- as(real$resp_mem, "CsparseMatrix"); all_genes <- rownames(resp_full)
  on_genes  <- intersect(unique(tdf$grna_target), all_genes)
  set.seed(1); rand_genes <- sample(setdiff(all_genes, on_genes), min(3000, length(all_genes)-length(on_genes)))
  cat(sprintf("  response subset: %d genes (%d on-target + random) [was %d]\n",
              length(unique(c(on_genes,rand_genes))), length(on_genes), nrow(resp_full)))
  resp <- resp_full[unique(c(on_genes, rand_genes)), , drop=FALSE]
  rm(resp_full); gc(FALSE)
} else {
  scep_obj <- real$obj; tdf <- scep_obj@grna_target_data_frame[, c("grna_id","grna_target")]
  resp_odm <- scep_obj@response_matrix[[1]]; all_genes <- rownames(resp_odm)
  on_genes <- intersect(unique(tdf$grna_target), all_genes)
  set.seed(1); rand_genes <- sample(setdiff(all_genes, on_genes), min(3000, length(all_genes)-length(on_genes)))
  cat(sprintf("  pulling %d genes (%d on-target + random) for %d cells\n",
              length(unique(c(on_genes,rand_genes))), length(on_genes), length(real$cells_idx)))
  resp <- .pull_response(resp_odm, unique(c(on_genes, rand_genes)), real$cells_idx)
}
cov_full <- data.frame(g_umis = Matrix::colSums(gm))   # method-independent grna covariate

## ============================ build A for an operating point =================
build_A <- function(corr, q){
  A <- if(corr=="perguideBH") apply_perguide(fit, q) else apply_corr(fit, corr, q)
  A <- as(A, "CsparseMatrix"); dimnames(A) <- list(GN, CN); A
}

## ============================ sceptre DE (faithful to de_native.R) ===========
run_de <- function(A, keep, op_dir){
  dir.create(op_dir, recursive=TRUE, showWarnings=FALSE)
  # cell subset (singleton QC drops multiplet cells; else keep all)
  resp_k <- resp[, keep, drop=FALSE]; A_k <- A[, keep, drop=FALSE]; cov_k <- cov_full[keep,,drop=FALSE]
  common_g <- intersect(rownames(A_k), tdf$grna_id)
  A_k <- A_k[common_g,,drop=FALSE]; tdf_k <- tdf[match(common_g, tdf$grna_id),]
  A_k <- as(A_k,"CsparseMatrix")*1
  stopifnot(ncol(A_k)==ncol(resp_k))
  imp <- function() import_data(response_matrix=as(resp_k,"dgCMatrix"), grna_matrix=A_k,
                                grna_target_data_frame=tdf_k, moi=moi, extra_covariates=cov_k)
  pcp <- tryCatch(construct_positive_control_pairs(imp()), error=function(e) NULL)
  pcp0 <- pcp %||% data.frame(grna_target=character(0), response_id=character(0))
  scep <- NULL
  for(fm in list(~ log(response_n_umis) + log(response_n_nonzero) + log(g_umis+1),
                 ~ log(response_n_umis) + log(g_umis+1), ~ log(response_n_umis))){
    scep <- tryCatch(set_analysis_parameters(imp(), positive_control_pairs=pcp0,
              formula_object=fm, control_group="default"),
              error=function(e){cat("   formula failed:",conditionMessage(e),"\n"); NULL})
    if(!is.null(scep)) break
  }
  if(is.null(scep)) return(NULL)
  scep <- assign_grnas(scep, method="thresholding", threshold=1, print_progress=FALSE)
  scep <- tryCatch(run_qc(scep, n_nonzero_trt_thresh=0, n_nonzero_cntrl_thresh=0,
            response_n_umis_range=c(0,1), response_n_nonzero_range=c(0,1), p_mito_threshold=1),
            error=function(e) scep)
  par <- ncpus>=2
  nt_ids <- tdf_k$grna_id[tdf_k$grna_target=="non-targeting"]
  gpt <- table(tdf_k$grna_target[tdf_k$grna_target!="non-targeting"])
  cgs <- max(1, min(round(median(as.numeric(gpt))), length(nt_ids)-1))
  scep <- tryCatch(run_calibration_check(scep, output_amount=2, print_progress=FALSE,
            n_calibration_pairs=NCAL, calibration_group_size=cgs,
            parallel=par, n_processors=if(par) ncpus else "auto"),
            error=function(e){cat("   CAL ERR:",conditionMessage(e),"\n"); scep})
  cal <- tryCatch(get_result(scep,"run_calibration_check"), error=function(e) NULL)
  if(!is.null(cal)) write.csv(cal, file.path(op_dir,"calibration.csv"), row.names=FALSE)
  pow <- NULL
  if(!is.null(pcp) && nrow(pcp)>0){
    scep <- tryCatch(run_power_check(scep, output_amount=2, print_progress=FALSE),
              error=function(e){cat("   POW ERR:",conditionMessage(e),"\n"); scep})
    pow <- tryCatch(get_result(scep,"run_power_check"), error=function(e) NULL)
    if(!is.null(pow)) write.csv(pow, file.path(op_dir,"power.csv"), row.names=FALSE)
  }
  cp <- if(!is.null(cal)) cal$p_value[!is.na(cal$p_value)] else numeric(0)
  pp <- if(!is.null(pow)) pow$p_value[!is.na(pow$p_value)] else numeric(0)
  list(n_cells=ncol(A_k), n_calls=sum(A_k>0), cal_pairs=length(cp),
       t1e_05=if(length(cp))mean(cp<0.05)else NA, t1e_01=if(length(cp))mean(cp<0.01)else NA,
       ks=if(length(cp)>2) as.numeric(ks.test(cp,"punif")$statistic) else NA,
       pow_pairs=length(pp), pow_frac_sig=if(length(pp))mean(pp<0.05)else NA,
       pow_med_p=if(length(pp))median(pp)else NA)
}

## ============================ operating-point loop (checkpointed) ============
sum_fp <- file.path(out_dir,"summaries.csv")
done <- if(file.exists(sum_fp)) read.csv(sum_fp, stringsAsFactors=FALSE) else NULL
for(r in seq_len(nrow(GRID))){
  op <- GRID[r,]; key <- sprintf("%s_q%.2f_%s", op$corr, op$q, op$qc)
  if(!is.null(done) && any(done$op==key)){ cat("  skip done:", key,"\n"); next }
  cat(sprintf("\n---- [%s] %s ----\n", dataset, key)); t0 <- Sys.time()
  A <- build_A(op$corr, op$q)
  keep <- if(op$qc=="singleton") (Matrix::colSums(A)<=1) else rep(TRUE, ncol(A))
  res <- tryCatch(run_de(A, keep, file.path(out_dir, key)),
                  error=function(e){cat("  RUN ERR:",conditionMessage(e),"\n"); NULL})
  dt <- as.numeric(difftime(Sys.time(),t0,units="mins"))
  row <- data.frame(dataset=dataset, op=key, corr=op$corr, q=op$q, qc=op$qc, moi=moi,
    n_cells_kept=if(is.null(res)) NA else res$n_cells,
    n_calls=if(is.null(res)) NA else res$n_calls,
    cal_pairs=if(is.null(res)) NA else res$cal_pairs,
    t1e_05=if(is.null(res)) NA else res$t1e_05, t1e_01=if(is.null(res)) NA else res$t1e_01,
    ks=if(is.null(res)) NA else res$ks,
    pow_pairs=if(is.null(res)) NA else res$pow_pairs,
    pow_frac_sig=if(is.null(res)) NA else res$pow_frac_sig,
    pow_med_p=if(is.null(res)) NA else res$pow_med_p, mins=round(dt,1))
  done <- rbind(done, row); write.csv(done, sum_fp, row.names=FALSE)
  cat(sprintf("  %s: calls=%s t1e05=%s ks=%s pow_sig=%s (%.1f min)\n", key,
      row$n_calls, signif(row$t1e_05,3), signif(row$ks,3), signif(row$pow_frac_sig,3), dt))
}
cat("\n==== DONE", dataset, "====\n"); print(done[,c("op","n_calls","t1e_05","ks","pow_frac_sig","pow_med_p","mins")], row.names=FALSE)
