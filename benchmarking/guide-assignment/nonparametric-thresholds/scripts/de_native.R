## Self-contained per-(dataset, method) DE runner using sceptre's NATIVE
## run_calibration_check (negative-control calibration) + run_power_check
## (positive-control power). This is sceptre's intended calibration API and
## picks the correct control group by MOI (nt_cells for low, complement for high).
## Saves calibration + power p-values for downstream analysis.
##
## Usage: Rscript scripts/de_native.R <dataset> <method> <max_cells> <out_dir>
suppressPackageStartupMessages({library(Matrix); library(sceptre)})
source(file.path(getwd(),"scripts","sim_de_lib.R"))
source(file.path(getwd(),"scripts","contingency_method.R"))
`%||%` <- function(a,b) if(is.null(a)) b else a

args <- commandArgs(trailingOnly=TRUE)
dataset <- args[1]; method <- args[2]
max_cells <- if(length(args)>=3) as.integer(args[3]) else 30000
out_dir <- args[4]; dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
ncpus <- as.integer(Sys.getenv("NCPUS","0"))

spec <- REAL_DATASETS[[dataset]]
real <- load_real_dataset(dataset, max_cells=max_cells, seed=1)
## MOI determined from RAW data (method-independent, registry flag is unreliable):
## per-cell adaptive threshold count >= max(3, 0.2*cellmax); high if median >=2.
.gm <- as(real$grna_mat,"CsparseMatrix"); .C <- ncol(.gm)
.mx <- sparseMatrixStats::colMaxs(.gm)
.gv <- .gm@i+1L; .cv <- rep.int(seq_len(.C), diff(.gm@p)); .xv <- .gm@x
.thr <- pmax(3, 0.2*.mx); .above <- .xv >= .thr[.cv]
.gpc <- tabulate(.cv[.above], nbins=.C)
moi <- if(median(.gpc) >= 2) "high" else "low"
cat(sprintf("  measured guides/cell median=%d mean=%.2f -> MOI=%s\n", median(.gpc), mean(.gpc), moi))

## assignment matrix for this method
run_external <- function(gmat, method){
  # export guides x cells matrix, run crispat/cleanser via python, read calls
  gmat <- as(gmat, "CsparseMatrix")
  td <- tempfile("ext_"); dir.create(td)
  Matrix::writeMM(gmat, file.path(td,"m.mtx"))
  writeLines(rownames(gmat), file.path(td,"guides.txt"))
  writeLines(colnames(gmat), file.path(td,"cells.txt"))
  py <- if(method=="cleanser" && nzchar(Sys.getenv("CLEANSER_PY"))) Sys.getenv("CLEANSER_PY")
        else file.path(getwd(), ".venv_crispat","bin","python")
  ec <- system2(py, c(file.path(getwd(),"scripts","run_external_assign.py"), method,
                      file.path(td,"m.mtx"), file.path(td,"guides.txt"),
                      file.path(td,"cells.txt"), file.path(td,"calls.csv"),
                      Sys.getenv("EXT_NITER","500")),
                stdout="", stderr="")
  calls <- if(file.exists(file.path(td,"calls.csv"))) read.csv(file.path(td,"calls.csv"), stringsAsFactors=FALSE) else data.frame(guide=character(0),cell=character(0))
  gi <- match(calls$guide, rownames(gmat)); ci <- match(calls$cell, colnames(gmat))
  ok <- !is.na(gi) & !is.na(ci)
  A <- Matrix::sparseMatrix(i=gi[ok], j=ci[ok], x=TRUE, dims=dim(gmat), dimnames=dimnames(gmat))
  unlink(td, recursive=TRUE); A
}
if(method=="depth_fix_nb"){
  A <- as(contingency_assign(real$grna_mat, cell_margin="ambient", tail="nb", fdr="GS")$assigned, "lgCMatrix")
} else if(method %in% c("crispat","cleanser")){
  A <- run_external(real$grna_mat, method)
} else {
  A <- run_methods_on_real(real, methods=method, verbose=TRUE)[[method]]
}
A <- as(A, "CsparseMatrix") * 1                              # logical -> 0/1 numeric dgCMatrix
cat(sprintf("[%s/%s] %d calls over %d guides x %d cells\n", dataset, method, sum(A>0), nrow(A), ncol(A)))

## response + target df + covariates
if(real$mode=="survey"){
  resp <- real$resp_mem
  tdf  <- read.csv(file.path(spec$dir,"grna_target_data_frame.csv"), stringsAsFactors=FALSE)
  tdf  <- tdf[, c("grna_id","grna_target")]
} else {
  ## ondisc (gasperini/replogle): pull on-target genes + random genes into memory
  scep_obj <- real$obj
  tdf  <- scep_obj@grna_target_data_frame[, c("grna_id","grna_target")]
  resp_odm <- scep_obj@response_matrix[[1]]
  all_genes <- rownames(resp_odm)
  on_genes  <- intersect(unique(tdf$grna_target), all_genes)
  set.seed(1); rand_genes <- sample(setdiff(all_genes, on_genes), min(3000, length(all_genes)-length(on_genes)))
  genes <- unique(c(on_genes, rand_genes))
  cat(sprintf("  pulling %d genes (%d on-target + random) from ondisc for %d cells\n",
              length(genes), length(on_genes), length(real$cells_idx)))
  # .pull_response(response_odm, gene_ids, cell_idx) from sim_de_lib
  resp <- .pull_response(resp_odm, genes, real$cells_idx)
}
# keep only guides present in both A and tdf
common_g <- intersect(rownames(A), tdf$grna_id)
A <- A[common_g,,drop=FALSE]; tdf <- tdf[match(common_g, tdf$grna_id),]
stopifnot(ncol(A)==ncol(resp))

# grna covariates from RAW counts (method-independent; the binary assignment A would
# make grna_n_umis==grna_n_nonzero -> collinear formula). Response covariates are
# auto-computed by sceptre from the GEX matrix.
rawg <- as(real$grna_mat,"CsparseMatrix")
cov <- data.frame(g_umis = Matrix::colSums(rawg))
scep <- import_data(response_matrix=as(resp,"dgCMatrix"), grna_matrix=A,
                    grna_target_data_frame=tdf, moi=moi, extra_covariates=cov)
pcp <- tryCatch(construct_positive_control_pairs(scep), error=function(e) NULL)
cat("  positive control pairs:", if(is.null(pcp)) 0 else nrow(pcp), "\n")
pcp0 <- pcp %||% data.frame(grna_target=character(0), response_id=character(0))
fmla_try <- list(~ log(response_n_umis) + log(response_n_nonzero) + log(g_umis+1),
                 ~ log(response_n_umis) + log(g_umis+1),
                 ~ log(response_n_umis))
scep <- NULL
for(fm in fmla_try){
  scep <- tryCatch(set_analysis_parameters(sc0 <- import_data(response_matrix=as(resp,"dgCMatrix"),
              grna_matrix=A, grna_target_data_frame=tdf, moi=moi, extra_covariates=cov),
              positive_control_pairs=pcp0, formula_object=fm, control_group="default"),
            error=function(e){cat("  formula failed:",conditionMessage(e),"\n"); NULL})
  if(!is.null(scep)){ cat("  formula:", deparse(fm), "\n"); break }
}
if(is.null(scep)) stop("all formulas failed")
scep <- assign_grnas(scep, method="thresholding", threshold=1, print_progress=FALSE)
scep <- tryCatch(run_qc(scep, n_nonzero_trt_thresh=0, n_nonzero_cntrl_thresh=0,
              response_n_umis_range=c(0,1), response_n_nonzero_range=c(0,1), p_mito_threshold=1),
          error=function(e){cat("  (run_qc skipped:",conditionMessage(e),")\n"); scep})

par <- ncpus>=2
## calibration group size: mimic real targets (median guides per targeting target),
## capped by available NT guides.
nt_ids <- tdf$grna_id[tdf$grna_target=="non-targeting"]
gpt <- table(tdf$grna_target[tdf$grna_target!="non-targeting"])
cgs <- max(1, min(round(median(as.numeric(gpt))), length(nt_ids)-1))
cat(sprintf("  NT guides=%d, calibration_group_size=%d\n", length(nt_ids), cgs))
## CALIBRATION (negative control)
scep <- tryCatch(run_calibration_check(scep, output_amount=2, print_progress=FALSE,
                   n_calibration_pairs=2000, calibration_group_size=cgs,
                   parallel=par, n_processors=if(par) ncpus else "auto"),
          error=function(e){cat("  CAL ERR:",conditionMessage(e),"\n"); scep})
cal <- tryCatch(get_result(scep, "run_calibration_check"), error=function(e) NULL)
if(!is.null(cal)) write.csv(cal, file.path(out_dir,"calibration.csv"), row.names=FALSE)

## POWER (positive control)
if(!is.null(pcp) && nrow(pcp)>0){
  scep <- tryCatch(run_power_check(scep, output_amount=2, print_progress=FALSE),
            error=function(e){cat("  POW ERR:",conditionMessage(e),"\n"); scep})
  pow <- tryCatch(get_result(scep, "run_power_check"), error=function(e) NULL)
  if(!is.null(pow)) write.csv(pow, file.path(out_dir,"power.csv"), row.names=FALSE)
}

## quick summary line
smry <- function(){
  cl <- if(file.exists(file.path(out_dir,"calibration.csv"))) read.csv(file.path(out_dir,"calibration.csv")) else NULL
  pw <- if(file.exists(file.path(out_dir,"power.csv"))) read.csv(file.path(out_dir,"power.csv")) else NULL
  cp <- if(!is.null(cl)) cl$p_value[!is.na(cl$p_value)] else numeric(0)
  pp <- if(!is.null(pw)) pw$p_value[!is.na(pw$p_value)] else numeric(0)
  data.frame(dataset=dataset, method=method, moi=moi, n_calls=sum(A>0),
    cal_pairs=length(cp), t1e_05=mean(cp<0.05), t1e_01=mean(cp<0.01),
    ks=if(length(cp)>2) as.numeric(ks.test(cp,"punif")$statistic) else NA,
    pow_pairs=length(pp), pow_med_p=if(length(pp)) median(pp) else NA,
    pow_frac_sig=if(length(pp)) mean(pp<0.05) else NA)
}
row <- smry(); write.csv(row, file.path(out_dir,"summary.csv"), row.names=FALSE)
cat("SUMMARY: "); print(row, row.names=FALSE)
