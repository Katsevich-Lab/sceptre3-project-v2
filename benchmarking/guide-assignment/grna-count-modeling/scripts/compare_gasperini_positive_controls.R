#!/usr/bin/env Rscript
# Paired full-cell SCEPTRE CRT analysis: threshold>=5 versus retained R Fishash+.
# Usage: Rscript compare_gasperini_positive_controls.R <output-directory> [smoke|full] [threshold5,mixture,fishashplus]
suppressPackageStartupMessages({library(Matrix); library(ondisc); library(sceptre)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Supply an output directory; see the usage comment.")
out <- args[1]
mode <- if (length(args) > 1) args[2] else "full"
arms <- if (length(args) > 2) strsplit(args[3],',',fixed=TRUE)[[1]] else c('threshold5','fishashplus')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
base <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
input <- file.path(base,'input_data/gasperini')
plus_fp <- file.path(base,'outputs/gasperini/fishashplus-poisson-20261008/grna_assignment_matrix.rds')
cache <- file.path(out,'paired_inputs.rds')
if (!file.exists(cache)) {
  cat('Preparing common input at',format(Sys.time()),'\n')
  original <- readRDS(file.path(input,'sceptre-pipeline/sceptre_object.rds'))
  response <- initialize_odm_from_backing_file(file.path(input,'sceptre-pipeline/response.odm'))
  target_df <- as.data.frame(original@grna_target_data_frame)
  targets <- unique(target_df$grna_target[target_df$grna_target %in% rownames(response)])
  rows <- which(target_df$grna_target %in% targets)
  counts <- as(readRDS(file.path(input,'sceptre/grna_matrix.rds')),'CsparseMatrix')
  plus <- readRDS(plus_fp)
  stopifnot(identical(dimnames(counts),dimnames(plus)),
            identical(rownames(counts),target_df$grna_id),
            ncol(counts)==ncol(response))
  ii <- jj <- xx <- vector('list',length(targets))
  for (k in seq_along(targets)) {
    y <- response[targets[k], ]; nz <- which(y>0)
    ii[[k]]<-rep.int(k,length(nz));jj[[k]]<-nz;xx[[k]]<-y[nz]
  }
  expr <- sparseMatrix(i=unlist(ii),j=unlist(jj),x=unlist(xx),
                       dims=c(length(targets),ncol(counts)),
                       dimnames=list(targets,colnames(counts)))
  covariates <- original@covariate_data_frame
  full_names <- c('response_n_nonzero','response_n_umis','grna_n_nonzero','grna_n_umis')
  extra <- covariates[,c(full_names,'prep_batch')]
  names(extra)[seq_along(full_names)]<-paste0(full_names,'_full')
  rownames(extra)<-colnames(counts)
  dat <- list(response=expr,targets=target_df[rows,c('grna_id','grna_target')],
              threshold5=as(counts[rows,]>=5,'dgCMatrix'),
              fishashplus=as(plus[rows,],'dgCMatrix'),extra=extra,
              source_input=input,plus_fp=plus_fp)
  stopifnot(nrow(expr)==377L,nrow(dat$targets)==754L,
            sum(dat$threshold5*dat$fishashplus)==sum(dat$threshold5))
  saveRDS(dat,cache)
  cat('Saved common input:',nrow(expr),'genes;',ncol(expr),'cells\n')
} else dat <- readRDS(cache)
rownames(dat$extra)<-colnames(dat$response)
if ('mixture' %in% arms) {
  mixture_fp <- file.path(base,'outputs/gasperini/fishashplus-poisson-20261008/mixture_assignment_matrix_aligned.rds')
  mixture <- readRDS(mixture_fp)
  stopifnot(identical(colnames(mixture),colnames(dat$response)),
            all(dat$targets$grna_id %in% rownames(mixture)))
  dat$mixture <- as(mixture[match(dat$targets$grna_id,rownames(mixture)),], 'dMatrix')
  stopifnot(identical(rownames(dat$mixture),dat$targets$grna_id))
  rm(mixture)
}
formula <- ~log(response_n_nonzero_full+1)+log(response_n_umis_full+1)+
  log(grna_n_umis_full+1)+log(grna_n_nonzero_full+1)+prep_batch
RNGkind("L'Ecuyer-CMRG")
for (arm in arms) {
  result_fp <- file.path(out,paste0(mode,'_',arm,'_results.rds'))
  if (file.exists(result_fp)) {cat('Using completed result',result_fp,'\n');next}
  cat('Starting',mode,arm,'at',format(Sys.time()),'\n')
  s <- import_data(response_matrix=dat$response,grna_matrix=dat[[arm]],
                   grna_target_data_frame=dat$targets,moi='high',extra_covariates=dat$extra)
  pc <- construct_positive_control_pairs(s)
  if (mode=='smoke')pc<-pc[1,,drop=FALSE]
  s <- set_analysis_parameters(s,positive_control_pairs=pc,formula_object=formula,
    side='both',grna_integration_strategy='union',resampling_mechanism='crt',
    resampling_approximation='skew_normal',control_group='complement',
    multiple_testing_method='BH',multiple_testing_alpha=.1)
  s <- assign_grnas(s,method='thresholding',threshold=1)
  # Match the repository positive-control runner: no additional cell/pair filtering.
  s <- run_qc(s,n_nonzero_trt_thresh=0,n_nonzero_cntrl_thresh=0,
    response_n_umis_range=c(0,1),response_n_nonzero_range=c(0,1),p_mito_threshold=1)
  stopifnot(length(s@cells_in_use)==ncol(dat$response))
  cat('QC:',s@n_ok_positive_control_pairs,'of',nrow(pc),'pairs;',length(s@cells_in_use),'cells\n')
  set.seed(20261008)
  t0 <- Sys.time()
  log_dir <- file.path(out,paste0(mode,'_',arm,'_logs'))
  dir.create(log_dir,showWarnings=FALSE)
  s <- run_power_check(s,output_amount=2,parallel=(mode=='full'),n_processors=6,
                       print_progress=TRUE,log_dir=log_dir)
  result <- get_result(s,'run_power_check')
  saveRDS(result,result_fp)
  write.csv(result,sub('rds$','csv',result_fp),row.names=FALSE)
  meta <- list(arm=arm,mode=mode,formula=formula,seed=20261008,
    resampling='CRT',approximation='skew_normal',side='both',grouping='union',
    B=c(s@B1,s@B2,s@B3),qc='No additional cell or pair filtering; same full-cell universe',
    n_cells=length(s@cells_in_use),n_pairs=nrow(pc),
    fit_seconds=as.numeric(difftime(Sys.time(),t0,units='secs')),
    package=packageDescription('sceptre'),session=sessionInfo())
  if (arm=='mixture')meta$assignment_source<-mixture_fp
  saveRDS(meta,file.path(out,paste0(mode,'_',arm,'metadata.rds')))
  cat('Completed',arm,'in',meta$fit_seconds,'seconds\n');print(head(result))
  rm(s);invisible(gc())
}
