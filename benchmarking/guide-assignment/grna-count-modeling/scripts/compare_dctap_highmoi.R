#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(sceptre);library(Matrix);library(ggplot2)})
self<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
source(file.path(dirname(normalizePath(self)),'fishash_repro_lib.R'))
source(file.path(dirname(normalizePath(self)), 'datasets.R'))
input <- dctap_sceptre_dir("highmoi")
out <- paste0(.get_config_path("LOCAL_BENCHMARKING_DIR"),
              "association/pos-control/dctap-highmoi-mixture-fishashplus-20261008")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-as(readRDS(file.path(input,'grna_matrix_aligned.rds')),'CsparseMatrix')
y<-as(readRDS(file.path(input,'response_matrix.rds')),'CsparseMatrix')
tdf<-read.csv(file.path(input,'grna_target_data_frame.csv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(x),colnames(y)),identical(rownames(x),tdf$grna_id))
extra<-data.frame(response_n_nonzero_full=as.numeric(colSums(y>0)),
                  response_n_umis_full=as.numeric(colSums(y)),
                  grna_n_nonzero_full=as.numeric(colSums(x>0)),
                  grna_n_umis_full=as.numeric(colSums(x)),row.names=colnames(x))
formula<-~log(response_n_nonzero_full+1)+log(response_n_umis_full+1)+
  log(grna_n_nonzero_full+1)+log(grna_n_umis_full+1)
RNGkind("L'Ecuyer-CMRG");set.seed(20261008)
mixfp<-file.path(out,'assignments_mixture.rds')
if(!file.exists(mixfp)){
  cat('Fitting mixture:',dim(x),format(Sys.time()),'\n')
  s<-import_data(y,x,tdf[,c('grna_id','grna_target')],moi='high',extra_covariates=extra)
  s<-set_analysis_parameters(s,formula_object=formula,resampling_mechanism='crt')
  s<-assign_grnas(s,method='mixture',parallel=TRUE,n_processors=4,
     log_dir=out,formula_object=formula,n_em_rep=5L,n_nonzero_cells_cutoff=10L,
     backup_threshold=5L,probability_threshold=.8)
  m<-as(get_grna_assignments(s),'CsparseMatrix')
  m<-m[match(rownames(x),rownames(m)),]
  colnames(m)<-colnames(x)
  stopifnot(identical(dimnames(m),dimnames(x)),sum(m*(x>0))==sum(m))
  saveRDS(m,mixfp)
  saveRDS(list(method='mixture',parameters=s@grna_assignment_hyperparameters,
               formula=formula,seed=20261008,session=sessionInfo()),
          file.path(out,'mixture_assignment_metadata.rds'))
  rm(s);invisible(gc())
}else m<-readRDS(mixfp)
plusfp<-file.path(out,'assignments_fishashplus.rds')
if(!file.exists(plusfp)){
  cat('Fitting Poisson Fishash+ at',format(Sys.time()),'\n')
  f<-fishash_plus_poisson(x,q=.05,refit=10,min_count=2,cell_margin='ambient',exclude_empty=TRUE)
  saveRDS(f,plusfp)
  extended<-fishash_plus_poisson(x,q=.05,refit=30,min_count=2,cell_margin='ambient',exclude_empty=TRUE)
  delta<-sum(f!=extended)
  cat('Fishash+ extended-refit changed calls:',delta,'\n')
  saveRDS(list(refit=10,validation_refit=30,changed_calls=delta),file.path(out,'fishashplus_stability.rds'))
  if(delta>0)saveRDS(extended,file.path(out,'assignments_fishashplus_refit30.rds'))
}else f<-readRDS(plusfp)
shared<-sum(m*f)
overlap<-data.frame(mixture=sum(m),fishashplus=sum(f),shared=shared,
                   mixture_only=sum(m)-shared,fishashplus_only=sum(f)-shared)
write.csv(overlap,file.path(out,'call_overlap.csv'),row.names=FALSE);print(overlap)
pcgenes<-unique(tdf$grna_target[tdf$grna_target %in% rownames(y)])
pc<-data.frame(grna_target=pcgenes,response_id=pcgenes)
for(scope in c('promoter_union','individual_guides'))for(arm in c('mixture','fishashplus')){
  resultfp<-file.path(out,paste0(scope,'_',arm,'_results.rds'))
  if(file.exists(resultfp))next
  calls<-if(arm=='mixture')m else f
  keep<-if(scope=='promoter_union')grepl('TSS',tdf$grna_id)|tdf$grna_target=='non-targeting' else rep(TRUE,nrow(tdf))
  s<-import_data(y,as(calls[keep,],'dMatrix'),tdf[keep,c('grna_id','grna_target')],
                  moi='high',extra_covariates=extra)
  s<-set_analysis_parameters(s,positive_control_pairs=pc,formula_object=formula,
    grna_integration_strategy=if(scope=='promoter_union')'union' else 'singleton',
    side='both',resampling_mechanism='crt',resampling_approximation='skew_normal',
    control_group='complement',multiple_testing_method='BH',multiple_testing_alpha=.1)
  s<-assign_grnas(s,method='thresholding',threshold=1)
  s<-run_qc(s,n_nonzero_trt_thresh=0,n_nonzero_cntrl_thresh=0,
    response_n_umis_range=c(0,1),response_n_nonzero_range=c(0,1),p_mito_threshold=1)
  stopifnot(length(s@cells_in_use)==ncol(x))
  cat('Running',scope,arm,'pairs',s@n_ok_positive_control_pairs,format(Sys.time()),'\n')
  set.seed(20261008)
  s<-run_power_check(s,output_amount=2,parallel=FALSE)
  r<-as.data.frame(get_result(s,'run_power_check'))
  stopifnot(all(is.finite(r$p_value)),all(is.finite(r$log_2_fold_change)),all(r$pass_qc))
  r$bh_adjusted_p<-p.adjust(r$p_value,'BH')
  saveRDS(r,resultfp);write.csv(r,sub('rds$','csv',resultfp),row.names=FALSE)
  print(r)
  rm(s);invisible(gc())
}
saveRDS(list(input=input,dim_counts=dim(x),dim_expression=dim(y),formula=formula,
             covariates=extra,target_df=tdf,seed=20261008,session=sessionInfo(),
             source_md5=tools::md5sum(file.path(input,c('grna_matrix_aligned.rds','response_matrix.rds','grna_target_data_frame.csv')))),
        file.path(out,'analysis_metadata.rds'))
cat('Finished',format(Sys.time()),'\n')
