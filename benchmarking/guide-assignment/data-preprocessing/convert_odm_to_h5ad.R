



#' 
#' conda_create("r-anndata", packages = c("python=3.12"))
#' use_condaenv("r-anndata", required = TRUE)
#' py_install(c("numpy", "scipy", "h5py", "anndata"), #"zarr")
#'            envname = "r-anndata", pip = TRUE)
#' py_config()
#' 
#' 
#' 
#' #' write_h5ad_from_odm_in_memory
#' #' 
#' #' Note that this function loads the entire odm file into memory.
#' #'
#' #' @param odm : an odm object as defined by the ondisc package
#' #' @param path_to_write : a string giving the path to write this object
#' #' @param chunk_size : the rows of `odm` are initially grouped into batches of this size
#' 
#' write_h5ad_from_odm_in_memory <- function(odm, path_to_write, chunk_size = 10000) {
#'   if(!inherits(odm, "odm")) {
#'     stop("`odm` must be an object of class odm as defined in the ondisc package.")
#'   }
#'   num_rows <- nrow(odm)
#'   num_cols <- ncol(odm)
#'   num_chunks <- ceiling(num_rows / chunk_size)
#'   chunk_list <- vector("list", num_chunks)
#'   for(i in 1:num_chunks) {
#'     idx_start <- (i - 1) * chunk_size + 1
#'     idx_end   <- min(i * chunk_size, num_rows)
#'     chunk_list[[i]] <- lapply(idx_start:idx_end, function(j) odm[j,]) |>
#'       do.call(what = rbind) |>
#'       t() |>  # h5ad needs cells x genes so we transpose now
#'       as("CsparseMatrix")
#'   }
#'   out <- do.call(cbind, chunk_list)
#'   colnames(out) <- rownames(odm) # we have transposed
#'   
#'   out_py  <- r_to_py(out)                     # dgCMatrix -> scipy.sparse.csc_matrix
#'   adata <- anndata::AnnData(X = out_py)
#'   adata$write_h5ad(path_to_write)
#' } 

row_to_sparse_matrix <- function(values) {
  nonzero <- which(values != 0)
  Matrix::sparseMatrix(
    i = rep(1L, length(nonzero)),  # it's all one row
    j = nonzero,          
    x = values[nonzero],                
    dims = c(1L, length(values)),
    repr = "R"
  )
}

odm_to_R <- function(odm, num_rows=NULL) {
  if(!inherits(odm, "odm")) {
    stop("`odm` must be an object of class odm as defined in the ondisc package.")
  }
  num_rows_to_read <- ifelse(is.null(num_rows), nrow(odm), num_rows)
  out <- lapply(1:num_rows_to_read, function(i) odm[i,] |> row_to_sparse_matrix()) |>
    do.call(what = rbind)
  rownames(out) <- rownames(odm)[1:num_rows_to_read]
  return(out)
}



R_to_h5ad <- function(mat, path_to_write) {
  library(reticulate)
  library(SingleCellExperiment)
  library(zellkonverter)
  
  if(is.null(colnames(mat))) {
    colnames(mat) <- paste0("CELL_", seq_len(ncol(mat)))
  }
  
  # Wrap & write (zellkonverter streams blocks to .h5ad)
  sce <- SingleCellExperiment(assays = list(counts = mat))
  writeH5AD(sce, path_to_write)
}



#' 
#' library(reticulate)
#' 
#' conda_create("r-anndata", packages = c("python=3.12"))
#' use_condaenv("r-anndata", required = TRUE)
#' py_install(c("numpy", "scipy", "h5py", "anndata"), #"zarr")
#'            envname = "r-anndata", pip = TRUE)
#' py_config()
#' 
#' #' write_h5ad_from_odm_in_memory
#' #' 
#' #' Note that this function loads the entire odm file into memory.
#' #'
#' #' @param odm : an odm object as defined by the ondisc package
#' #' @param path_to_write : a string giving the path to write the resulting h5ad
#' #' @param chunk_size : the rows of `odm` are initially grouped into batches of this size
#' write_h5ad_from_odm_in_memory <- function(odm, path_to_write, chunk_size = 100) {
#'   if(!inherits(odm, "odm")) {
#'     stop("`odm` must be an object of class odm as defined in the ondisc package.")
#'   }
#'   suppressPackageStartupMessages({
#'     library(Matrix)
#'     library(SingleCellExperiment)
#'     library(zellkonverter)
#'   })
#'   
#'   # helper func to make each row of the odm into a 1 x num_cells dgRMatrix
#'   row_to_sparse_matrix <- function(values) {
#'     nonzero <- which(values != 0)
#'     Matrix::sparseMatrix(
#'       i = rep(1L, length(nonzero)),  # it's all one row
#'       j = nonzero,          
#'       x = values[nonzero],                
#'       dims = c(1L, length(values)),
#'       repr = "R"
#'     )
#'   }
#'   
#'   n_vars <- nrow(odm)   
#'   n_obs <- ncol(odm)  
#'   n_chunks <- ceiling(n_vars / chunk_size)
#'   chunk_list <- vector("list", n_chunks)
#'   # each chunk gets turned into a dgRMatrix before they all get rbind'd together
#'   for (i in seq_len(n_chunks)) {
#'     idx_start <- (i - 1) * chunk_size + 1
#'     idx_end   <- min(i * chunk_size, n_vars)
#'     chunk_list[[i]] <- lapply(idx_start:idx_end, function(j) odm[j,] |> row_to_sparse_matrix()) |>
#'       do.call(what = rbind)
#'   }
#'   
#'   out_delayed <- do.call(rbind, chunk_list)
#'   rownames(out_delayed) <- rownames(odm)
#'   # ondisc doesn't let these be accessed, so we will make our own
#'   colnames(out_delayed) <- paste0("CELL_", seq_len(n_obs))
#'   
#'   # Wrap & write (zellkonverter streams blocks to .h5ad)
#'   sce <- SingleCellExperiment(assays = list(counts = out_delayed))
#'   writeH5AD(sce, path_to_write)
#' }

# for testing
# my_odm <- ondisc::initialize_odm_from_backing_file("~/katsevich-lab/code/sceptre-examples/nextflow-example-data/response.odm")
# write_h5ad_from_odm_in_memory(my_odm, "~/katsevich-lab/code/sceptre-examples/nextflow-example-data/response.h5ad", 75)


# for real
# gasperini_odm <- ondisc::initialize_odm_from_backing_file("~/katsevich-lab/data/external/gasperini-2019-v3/at-scale/processed/grna.odm")
# dim(gasperini_odm)
# write_h5ad_from_odm_in_memory(gasperini_odm, "~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini/grna.h5ad", 50)
