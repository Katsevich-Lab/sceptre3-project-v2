this_file <- normalizePath(sys.frame(1)$ofile)
script_dir <- dirname(this_file)
source(file.path(script_dir, "convert_odm_to_h5ad.R"))
source("~/.Rprofile")

make_data <- function(data_dir, write_crispat_fp=NULL, write_pertpy_fp=NULL, write_sceptre_fp=NULL,
                      write_cleanser_fp=NULL) {
  make_h5ad <- !is.null(write_crispat_fp) || !is.null(write_pertpy_fp)  # only need to do one, and copy the other if both needed
  
  # loading grna_matrix into R
  grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(data_dir, "grna.odm"))
  cell_names <- paste0("CELL_", 1:ncol(grna_odm))
  # grna_mat <- odm_to_R(grna_odm)  # sparse matrix
  # colnames(grna_mat) <- cell_names
  
  if(!is.null(write_cleanser_fp)) {
    Matrix::writeMM(grna_mat, write_cleanser_fp)
  } 
  if(make_h5ad) {
    library(reticulate)
    conda_create("r-anndata", packages = c("python=3.12"))
    use_condaenv("r-anndata", required = TRUE)
    py_install(c("numpy", "scipy", "h5py", "anndata"), 
               envname = "r-anndata", pip = TRUE)
    py_config()
    
    # if only doing one, just write that one. Otherwise write pertpy first
    # and then cp that to make crispat's
    if(!is.null(write_crispat_fp) && is.null(write_pertpy_fp)) {
      R_to_h5ad(grna_mat, write_crispat_fp)
    } else if(is.null(write_crispat_fp) && !is.null(write_pertpy_fp)) {
      R_to_h5ad(grna_mat, write_pertpy_fp)
    } else {
      R_to_h5ad(grna_mat, write_pertpy_fp)
      cmd <- paste("cp", write_pertpy_fp, write_crispat_fp)
      system(cmd)
    }
  }
  
  if(!is.null(write_sceptre_fp)) {
    # saveRDS(grna_mat, file.path(write_sceptre_fp, "grna_matrix.rds"))
    
    # rm(grna_mat)
    
    #response_odm <- ondisc::initialize_odm_from_backing_file(file.path(data_dir, "response.odm"))
    print("loading...")
    #response_mat <- odm_to_R(response_odm)
    print("loaded")
    #colnames(response_mat) <- cell_names
    #saveRDS(response_mat, file.path(write_sceptre_fp, "response_matrix.rds"))
    #rm(response_mat)
    
    scep <- readRDS(file.path(data_dir, "sceptre_object.rds"))
    write.csv(scep@grna_target_data_frame, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names=FALSE)
  }
}

# paths
gasperini_base_fp <- .get_config_path("LOCAL_GASPERINI_2019_V3_DATA_DIR")

data_dir <- file.path(gasperini_base_fp, "at-scale/processed")

write_cleanser_fp <- NULL # file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/grna_counts_cleanser.mtx")
write_crispat_fp <- NULL # file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/grna_counts_crispat.h5ad")
write_pertpy_fp <- NULL # file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/grna_counts_pertpy.h5ad")
write_sceptre_fp <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/sceptre_input")


make_data(data_dir = data_dir, write_cleanser_fp = write_cleanser_fp, write_crispat_fp = write_crispat_fp, 
          write_pertpy_fp = write_pertpy_fp, write_sceptre_fp = write_sceptre_fp)
