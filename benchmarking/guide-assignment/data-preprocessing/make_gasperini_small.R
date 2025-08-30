source("~/katsevich-lab/code/sceptre3-project-v2/benchmarking/guide-assignment/data-preprocessing/convert_odm_to_h5ad.R")

# making a small piece of gasperini
grna_odm <- ondisc::initialize_odm_from_backing_file("~/katsevich-lab/data/external/gasperini-2019-v3/at-scale/processed/grna.odm")
num_grnas <- 200
num_cells <- 2000
cell_names <- paste0("CELL_", 1:num_cells)

grna_mat <- odm_to_R(grna_odm, num_rows = num_grnas)
dim(grna_mat)

grna_small <- grna_mat[,1:num_cells]
colnames(grna_small) <- cell_names

# making the h5ad files for pertpy and crispat

# library(reticulate)
# 
# conda_create("r-anndata", packages = c("python=3.12"))
# use_condaenv("r-anndata", required = TRUE)
# py_install(c("numpy", "scipy", "h5py", "anndata"), #"zarr")
#            envname = "r-anndata", pip = TRUE)
# py_config()
# R_to_h5ad(grna_small, "~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/grna_counts_crispat.h5ad")
# R_to_h5ad(grna_small, "~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/grna_counts_pertpy.h5ad")


# making the mtx file for CLEANSER
Matrix::writeMM(grna_small, "~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/grna_counts_cleanser.mtx")


# making sceptre's files
grna_small |>
  as.matrix() |>
  as.data.frame() |>
  `colnames<-`(paste0("CELL_", 1:ncol(grna_small))) |>
  write.csv("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/sceptre_input/grna_matrix.csv")


response_odm <- ondisc::initialize_odm_from_backing_file("~/katsevich-lab/data/external/gasperini-2019-v3/at-scale/processed/response.odm")
response_small <- odm_to_R(response_odm, num_rows = 100)[,1:num_cells]
colnames(response_small) <- cell_names
response_small |>
  as.matrix() |>
  as.data.frame() |>
  `colnames<-`(paste0("CELL_", 1:ncol(grna_small))) |>
  write.csv("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/sceptre_input/response_matrix.csv")

scep <- readRDS("~/katsevich-lab/data/external/gasperini-2019-v3/at-scale/processed/sceptre_object.rds")
write.csv(scep@grna_target_data_frame[1:num_grnas, ], "~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/gasperini_small/sceptre_input/grna_target_data_frame.csv")

response_matrix <- response_small
grna_matrix <- grna_small
