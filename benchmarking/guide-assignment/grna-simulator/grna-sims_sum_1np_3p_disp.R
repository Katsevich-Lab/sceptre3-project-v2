rm(list=ls())

source("~/.Rprofile")
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "../../association/data-preparation/utils_data_prep.R"))
source(file.path(script_dir, "/sims-for-paper.R"))



path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data",
  "replogle-rd7",
  "sceptre-pipeline"
)
scep_repl <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)


params_np_highvar <- list(
  prob_endo = 0.002,
  mu_endo = 10,
  theta_endo = 1,
  mu_exo = 0.01
)

params_p_highdisp <- list(
  prob_pert = 0.01,
  mu_pert = 1000,
  theta_pert = .1
)

params_p_meddisp <- list(
  prob_pert = 0.01,
  mu_pert = 1000,
  theta_pert = 1
)

params_p_lowdisp <- list(
  prob_pert = 0.01,
  mu_pert = 1000,
  theta_pert = 10
)

params_list <- list(
  NPhighvar_Plowdisp = c(params_np_highvar, params_p_lowdisp),
  NPhighvar_Pmeddisp = c(params_np_highvar, params_p_meddisp),
  NPhighvar_Phighdisp = c(params_np_highvar, params_p_highdisp)
)

dataset_name = "sims_sum_1np_3p_disp"
make_grna_simulation_sum_process(
  dataset_name = dataset_name,
  params_list = params_list,
  num_cells = 50000,
  num_guides_per_param = 100,
  scep = scep_repl,
  seed = 654,
  use_response_covariates = TRUE
)


# load the grna matrix and odm, make side-by-side 
path_to_new_data = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data", dataset_name)
grna_mat = file.path(path_to_new_data,  "sceptre/grna_matrix.rds") |>
  read_rds()
true_perts = file.path(path_to_new_data,  "true_pert_matrix.rds") |>
  read_rds()

grna_odm = ondisc::initialize_odm_from_backing_file(
  file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
            "guide_assignment/input_data/replogle-rd7/sceptre-pipeline/grna.odm")
)

i=201
j=1
plot_umi_histogram_real_vs_sim(
  umis_real = grna_odm[j,],
  umis_sim = grna_mat[i,],
  is_pert = true_perts[i,],
  title = rownames(grna_mat)[i]
)



