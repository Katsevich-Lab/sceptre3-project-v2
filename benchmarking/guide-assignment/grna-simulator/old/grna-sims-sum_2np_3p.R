rm(list=ls())

source("~/.Rprofile")
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "../../association/data-preparation/utils_data_prep.R"))
source(file.path(script_dir, "/sims-sum-of-three.R"))



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

params_np_lowvar <- list(
  prob_endo = 0.001,
  mu_endo = 5,
  theta_endo = 10,
  mu_exo = 0.01
)

params_p_3 <- list(
  prob_pert = 0.01,
  mu_pert = 1000,
  theta_pert = 5
)

params_p_2 <- list(
  prob_pert = 0.01,
  mu_pert = 750,
  theta_pert = 2.5
)

params_p_1 <- list(
  prob_pert = 0.01,
  mu_pert = 500,
  theta_pert = 1
)

params_list <- list(
  NPlowvar_P1 = c(params_np_lowvar, params_p_1),
  NPlowvar_P2 = c(params_np_lowvar, params_p_2),
  NPlowvar_P3 = c(params_np_lowvar, params_p_3),
  NPhighvar_P1 = c(params_np_highvar, params_p_1),
  NPhighvar_P2 = c(params_np_highvar, params_p_2),
  NPhighvar_P3 = c(params_np_highvar, params_p_3)
)

dataset_name = "sims_sum_2np_3p"
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

i=1
j=1
plot_umi_histogram_real_vs_sim(
  umis_real = grna_odm[j,],
  umis_sim = grna_mat[i,],
  is_pert = true_perts[i,],
  title = rownames(grna_mat)[i]
)



