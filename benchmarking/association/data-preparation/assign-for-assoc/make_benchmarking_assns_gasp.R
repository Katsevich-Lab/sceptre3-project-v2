## this script makes the gasperini assignments that i will use for
## making my association benchmarking datasets

## it uses threshold=5, which is what the Gasperini authors use, and
## EDA confirms that this is reasonable.

rm(list=ls())
library(sceptre)
source("~/.Rprofile")

dataset <- "gasperini"
THRESHOLD = 5

path_to_data <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data",
  dataset,
  "sceptre-pipeline"
)

path_to_outputs <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs",
  dataset,
  paste0("thresholding-", THRESHOLD)
)
dir.create(path_to_outputs, showWarnings = TRUE, recursive = TRUE)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)
if(scep@low_moi) {
  stop("The loaded sceptre_object has low MOI but this is Gasperini and it should be high MOI!")
}

scep = scep |>
  set_analysis_parameters(resampling_mechanism = "permutations", formula_object = ~1) |>
  assign_grnas(method = "thresholding", threshold = THRESHOLD)

assn_mat = get_grna_assignments(scep)

## some checks before writing ~~~~~~~~~~~~~~~~

grna_odm <- file.path(path_to_data, "grna.odm") |> ondisc::initialize_odm_from_backing_file()

# check that guides haven't been shuffled
if(!identical(rownames(assn_mat), rownames(grna_odm))) {
  stop("assignment matrix rownames and odm rownames do not match.")
}

# check that cell names agree [Gasperini only]
if(!identical(colnames(assn_mat), colnames(grna_odm))) {
  stop("assignment matrix and grna ODM cell names do not match.")
}

# checking that cells haven't been shuffled
guides_to_check <- seq(1, nrow(grna_odm), by = 100)
guides_agree <- sapply(guides_to_check, function(idx) all((grna_odm[idx,] >= THRESHOLD) == assn_mat[idx,])) |>
  all()
if(!guides_agree) {
  stop("directly computing the thresholding on the ODM leads to different cell-wise results than in the assignment matrix.")
}

cat("Average MOI:", mean(Matrix::colSums(assn_mat)) |> round(2), "\n")

## write ~~~~~~~~~~~~~~~~~~~
saveRDS(assn_mat, file.path(path_to_outputs, "grna_assignment_matrix.rds"))
saveRDS(
  list(method="thresholding", threshold=THRESHOLD, n_cells=ncol(assn_mat), date = Sys.time()),
  file.path(path_to_outputs, "assignment_metadata.rds")
)
