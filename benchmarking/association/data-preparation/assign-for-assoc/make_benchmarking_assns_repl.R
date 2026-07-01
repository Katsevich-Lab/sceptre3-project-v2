## this script makes the replogle assignments that i will use for
## making my association benchmarking datasets

## it uses maximum assignment
## The Replogle authors use a Pois-Gauss mixture on the log counts, so that's different
## We will have an issue they do not: we will assign lowly expressed cells to their argmax
## while a mixture would always put those in the NP component.
## These assignments need to be paired with explicit removal of these junk cells. 

rm(list=ls())
library(sceptre)
source("~/.Rprofile")

dataset <- "replogle-rd7"

# shared normalizer for the cell-universe fingerprint (see fingerprint contract);
# MUST match the definition in remove_low_umi_cells_replogle.R and the data-prep loader.
stopifnot(requireNamespace("digest", quietly = TRUE))
fingerprint_cells <- function(x) digest::digest(as.integer(x))

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
  paste0("maximum")
)
dir.create(path_to_outputs, showWarnings = TRUE, recursive = TRUE)

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_data, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_data, "response.odm"),
  grna_odm_file_fp = file.path(path_to_data, "grna.odm")
)
if(!scep@low_moi) {
  stop("The loaded sceptre_object has high MOI but this is Replogle and it should be low MOI!")
}

scep = scep |>
  set_analysis_parameters(resampling_mechanism = "permutations", formula_object = ~1) |>
  assign_grnas(method = "maximum")

assn_mat = get_grna_assignments(scep)

## some checks before writing ~~~~~~~~~~~~~~~~

grna_odm <- file.path(path_to_data, "grna.odm") |> ondisc::initialize_odm_from_backing_file()

# check that guides haven't been shuffled
if(!identical(rownames(assn_mat), rownames(grna_odm))) {
  stop("assignment matrix rownames and odm rownames do not match.")
}

# checking that cells haven't been shuffled

# idea: for perturbed cells for a given guide y,
# we should have `y / grna_n_umis` pretty large,
# because we are looking at the cell's max here
dominance <- function(assn) {
  guides_to_check <- seq(1, nrow(grna_odm), by = 50)
  sapply(guides_to_check, function(i) {     
    is_pert <- which(assn[i, ])
    if (length(is_pert) == 0) return(NA_real_)
    umis <- grna_odm[i, ]
    median(umis[is_pert] / scep@covariate_data_frame$grna_n_umis[is_pert])  # ~1 aligned, ~0 shuffled
  })
}

real <- dominance(assn_mat)
shuf <- dominance(assn_mat[, sample(ncol(assn_mat))])       # negative control

cat("median dominance — real:", round(median(real, na.rm=TRUE), 3),
    " shuffled:", round(median(shuf, na.rm=TRUE), 3), "\n")
stopifnot(median(real, na.rm=TRUE) > 10 * (median(shuf, na.rm=TRUE) + .01))


cat("Average MOI:", mean(Matrix::colSums(assn_mat)) |> round(3), "\n")

guides_per_cell <- Matrix::colSums(assn_mat)
stopifnot(all(guides_per_cell <= 1))                  # low-MOI property of maximum
cat("Cells with 0 assigned guides:", sum(guides_per_cell == 0), "\n")

## write ~~~~~~~~~~~~~~~~~~~
saveRDS(assn_mat, file.path(path_to_outputs, "grna_assignment_matrix.rds"))
saveRDS(
  list(method="maximum", n_cells=ncol(assn_mat),
       # order-sensitive hash to confirm that we haven't shuffled cells when we don't have cell names
       cell_fingerprint = fingerprint_cells(scep@covariate_data_frame$grna_n_umis),  
       date = Sys.time()),
  file.path(path_to_outputs, "assignment_metadata.rds")
)
