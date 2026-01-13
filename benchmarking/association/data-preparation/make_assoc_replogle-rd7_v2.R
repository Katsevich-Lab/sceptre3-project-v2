rm(list=ls())
source("~/.Rprofile")

library(Matrix)
library(tidyverse)

source_data <- "replogle-rd7"

input_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data", source_data, "sceptre-pipeline/"
)
output_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/outputs", source_data, "sceptre-pipeline/"
)

grna_odm <- ondisc::initialize_odm_from_backing_file(paste0(input_fp, "grna.odm"))
response_odm <- ondisc::initialize_odm_from_backing_file(paste0(input_fp, "response.odm"))

scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(input_fp, "sceptre_object.rds"),
  response_odm_file_fp = paste0(input_fp, "response.odm"),
  grna_odm_file_fp = paste0(input_fp, "grna.odm")
)

grna_assign_mat <- readRDS(file.path(output_fp, "grna_assignment_matrix.rds"))

# inputs
num_pc_cells = 35000
num_nt_cells = 15000

grna_target_df = scep@grna_target_data_frame
cell_covariates = scep@covariate_data_frame
NT_name <- "non-targeting"

dataset_label <- "medium"
dataset_name <- paste0(source_data, "_", dataset_label) 


# NOTE
# this has low MOI pretty hard-coded in, sicne we only take the 
# cells expressing exactly one gRNA

## 0. determining which cells express exactly one gRNA. 
##    This is the population we subsample from.

num_grnas_expressed_per_cell <- Matrix::colSums(grna_assign_mat) # potentially big operation
expresses_exactly_one_grna <- which(num_grnas_expressed_per_cell == 1)
grna_assign_mat_to_use <- grna_assign_mat[,expresses_exactly_one_grna]

## 1. getting  pos control cells ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# key idea: for each response gene name, in the order 
# they appear in the data, get the cells expressing any gRNA
# targeting that gene. Keep doing this until we have at least
# `num_pc_cells` many cells.
all_response_gene_names <- rownames(response_odm) 
target_to_grna_ids <- split(grna_target_df$grna_id, grna_target_df$grna_target)
num_cells_per_grna <- list()
cells_per_grna <- list()
genes_kept <- c()
for(gene in all_response_gene_names) {
  # a. get the corresponding guide rnas
  curr_grna_ids <- target_to_grna_ids[[gene]]
  if(is.null(curr_grna_ids)) {
    next  # skip this gene if it's not a target
  } 
  genes_kept <- c(genes_kept, gene)
  # b. get the cells that express each of these gRNAs
  
  for(grna in curr_grna_ids) {
    cells_per_grna[[grna]] <- which(grna_assign_mat_to_use[grna, ])
    num_cells_per_grna[[grna]] <- length(cells_per_grna[[grna]])
  }

  if(sum(unlist(num_cells_per_grna)) >= num_pc_cells) {
    break
  }
}

## 2. adding in the NT cells
nt_grna_ids <- target_to_grna_ids[[NT_name]]
num_nt_so_far <- 0
for(grna in nt_grna_ids) {
  cells_per_grna[[grna]] <- which(grna_assign_mat_to_use[grna, ])
  num_cells_per_grna[[grna]] <- length(cells_per_grna[[grna]])
  num_nt_so_far <- num_nt_so_far + num_cells_per_grna[[grna]] 
  if(num_nt_so_far >= num_nt_cells) {
    break
  }
}



## 3. making cell info dataframe 
all_cell_idx <- unlist(cells_per_grna) |> setNames(NULL)
cell_info <- data.frame(stringsAsFactors = FALSE,
  cell_idx = all_cell_idx,
  grna_id = rep(names(cells_per_grna), num_cells_per_grna)
) |>
  dplyr::left_join(grna_target_df, by = "grna_id")
cell_info <- cbind(
  cell_info,
  cell_covariates[all_cell_idx, ] |>
    setNames(paste0(names(cell_covariates), "_full"))
)
  
## 4. getting the final data and writing
response_mat <- lapply(genes_kept,
       function(gene) response_odm[gene, ][all_cell_idx]) |>
  do.call(what = rbind)


write_fp <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/input_data", dataset_name
)
dir.create(write_fp, showWarnings = FALSE, recursive = TRUE)
write.csv(cell_info, file.path(write_fp, "cell_info.csv"), row.names = FALSE)

## 4a. sceptre
write_sceptre_fp <- file.path(write_fp, "sceptre")
dir.create(write_sceptre_fp, showWarnings = FALSE, recursive = TRUE)

grna_target_df_kept <- grna_target_df |>
  dplyr::filter(grna_id %in% names(cells_per_grna))
grna_matrix <- lapply(names(cells_per_grna), function(grna) {
  x <- grna_odm[grna, ][all_cell_idx] 
  j <- which(x != 0)
  sparseMatrix(
    i = rep.int(1L, length(j)),  # row index (all 1s)
    j = j,                       # column indices
    x = x[j],                    # nonzero values
    dims = c(1L, length(x))
  )
}) |>
  do.call(what = rbind)

saveRDS(response_mat, file.path(write_sceptre_fp, "response_matrix.rds"))
saveRDS(grna_matrix, file.path(write_sceptre_fp, "grna_matrix.rds"))
write.csv(grna_target_df_kept, file.path(write_sceptre_fp, "grna_target_data_frame.csv"), row.names = FALSE)


## 4b. Mixscale
write_mixscale_fp <- file.path(write_fp, "mixscale")
dir.create(write_mixscale_fp, showWarnings = FALSE, recursive = TRUE)
cell_names <- paste0("CELL_", cell_info$cell_idx)

saveRDS(
  response_mat |> `colnames<-`(cell_names),
  file.path(write_mixscale_fp, "response_matrix.rds")
)
saveRDS(
  setNames(cell_info$grna_target, cell_names),
  file.path(write_mixscale_fp, "assignments.rds")
)


# what do we need?
# - the response matrix
# - cell_info
# we then make the vector as follows: ~~~~~~~~~~~~~~~~~~~~~~~~
# cell_names <- paste0("CELL_", cell_info$cell_idx)
# colnames(response_mat) <- cell_names
# setNames(cell_info$grna_target, cell_names)
