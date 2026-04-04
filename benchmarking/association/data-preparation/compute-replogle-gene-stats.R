# let's determine the 100 most expressed genes
# even better would be the 100 most expressed for just the cells of interest

#ok or i can load the cell info, do this for just those cells, 

source("~/.Rprofile")

path_to_odm <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "guide_assignment/input_data/replogle-rd7/sceptre-pipeline"
)

response_odm = ondisc::initialize_odm_from_backing_file(file.path(path_to_odm, "response.odm"))


# getting the seleted cells
dataset_name = "replogle-rd7_comp_ngenes=100_ntargets=100_ncells=50k_bigger_genes"
path_to_subset <- file.path(
  .get_config_path("LOCAL_BENCHMARKING_DIR"),
  "association/computational/input_data",
  dataset_name
)
cell_info <- read_csv(file.path(path_to_subset, "cell_info.csv"))
cell_info$cell_idx |> unique() |> length()

cell_idx = cell_info$cell_idx

# getting sums
# out = data.frame(
#   gene_total = numeric(nrow(response_odm)),
#   gene_n_nonzero = 0,
#   gene_at_cell_idx_total = 0,
#   gene_at_cell_idx_n_nonzero = 0
# )
# for(i in 1:nrow(response_odm)) {
#   curr = response_odm[i,]
#   out$gene_total[i] = sum(curr)
#   out$gene_n_nonzero[i] = sum(curr > 0)
#   out$gene_at_cell_idx_total[i] = sum(curr[cell_idx])
#   out$gene_at_cell_idx_n_nonzero[i] = sum(curr[cell_idx] > 0)
# }
# out <- out |>
#   mutate(
#     gene = rownames(response_odm)
#   )
# 
# write_csv(out, file  ="~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline/gene_summary_stats.csv")

out = read_csv("~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/replogle-rd7/sceptre-pipeline/gene_summary_stats.csv")


dim(out)
head(out)
out |>
  mutate(
    gene = rownames(response_odm)
  ) |> 
  pivot_longer(!gene, names_to = "stat", values_to = "umi") |>
  mutate(
    log1p_umi = log1p(umi)
  ) |>
  ggplot(aes(x = log1p_umi, fill = stat)) +
  geom_histogram(bins=70) +
  facet_wrap(~stat, scale="free") +
  scale_y_log10() +
  theme_bw()

GGally::ggcorr(out |> log1p())
GGally::ggpairs(out |> sample_n(10000) |> log1p())
# pretty tight correlation, so i think we're taking enough cells that we can just work with the overall
# top umi ones vs needing to use the specific subset under consideration. That makes things easier.

head(out)


out |>
  filter(gene_at_cell_idx_n_nonzero == 0) |>
  dplyr::select(-gene_at_cell_idx_total, -gene_at_cell_idx_n_nonzero) |>
  pivot_longer(!gene) |>
  ggplot(aes(x = value, fill = name)) +
  geom_histogram(bins=50) +
  facet_wrap(~name) +
  scale_y_log10() +
  theme_bw()


out |>
  arrange(-gene_at_cell_idx_n_nonzero) |>
  head(100) |>
  pull(gene)


# Load data
cat("Loading Replogle RD7 data...\n")
scep <- sceptre::read_ondisc_backed_sceptre_object(
  sceptre_object_fp = file.path(path_to_odm, "sceptre_object.rds"),
  response_odm_file_fp = file.path(path_to_odm, "response.odm"),
  grna_odm_file_fp = file.path(path_to_odm, "grna.odm")
)


