library(sceptre)
library(sceptredata)
library(ggplot2)

# "grna_GACCTCC", "grna_TCCATAG"
data(highmoi_example_data)
data(grna_target_data_frame_highmoi)
response_matrix <- highmoi_example_data$response_matrix
grna_matrix <- highmoi_example_data$grna_matrix
rownames(grna_matrix)[rownames(grna_matrix) == "grna_GACCTCC"] <- "gRNA 1"
rownames(grna_matrix)[rownames(grna_matrix) == "grna_TCCATAG"] <- "gRNA 2"
grna_target_data_frame_highmoi$grna_id[grna_target_data_frame_highmoi$grna_id == "grna_GACCTCC"] <- "gRNA 1"
grna_target_data_frame_highmoi$grna_id[grna_target_data_frame_highmoi$grna_id == "grna_TCCATAG"] <- "gRNA 2"

sceptre_object <- import_data(
  response_matrix = highmoi_example_data$response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_data_frame_highmoi,
  moi = "high",
  extra_covariates = highmoi_example_data$extra_covariates,
  response_names = highmoi_example_data$gene_names
)

# 2. set the analysis parameters
positive_control_pairs <- construct_positive_control_pairs(sceptre_object)
discovery_pairs <- construct_cis_pairs(sceptre_object,
                                       positive_control_pairs = positive_control_pairs,
                                       distance_threshold = 5e6
)

sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = "left"
)
print(sceptre_object)

# 3. assign grnas
plot_grna_count_distributions(sceptre_object)
sceptre_object <- sceptre_object |> assign_grnas(parallel = TRUE, n_processors = 2)
p_1 <- plot_assign_grnas(sceptre_object,
                         return_indiv_plots = TRUE,
                         grnas_to_plot = c("gRNA 1", "gRNA 2"),
                         point_size = 0.3)[[1]] + theme(axis.title.y = element_blank())
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/aasign_grnas.png",
       plot = p_1, device = "png", scale = 0.6, width = 4, height = 3, dpi = 330)

# 4. run qc
plot_covariates(sceptre_object, p_mito_threshold = 0.075)
sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)
plot(sceptre_object)
print(sceptre_object)

# 5. run the calibration check
sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)
plot(sceptre_object)
print(sceptre_object)

# 6. run the power check
sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)
plot(sceptre_object)
print(sceptre_object)

# 7. run discovery analysis
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)
plot(sceptre_object)
print(sceptre_object)