library(sceptre)
library(sceptredata)
library(katlabutils)
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
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/assign_grnas.png",
       plot = p_1, device = "png", scale = 0.6, width = 4, height = 3.5, dpi = 330)

# 4. run qc
sceptre_object <- sceptre_object |> run_qc(p_mito_threshold = 0.075)
ps <- plot_run_qc(sceptre_object, return_indiv_plots = TRUE)
p_3 <- ps[[1]] + theme(plot.title = element_blank())
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/run_qc.png",
       plot = p_3, device = "png", scale = 0.6, width = 4, height = 3.5, dpi = 330)
 
# 5. run the calibration check
sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = 2)
ps <- plot_run_calibration_check(sceptre_object = sceptre_object,
                                 return_indiv_plots = TRUE)
p_4 <- ps[[2]] + theme(axis.title.y = element_blank(),
                       axis.title.x = element_blank(),
                       plot.title = element_blank())
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/run_calibration_check.png",
       plot = p_4, device = "png", scale = 0.6, width = 4, height = 3.5, dpi = 330)

# 6. run the power check
load_all("~/research_code/sceptre")
plot_run_power_check <- function(sceptre_object, point_size = 1, transparency = 0.8, clip_to = 1e-20) {
  if (!sceptre_object@functs_called["run_power_check"]) {
    stop("This `sceptre_object` has not yet had `run_power_check` called on it.")
  }
  my_theme <- get_my_theme()
  set.seed(3)
  my_cols <- c("mediumseagreen", "firebrick1")
  
  pos_ctrl_pvals <- sceptre_object@power_result$p_value |> stats::na.omit()
  neg_ctrl_pval_sub <- if (nrow(sceptre_object@calibration_result) >= 1) {
    downsample_result_data_frame(
      result_df = sceptre_object@calibration_result
    ) |> dplyr::pull(p_value)
  } else {
    numeric(0)
  }
  group_names <- c("Pos. Control", "Neg. Control")
  df <- data.frame(
    lab = rep(
      group_names,
      c(length(pos_ctrl_pvals), length(neg_ctrl_pval_sub))
    ) |>
      factor(levels = group_names),
    p_values = c(pos_ctrl_pvals, neg_ctrl_pval_sub) |>
      pmax(clip_to)
  )
  
  p <- ggplot2::ggplot(
    data = df,
    mapping = ggplot2::aes(x = lab, y = p_values, color = lab)
  )
  my_breaks <- 10^(seq(-2, -40, by = -4))
  p <- p +
    ggplot2::geom_jitter(width = .25, height = 0, size = point_size, alpha = transparency) +
    ggplot2::scale_y_continuous(trans = revlog_trans(base = 10), expand = c(0.01, 0), breaks = my_breaks) +
    ggplot2::scale_color_manual(values = my_cols, guide = "none") +
    ggplot2::labs(
      x = "Pair type",
      y = "p-value",
      title = "Positive and negative control p-values"
    ) +
    my_theme
  
  return(p)
}

sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = 2)
p_5 <- plot_run_power_check(sceptre_object = sceptre_object, clip_to = 0) +
  scale_y_continuous(breaks = 10^(-seq(10, 250, length.out = 5)),
                     trans = revlog_trans(base = 10)) +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank())
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/run_power_check.png",
       plot = p_5, device = "png", scale = 0.6, width = 4, height = 3.5, dpi = 330)

# 7. run discovery analysis
sceptre_object <- run_discovery_analysis(sceptre_object, parallel = TRUE, n_processors = 2)
ps <- plot_run_discovery_analysis(sceptre_object, return_indiv_plots = TRUE)
p_6 <- ps[[2]] + theme(axis.title.y = element_blank(),
                       axis.title.x = element_blank(),
                       plot.title = element_blank())
ggsave(filename = "~/research_code/sceptre3-project-v2/figs/run_discovery_analysis.png",
       plot = p_6, device = "png", scale = 0.6, width = 4, height = 3.5, dpi = 330)
