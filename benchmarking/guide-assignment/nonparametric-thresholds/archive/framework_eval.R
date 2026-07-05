#!/usr/bin/env Rscript
# Insert the two nonparametric thresholding methods (Otsu on log1p, smoothed
# valley) into the formal precision/recall/Jaccard comparison used in
# assignment-analysis/assignment-writeup.Rmd (simulated-data section).
#
# The four baselines are scored with the framework's OWN code (metric +
# process_* functions copied verbatim from the Rmd) so their numbers match the
# canonical plots exactly; the two new methods are scored identically.

suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
HERE  <- dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE))))
BENCH <- path.expand("~/data/projects/sceptre3/benchmarking")
source(file.path(HERE, "..", "guide-assignment-pipeline", "bin", "script", "lib",
                 "threshold_methods.R"))

# ---- metric functions (verbatim from assignment-writeup.Rmd:118-149) --------
jaccard <- function(s1, s2) { i <- length(intersect(s1, s2)); u <- length(union(s1, s2))
  if (u == 0) return(NA_real_); i / u }
precision <- function(s1, s2) { i <- length(intersect(s1, s2)); d <- length(s1)
  if (d == 0) return(NA_real_); i / d }
recall <- function(s1, s2) { i <- length(intersect(s1, s2)); d <- length(s2)
  if (d == 0) return(NA_real_); i / d }
F1 <- function(s1, s2) { i <- length(intersect(s1, s2)); d <- length(s1) + length(s2)
  if (d == 0) return(NA_real_); 2 * i / d }
metrics <- list(Jaccard = jaccard, Precision = precision, Recall = recall, F1 = F1)

# ---- method loaders (verbatim from assignment-writeup.Rmd:413-453) ----------
process_sceptre <- function(fp, assn_list) {
  raw <- readRDS(fp); a <- assn_list
  for (g in names(a)) a[[g]] <- which(raw[g, ]); a
}
process_crispat <- function(fp, assn_list) {
  raw <- read.csv(fp); ul <- split(raw$cell_id, raw$grna_id); a <- assn_list
  for (g in names(ul)) a[[g]] <- ul[[g]]; a
}
process_pert <- function(fp, assn_list) {
  raw <- read.csv(fp); a <- assn_list
  for (g in names(a)) a[[g]] <- raw$cell_id[grepl(g, raw$grna_id)]; a
}
process_clns <- function(fp, assn_list, ground_truth) {
  raw <- read.csv(fp); ul <- split(raw$cell_id, raw$grna_id); a <- assn_list
  for (i in names(ul)) a[[rownames(ground_truth)[[as.integer(i)]]]] <- ul[[i]]; a
}

# ---- new methods: per-guide cell-NAME sets, same shape as the baselines -----
process_threshold <- function(grna_matrix, fn, assn_list, cell_names) {
  am <- assign_by_threshold(grna_matrix, fn)$assignment_matrix   # grnas x cells
  a  <- assn_list
  for (g in names(a)) {
    idx <- if (g %in% rownames(am)) which(am[g, ]) else integer(0)
    a[[g]] <- cell_names[idx]
  }
  a
}

# ---- datasets (verbatim from assignment-writeup.Rmd:483-490) ----------------
sim_names <- list(
  list(dataset_name = "replogle-rd7_simulated_100k_0.015-pert",
       run_name = "run_all_sims_100k",            short_name = "Sim: 1.5% Pert."),
  list(dataset_name = "replogle-rd7_simulated_100k_0.005-pert",
       run_name = "run_all_sims_100k_0.005-pert", short_name = "Sim: 0.5% Pert.")
)

all_results   <- setNames(vector("list", length(sim_names)), sapply(sim_names, `[[`, "short_name"))
ground_truths <- all_results

for (sn in sim_names) {
  in_fp  <- file.path(BENCH, "guide_assignment/input_data", sn$dataset_name)
  out_fp <- file.path(BENCH, "guide_assignment/outputs",    sn$run_name)
  gt     <- readRDS(file.path(in_fp, "true_perturbation_status.rds"))
  cell_names <- colnames(gt)

  tmpl <- setNames(vector("list", nrow(gt)), rownames(gt))
  truth <- tmpl
  for (g in names(truth)) truth[[g]] <- names(which(gt[g, ] == 1))
  ground_truths[[sn$short_name]] <- truth

  grna_matrix <- readRDS(file.path(in_fp, "sceptre", "grna_matrix.rds"))

  scep <- process_sceptre(file.path(out_fp, paste0("assignment_matrix_sceptre_", sn$dataset_name, ".rds")), tmpl)
  clns <- process_clns(file.path(out_fp, paste0("assignments_cleanser_", sn$dataset_name, ".csv")), tmpl, gt)
  for (g in names(scep)) { scep[[g]] <- cell_names[scep[[g]]]; clns[[g]] <- cell_names[clns[[g]]] }
  crisp <- process_crispat(file.path(out_fp, paste0("assignments_crispat_", sn$dataset_name, ".csv")), tmpl)
  pert  <- process_pert(file.path(out_fp, paste0("assignments_pertpy_", sn$dataset_name, ".csv")), tmpl)

  otsu  <- process_threshold(grna_matrix, otsu_threshold_log1p,     tmpl, cell_names)
  vall  <- process_threshold(grna_matrix, smoothed_valley_threshold, tmpl, cell_names)

  all_results[[sn$short_name]] <- list(scep = scep, crisp = crisp, pert = pert,
                                       clns = clns, otsu = otsu, vall = vall)
  cat("Loaded", sn$dataset_name, "-", nrow(gt), "guides x", ncol(gt), "cells\n")
}

# ---- build metric_df (mirrors assignment-writeup.Rmd:588-632) ---------------
method_renamer <- c(scep = "SCEPTRE", crisp = "crispat", pert = "pertpy",
                    clns = "CLEANSER", otsu = "Otsu (log1p)", vall = "Smoothed valley")
dataset_names <- names(all_results)
method_names  <- names(all_results[[1]])

rows <- list(); idx <- 1L
for (mn in names(metrics)) for (d in seq_along(all_results)) {
  truth_d <- ground_truths[[d]]
  for (m in seq_along(all_results[[d]])) {
    preds <- all_results[[d]][[m]]
    vals  <- vapply(seq_along(preds), function(k) metrics[[mn]](preds[[k]], truth_d[[k]]), numeric(1))
    rows[[idx]] <- data.frame(dataset = names(all_results)[d],
                              method = names(all_results[[d]])[m],
                              value = vals, metric = mn); idx <- idx + 1L
  }
}
metric_df <- do.call(rbind, rows)
metric_df$value   <- ifelse(is.na(metric_df$value), 0, metric_df$value)
metric_df$dataset <- factor(metric_df$dataset, levels = dataset_names)
metric_df$method  <- factor(method_renamer[metric_df$method], levels = method_renamer)
metric_df$metric  <- factor(metric_df$metric, levels = names(metrics))

# numeric summary
summ <- aggregate(value ~ method + metric + dataset, metric_df, mean)
summ_wide <- reshape(summ, idvar = c("dataset", "method"), timevar = "metric", direction = "wide")
names(summ_wide) <- sub("value.", "", names(summ_wide))
cat("\n===== mean per-guide metrics =====\n"); print(summ_wide, row.names = FALSE, digits = 3)
write.csv(summ_wide, file.path(HERE, "results", "framework_eval_summary.csv"), row.names = FALSE)

# ---- plot: Jaccard / Precision / Recall, 6 methods, faceted by dataset ------
pd <- metric_df[metric_df$metric %in% c("Jaccard", "Precision", "Recall"), ]
pd$metric <- factor(pd$metric, levels = c("Precision", "Recall", "Jaccard"))
is_new <- pd$method %in% c("Otsu (log1p)", "Smoothed valley")
p <- ggplot(pd, aes(method, value, fill = is_new)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.9, linewidth = 0.3) +
  stat_summary(fun = mean, geom = "point", size = 1.1, color = "black") +
  facet_grid(dataset ~ metric) +
  scale_fill_manual(values = c(`FALSE` = "grey80", `TRUE` = "#3a7bd5"), guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Guide-level metrics vs. true perturbation status",
       subtitle = "blue = new nonparametric thresholds; dot = mean", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(HERE, "results", "framework_eval_metrics.png"), p,
       width = 10, height = 6, dpi = 120)
cat("\nWrote results/framework_eval_metrics.png\n")
