#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# bugfix-sims-metrics.R
#
# Self-contained analysis of the `run_bugfix_sims` guide-assignment run:
# two methods (glmpois_poisbug, glmpois_poisfix) on three simulation datasets.
# For each dataset it computes per-guide precision / recall / jaccard against
# the simulation ground truth (true_pert_matrix.rds) and plots the average
# metric vs the perturbed mean (mu_pert), faceted by NP tail type and metric.
#
# Adapted from sceptre-nb/additive-model.Rmd, but with all helpers inlined so
# it depends only on tidyverse + Matrix (no project sourcing).
#
# Usage:  Rscript bugfix-sims-metrics.R
# Output: PDFs in ./plots/run_bugfix_sims/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
})

## --- paths ---------------------------------------------------------------
source("~/.Rprofile")
BENCHMARK_DIR <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
RUN_NAME      <- "run_bugfix_sims"
OUTPUTS_DIR   <- file.path(BENCHMARK_DIR, "outputs", RUN_NAME)
INPUT_DIR     <- file.path(BENCHMARK_DIR, "input_data")

script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)))),
  error = function(e) getwd()
)
if (length(script_dir) == 0) script_dir <- getwd()
PLOT_DIR <- file.path(script_dir, "plots", RUN_NAME)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

## --- the datasets in this run -------------------------------------------
# Each has its own P-label -> mu_pert mapping (copied from additive-model.Rmd).
DATASETS <- list(
  sims_sum_2np_3p = tibble(
    p = c("P1", "P2", "P3"),
    mu_pert = c(500, 750, 1000)
  ),
  sims_sum_1np_3p_disp = tibble(
    p = c("Plowdisp", "Pmeddisp", "Phighdisp"),
    mu_pert = c(101000, 1001000, 10001000)
  ),
  sims_sum_repeat_old = tibble(
    p = c("Psmall", "Pmed", "Plarge"),
    mu_pert = c(970 / 8, 970 / 4, 970)
  )
)

## --- metric helpers (inlined from nb-bench_v3-helpers.R) -----------------
set_metrics <- function(pred, truth) {
  pred <- unique(pred); truth <- unique(truth)
  n_inter <- length(intersect(pred, truth))
  n_pred  <- length(pred); n_truth <- length(truth)
  n_union <- length(union(pred, truth))
  tibble(
    precision = if (n_pred  == 0) NA_real_ else n_inter / n_pred,
    recall    = if (n_truth == 0) NA_real_ else n_inter / n_truth,
    jaccard   = if (n_union == 0) NA_real_ else n_inter / n_union,
    n_pred = n_pred, n_truth = n_truth, n_inter = n_inter
  )
}

extract_part <- function(x, k) {
  vapply(strsplit(x, "_", fixed = TRUE),
         function(z) if (length(z) >= k) z[[k]] else NA_character_,
         character(1))
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_or_na   <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)

## --- load one dataset's run ---------------------------------------------
# Reads every assignment_matrix_*_<dataset>.rds in the run dir, converts each
# guide's assigned cells to integer indices, and derives the ground-truth
# assignment (perturbed AND observed nonzero) from true_pert_matrix.rds.
load_run <- function(dataset_name) {
  input_data_path <- file.path(INPUT_DIR, dataset_name)
  grna_matrix <- readRDS(file.path(input_data_path, "sceptre", "grna_matrix.rds"))
  true_perts  <- readRDS(file.path(input_data_path, "true_pert_matrix.rds"))
  guide_names <- rownames(grna_matrix)

  fnames <- dir(OUTPUTS_DIR)
  assn_fnames <- fnames[grepl("^assignment_matrix", fnames) &
                          grepl(dataset_name, fnames, fixed = TRUE)]

  all_assns <- list()
  for (fname in assn_fnames) {
    method_name <- fname |>
      sub(pattern = "^assignment_matrix_script_", replacement = "") |>
      sub(pattern = paste0("_", dataset_name, ".rds"), replacement = "", fixed = TRUE)
    mat <- readRDS(file.path(OUTPUTS_DIR, fname))
    all_assns[[method_name]] <- lapply(guide_names,
      function(g) as.integer(which(mat[g, ]))) |> setNames(guide_names)
  }

  true_assns <- lapply(guide_names, function(g)
    as.integer(which(true_perts[g, ] == 1 & grna_matrix[g, ] > 0))) |>
    setNames(guide_names)

  list(dataset_name = dataset_name, grna_matrix = grna_matrix,
       all_assns = all_assns, true_assns = true_assns)
}

## --- per-guide metric table for one run ---------------------------------
make_guide_metric_tbl <- function(run) {
  imap_dfr(run$all_assns, function(method_assns, method_name) {
    map_dfr(rownames(run$grna_matrix), function(grna_id) {
      set_metrics(method_assns[[grna_id]], run$true_assns[[grna_id]]) |>
        mutate(method_name = method_name, grna_id = grna_id, .before = 1)
    })
  })
}

summarize_metric_tbl <- function(tbl, group_vars) {
  tbl |>
    pivot_longer(c(precision, recall, jaccard),
                 names_to = "metric", values_to = "value") |>
    group_by(across(all_of(c(group_vars, "metric")))) |>
    summarize(mean = mean_or_na(value), sd = sd_or_na(value),
              n_nonmissing = sum(!is.na(value)), .groups = "drop")
}

## --- plotting ------------------------------------------------------------
make_avg_metrics_plot <- function(df, x_breaks = NULL) {
  if (is.null(x_breaks)) x_breaks <- c(200, 400, 800, 1600)
  df$metric <- factor(df$metric, levels = c("precision", "recall", "jaccard"))
  ggplot(df, aes(x = mu_pert, y = mean, group = method_name,
                 color = method_name, shape = method_name)) +
    geom_line() +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.04, alpha = 0.7) +
    geom_point(size = 2) +
    facet_grid(np ~ metric) +
    scale_x_continuous(trans = "log2", breaks = x_breaks) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8)) +
    labs(x = expression(mu[pert]), y = "Mean", color = "Method", shape = "Method") +
    theme_bw() +
    theme(legend.position = "bottom")
}

## --- run each dataset ----------------------------------------------------
analyze_dataset <- function(dataset_name, p_to_mean) {
  message("Loading ", dataset_name, " ...")
  run <- load_run(dataset_name)
  message("  methods: ", paste(names(run$all_assns), collapse = ", "))

  guide_meta <- tibble(
    grna_id = rownames(run$grna_matrix),
    np = extract_part(grna_id, 2),
    p  = extract_part(grna_id, 3)
  ) |> left_join(p_to_mean, by = "p")

  plot_df <- make_guide_metric_tbl(run) |>
    left_join(guide_meta, by = "grna_id") |>
    summarize_metric_tbl(group_vars = c("method_name", "np", "p", "mu_pert")) |>
    mutate(ci = 2 * sd / sqrt(n_nonmissing),
           ymin = pmax(0, mean - ci), ymax = pmin(1, mean + ci),
           method_name = sub("^script_", "", method_name))

  x_breaks <- sort(unique(p_to_mean$mu_pert))
  p <- make_avg_metrics_plot(plot_df, x_breaks = x_breaks) +
    ggtitle(paste0("Precision / recall / jaccard vs ground truth: ", dataset_name))

  out_fp <- file.path(PLOT_DIR, paste0("metrics_", dataset_name, ".pdf"))
  ggsave(out_fp, p, width = 8, height = 6)
  message("  wrote ", out_fp)
  invisible(plot_df)
}

for (ds in names(DATASETS)) analyze_dataset(ds, DATASETS[[ds]])
message("Done. Plots in ", PLOT_DIR)
