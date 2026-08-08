# =============================================================================
# association-helpers.R
# Shared setup + helper functions for the association-testing benchmarking
# chapter. Source this once at the top of association-for-thesis.Rmd; every
# figure chunk then stands alone.
#
# Adapted from ~/katsevich-lab/code/sceptre3-project-v2/benchmarking/association/
# writeups/assoc-for-defense.Rmd (and its association-analysis-functions.R). The
# statistical loaders are ported here so this file is self-contained; only the
# QQ helpers still reach into sceptre::: internals (stat_qq_points/band).
# =============================================================================

library(Matrix)
library(tidyverse)
library(scales)
library(pROC)

source("~/.Rprofile")   # defines .get_config_path()

# -----------------------------------------------------------------------------
# Output location + one place to control how every figure is written
# -----------------------------------------------------------------------------
# Figures are written to a thesis-figures/ folder next to this writeup (created if
# absent). Relative to the working directory, i.e. this writeup's directory when the
# .Rmd is knit.
FIG_DIR <- "thesis-figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Default is vector PDF (sharpest in LaTeX). Switch device/ext here to change all
# figures at once (e.g. device = "png", ext = "png", plus dpi = 300).
save_fig <- function(plot, name, width = 6.5, height = 4.2,
                     ext = "pdf", device = grDevices::cairo_pdf, ...) {
  path <- file.path(FIG_DIR, paste0(name, ".", ext))
  ggsave(path, plot = plot, width = width, height = height, device = device, ...)
  invisible(path)
}

# -----------------------------------------------------------------------------
# Palettes + pretty-name maps
# -----------------------------------------------------------------------------
# Default ggplot2 discrete hues, matching the defense deck: SCEPTRE red, Mixscale
# teal, FR-Perturb purple. "pipeline" (compute plots only) is bright red.
assoc_method_pal <- c(
  SCEPTRE       = "#F8766D",
  Mixscale      = "#00BFC4",
  `FR-Perturb`  = "#C77CFF",
  pipeline      = "#FF0000"
)

# THE dataset palette, shared across both dissertation chapters (kept identical in
# grna-assignment/assignment-helpers.R): periwinkle for Replogle, amber for
# Gasperini. Used wherever a figure colours by dataset.
dataset_colour <- c(Replogle = "#7480CC", Gasperini = "#E3A245")

# Map the raw p-value-column suffixes to display names. Gasperini has no Mixscale.
assoc_method_renamer <- c(scep = "SCEPTRE", mix = "Mixscale", frpert = "FR-Perturb")

# Compute-trace process names -> display names.
assoc_method_renamer_comp <- c(
  sceptre     = "SCEPTRE",
  mixscale    = "Mixscale",
  frperturb   = "FR-Perturb",
  `scep-pipe` = "pipeline"
)

# Rename every present `p_value_<key>` column to `p_value_<DisplayName>`.
rename_pval_cols <- function(df) {
  dplyr::rename_with(
    df,
    .fn = function(nm) paste0("p_value_", assoc_method_renamer[sub("^p_value_", "", nm)]),
    .cols = dplyr::starts_with("p_value_")
  )
}

# -----------------------------------------------------------------------------
# Per-dataset run/dataset names (the two association benchmarks)
# -----------------------------------------------------------------------------
# Both datasets now live in one shared neg/pos run dir (run_all_max_t5); the
# dataset-name suffix in each filename ("replogle-rd7_max_*" / "gasperini_t5_*")
# is what distinguishes them, so the loaders are unchanged -- only these names are.
assoc_config <- list(
  replogle = list(
    label            = "Replogle",
    run_name_neg     = "run_all_max_t5",
    dataset_name_neg = "replogle-rd7_max_neg_ctrl_ngenes=1000_gene_thresh=7",
    run_name_pos     = "run_all_max_t5",
    dataset_name_pos = "replogle-rd7_max_pos_ctrl_ntargets=750_gene_thresh=7"
  ),
  gasperini = list(
    label            = "Gasperini",
    run_name_neg     = "run_all_max_t5",
    dataset_name_neg = "gasperini_t5_neg_ctrl_ngenes=1000_gene_thresh=7",
    run_name_pos     = "run_all_max_t5",
    dataset_name_pos = "gasperini_t5_pos_ctrl_ntargets=377_ncells=100k_gene_thresh=7"
  )
)

# -----------------------------------------------------------------------------
# Print theme (shared with the assignment chapter)
# -----------------------------------------------------------------------------
theme_thesis <- function(base_size = 10) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = rel(1.05), face = "bold",
                                       hjust = 0, margin = margin(b = 6)),
      plot.subtitle    = element_text(size = rel(0.9), colour = "grey35",
                                       hjust = 0, margin = margin(b = 6)),
      axis.title       = element_text(size = rel(1.0)),
      axis.text        = element_text(size = rel(0.85), colour = "grey25"),
      legend.position  = "bottom",
      legend.title     = element_text(size = rel(0.9), face = "bold"),
      legend.text      = element_text(size = rel(0.85)),
      legend.key.size  = unit(0.9, "lines"),
      legend.margin    = margin(t = 0),
      strip.text       = element_text(size = rel(0.9), face = "bold",
                                       margin = margin(4, 4, 4, 4)),
      strip.background = element_rect(fill = "grey92", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.border     = element_rect(colour = "grey70", fill = NA, linewidth = 0.4),
      plot.margin      = margin(5, 7, 5, 5)
    )
}
theme_set(theme_thesis())

# -----------------------------------------------------------------------------
# EDA: feature sparsity (genes + gRNAs)
# -----------------------------------------------------------------------------
# Both response.odm (genes) and grna.odm (gRNAs) ship with a per-feature
# <feature>_summary_stats.csv giving n_nonzero cells per feature, so no matrix
# scan is needed: the fraction of zero-count cells is 1 - n_nonzero / n_cells,
# with n_cells = the odm's column count (an instant metadata read).
gene_odm_dirs <- c(Replogle = "replogle-rd7", Gasperini = "gasperini")

# modality -> (summary csv, its id + n_nonzero columns, backing odm for cell count)
sparsity_modalities <- list(
  Gene = list(csv = "gene_summary_stats.csv", id = "gene", nz = "gene_n_nonzero", odm = "response.odm"),
  gRNA = list(csv = "grna_summary_stats.csv", id = "grna", nz = "grna_nonzero",   odm = "grna.odm")
)

# Per-feature stats for one dataset x modality (uses every feature): the feature
# id, the raw n_nonzero count, and the zero fraction 1 - n_nonzero / n_cells.
feature_stats <- function(odm_dir, modality) {
  m     <- sparsity_modalities[[modality]]
  dir   <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                     "guide_assignment", "input_data", odm_dir, "sceptre-pipeline")
  ncell <- ncol(ondisc::initialize_odm_from_backing_file(file.path(dir, m$odm)))
  ss    <- readr::read_csv(file.path(dir, m$csv), show_col_types = FALSE)
  tibble::tibble(feature = ss[[m$id]], n_nonzero = ss[[m$nz]], zero_frac = 1 - ss[[m$nz]] / ncell)
}

# Long tibble (dataset, modality, n_nonzero, zero_frac) over every dataset x
# modality. Feeds both plot_feature_sparsity() and plot_feature_nnzero().
load_feature_sparsity_df <- function(odm_dirs = gene_odm_dirs,
                                     modalities = names(sparsity_modalities)) {
  purrr::imap_dfr(odm_dirs, function(dir, ds)
    purrr::map_dfr(modalities, function(mod)
      dplyr::mutate(feature_stats(dir, mod), dataset = ds, modality = mod, .before = 1)))
}

# Coerce dataset/modality to their canonical factor orders (rows / columns).
.feature_factor <- function(df)
  dplyr::mutate(df,
    dataset  = factor(dataset,  levels = names(gene_odm_dirs)),
    modality = factor(modality, levels = names(sparsity_modalities)))

# -----------------------------------------------------------------------------
# Gene-wise QC (the genes actually carried into the association analysis)
# -----------------------------------------------------------------------------
# Which guide-assignment feeds the gene-QC MOI calculation, per dataset.
gene_qc_thresh       <- 7L
qc_assignment_method <- c(Replogle = "maximum", Gasperini = "thresholding-5")

# gene ids passing compute_genes_passing_qc() -- copied verbatim (behaviour-wise)
# from association/data-preparation/lib/io.R: a threshold on the expected number
# of perturbed cells (for a typical target) in which the gene is nonzero. NT stays
# in the guide/target counts; only nt_off_target and unknown are dropped. Because
# gene_n_nonzero is the only per-gene term, this is algebraically a floor on it.
load_gene_qc_pass <- function(odm_dir, assignment_method, thresh = gene_qc_thresh,
                              ignore_targets = c("nt_off_target", "unknown")) {
  base    <- .get_config_path("LOCAL_BENCHMARKING_DIR")
  pipe_fp <- file.path(base, "guide_assignment", "input_data", odm_dir, "sceptre-pipeline")
  gss  <- readr::read_csv(file.path(pipe_fp, "gene_summary_stats.csv"), show_col_types = FALSE)
  assn <- readRDS(file.path(base, "guide_assignment", "outputs", odm_dir,
                            assignment_method, "grna_assignment_matrix.rds"))
  gtd  <- readRDS(file.path(pipe_fp, "sceptre_object.rds"))@grna_target_data_frame

  moi_observed <- mean(Matrix::colSums(assn))
  keep         <- !gtd$grna_target %in% ignore_targets
  total_num_guides <- length(unique(gtd$grna_id[keep]))
  avg_num_guides_per_target <- gtd[keep, , drop = FALSE] |>
    dplyr::group_by(grna_target) |> dplyr::summarize(n = dplyr::n(), .groups = "drop") |>
    dplyr::pull(n) |> mean()
  qc_stat <- gss$gene_n_nonzero * moi_observed / total_num_guides * avg_num_guides_per_target
  gss$gene[qc_stat >= thresh]
}

# load_feature_sparsity_df() plus, for each dataset, the QC-passing gene rows
# re-labelled modality = "Gene (QC pass)" (a subset of the "Gene" rows). Slow
# (reads assignment matrices + sceptre objects), so cache the calling chunk.
load_feature_sparsity_qc_df <- function(odm_dirs = gene_odm_dirs,
                                        assignment_methods = qc_assignment_method) {
  base_df <- load_feature_sparsity_df(odm_dirs)
  qc_rows <- purrr::imap_dfr(odm_dirs, function(dir, ds) {
    qc_genes <- load_gene_qc_pass(dir, assignment_methods[[ds]])
    base_df |>
      dplyr::filter(dataset == ds, modality == "Gene", feature %in% qc_genes) |>
      dplyr::mutate(modality = "Gene (QC pass)")
  })
  dplyr::bind_rows(base_df, qc_rows)
}

# 2x2 histogram grid: rows = dataset, cols = modality. y is the fraction of
# features per bin -- populations differ by orders of magnitude (tens of
# thousands of genes vs a few thousand gRNAs), so raw counts aren't comparable
# across panels; normalizing each panel to sum to 1 puts them on one axis.
plot_feature_sparsity <- function(df, title = "Feature sparsity: fraction of zero-count cells",
                                  bins = 30, fill = "#4C72B0") {
  ggplot(.feature_factor(df), aes(zero_frac, y = after_stat(count / sum(count)))) +
    geom_histogram(bins = bins, boundary = 0, fill = fill,
                   colour = "white", linewidth = 0.2) +
    facet_grid(dataset ~ modality) +
    scale_x_continuous(labels = scales::label_percent(), expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(labels = scales::label_percent()) +
    labs(title = title,
         x = "Fraction of cells with zero count (per feature)",
         y = "Fraction of features") +
    theme_thesis()
}

# 2x2 histogram grid of the raw number of nonzero cells per feature. Scales are
# free per panel: the four dataset x modality populations span very different
# count ranges (and the two datasets have different cell totals), so no shared
# axis is meaningful. Panel order is still rows = dataset, cols = modality.
plot_feature_nnzero <- function(df, title = "Number of nonzero cells per feature",
                                bins = 30, fill = "#4C72B0") {
  ggplot(.feature_factor(df), aes(n_nonzero)) +
    geom_histogram(bins = bins, boundary = 0, fill = fill,
                   colour = "white", linewidth = 0.2) +
    facet_wrap(vars(dataset, modality), scales = "free", ncol = 2,
               labeller = labeller(.multi_line = FALSE)) +
    scale_x_continuous(labels = scales::label_number(big.mark = ","),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = title,
         x = "Number of nonzero cells (per feature)",
         y = "Number of features") +
    theme_thesis()
}

# Overlaid density of the detection rate (fraction of cells in which a feature is
# nonzero), gene vs gRNA, one panel per dataset. This is the view that shows all
# three facts at once: genes span many orders of magnitude (broad, skewed) while
# gRNAs are tightly clustered; the gene tail reaches far to the right (densest
# gene >> densest guide); and the dashed medians show which modality is typically
# denser (gRNA in Replogle, gene in Gasperini). Density is computed in log10 space
# (x is log-transformed) so the KDE isn't distorted; features detected in zero
# cells have no log and are dropped (reported via the caller).
feature_density_pal <- c(Gene = "#4C72B0", gRNA = "#DD8452")
# Same view with the QC-passing gene subset added as a third (green) curve.
gene_qc_density_pal <- c(Gene = "#4C72B0", `Gene (QC pass)` = "#55A868", gRNA = "#DD8452")

# Shared core for the detection-rate overlays. Groups are the `modality` levels of
# `pal` (2 or 3 curves); density is computed in log10 space so the KDE isn't
# distorted. Dashed medians (labelled with the percent, over ALL features so
# dropping never-detected zeros doesn't inflate them) show the typical feature;
# never-detected features can't be logged and are noted per panel. Median labels
# are staggered vertically per panel so nearby lines don't collide.
.density_overlay <- function(df, pal, title, y_lab = "Density (per group)") {
  d_all <- df |> dplyr::mutate(
    dataset  = factor(dataset,  levels = names(gene_odm_dirs)),
    modality = factor(modality, levels = names(pal)),
    density  = 1 - zero_frac)
  d <- d_all |> dplyr::filter(density > 0) |> dplyr::mutate(log_density = log10(density))
  meds <- d_all |> dplyr::group_by(dataset, modality) |>
    dplyr::summarise(med = median(density), .groups = "drop") |>
    dplyr::filter(med > 0) |>
    dplyr::mutate(log_density = log10(med), lab = paste0(signif(med * 100, 2), "%")) |>
    dplyr::arrange(dataset, log_density) |>
    dplyr::group_by(dataset) |> dplyr::mutate(vj = 1.4 + (dplyr::row_number() - 1) * 1.4) |>
    dplyr::ungroup()
  dropped <- d_all |> dplyr::group_by(dataset, modality) |>
    dplyr::summarise(n_zero = sum(density <= 0), .groups = "drop") |>
    dplyr::filter(n_zero > 0) |>
    dplyr::group_by(dataset) |> dplyr::mutate(row = dplyr::row_number()) |> dplyr::ungroup() |>
    dplyr::mutate(lab = sprintf("%s %s never detected (excluded)",
                                scales::comma(n_zero),
                                ifelse(grepl("RNA", modality), "gRNAs", "genes")))
  pct_lab <- function(b) scales::label_percent(accuracy = NULL)(10^b)  # 10^k -> % label
  ggplot(d, aes(log_density, colour = modality, fill = modality)) +
    geom_density(alpha = 0.22, linewidth = 0.6) +
    geom_vline(data = meds, aes(xintercept = log_density, colour = modality),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    geom_text(data = meds, aes(x = log_density, label = lab, colour = modality, vjust = vj),
              y = Inf, hjust = 0, nudge_x = 0.1, size = 2.9, fontface = "bold",
              inherit.aes = FALSE, show.legend = FALSE) +
    geom_text(data = dropped, aes(label = lab, colour = modality, vjust = 1.6 + (row - 1) * 1.5),
              x = -Inf, y = Inf, hjust = -0.03, size = 2.6,
              inherit.aes = FALSE, show.legend = FALSE) +
    facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = -6:0, labels = pct_lab) +
    scale_colour_manual(values = pal, name = NULL, aesthetics = c("colour", "fill")) +
    labs(title = title,
         x = "Fraction of cells in which the feature is detected (log scale)",
         y = y_lab) +
    theme_thesis()
}

# Two-curve version: gene vs gRNA (the original §0 figure).
plot_feature_density_overlay <- function(df, title = "Feature detection rate: genes vs gRNAs",
                                         pal = feature_density_pal)
  .density_overlay(df, pal, title)

# Three-curve version: all genes, QC-passing genes, and gRNAs. Feed it the df from
# load_feature_sparsity_qc_df() (which appends the "Gene (QC pass)" rows).
plot_gene_qc_density_overlay <- function(df,
    title = "Feature detection rate: all genes, QC-passing genes, gRNAs",
    pal = gene_qc_density_pal)
  .density_overlay(df, pal, title)

# =============================================================================
# Loaders: negative- and positive-control association results -> one tibble.
# (Ported from association-analysis-functions.R.)
# =============================================================================
load_associations_neg <- function(base_fp, run_name, dataset_name) {
  curr_fp <- file.path(base_fp, run_name)
  scep <- read.csv(paste0(curr_fp, "/association_neg_control_sceptre_", dataset_name, ".csv"))

  mixscale_fp <- paste0(curr_fp, "/association_neg_control_mixscale_", dataset_name, ".csv")
  mix <- if (file.exists(mixscale_fp)) read.csv(mixscale_fp) |> dplyr::select(-gene_ID) else NULL

  # FR-Perturb: a gene x target matrix of p-values -> melt to long.
  frpert_pvals_raw <- read.table(paste0(curr_fp, "/frperturb/frperturb_results_", dataset_name, "_pvals.txt"))
  frpert <- frpert_pvals_raw |>
    mutate(gene = rownames(frpert_pvals_raw)) |>
    pivot_longer(cols = !gene, values_to = "p_value_frpert", names_to = "target") |>
    mutate(target = str_replace_all(target, "\\.", "-"))

  out <- list(scep = scep, mix = mix, frpert = frpert)
  out[!sapply(out, is.null)]
}

combine_and_prepare_results_neg <- function(results_list_neg, disc_pairs, eps = NULL) {
  results_neg <- disc_pairs |>
    transmute(target = grna_target, gene = response_id) |>
    left_join(
      results_list_neg$scep |>
        dplyr::select(target = grna_target, gene = response_id,
                      n_nonzero_trt, n_nonzero_cntrl, p_value_scep = p_value),
      by = c("target", "gene")
    ) |>
    left_join(
      results_list_neg$frpert |> dplyr::select(target, gene, p_value_frpert),
      by = c("target", "gene")
    ) |>
    mutate(pair_type = "negative control")

  if (is.numeric(eps)) {
    results_neg <- results_neg |>
      mutate(p_value_scep   = pmax(pmin(p_value_scep, 1), eps),
             p_value_frpert = pmax(pmin(p_value_frpert, 1), eps))
  }
  if ("mix" %in% names(results_list_neg)) {
    results_neg <- results_neg |>
      left_join(
        results_list_neg$mix |>
          dplyr::select(target = grna_target, gene = response_id, p_value_mix = p_weight),
        by = c("target", "gene")
      )
    if (is.numeric(eps)) {
      results_neg <- results_neg |> mutate(p_value_mix = pmax(pmin(p_value_mix, 1), eps))
    }
  }
  results_neg
}

load_associations_pos <- function(base_fp, run_name, dataset_name) {
  base_fp <- file.path(base_fp, run_name)
  scep <- read_csv(paste0(base_fp, "/association_on_target_sceptre_", dataset_name, ".csv"),
                   show_col_types = FALSE) |>
    dplyr::select(target = grna_target, gene = response_id,
                  n_nonzero_trt, n_nonzero_cntrl, p_value_scep = p_value)

  mixscale_fp <- paste0(base_fp, "/association_on_target_mixscale_", dataset_name, ".csv")
  mix <- if (file.exists(mixscale_fp)) {
    suppressMessages(read_csv(mixscale_fp, show_col_types = FALSE)) |>
      dplyr::select(target = ...1, gene = gene_ID, p_value_mix = p_weight)
  } else NULL

  frpert_pvals_raw <- read.table(paste0(base_fp, "/frperturb/frperturb_results_", dataset_name, "_pvals.txt"))
  frpert <- tibble(
    target = rownames(frpert_pvals_raw),
    gene   = target,
    p_value_frpert = sapply(rownames(frpert_pvals_raw), function(t) frpert_pvals_raw[t, t])
  )

  out <- list(scep = scep, mix = mix, frpert = frpert)
  out[!sapply(out, is.null)]
}

combine_and_prepare_results_pos <- function(results_list_pos, eps = NULL) {
  results_pos <- results_list_pos$scep |>
    left_join(results_list_pos$frpert, by = c("target", "gene"))
  if ("mix" %in% names(results_list_pos)) {
    results_pos <- results_pos |> left_join(results_list_pos$mix, by = c("target", "gene"))
  }
  if (is.numeric(eps)) {
    results_pos <- results_pos |>
      mutate(p_value_scep   = pmax(pmin(p_value_scep, 1), eps),
             p_value_frpert = pmax(pmin(p_value_frpert, 1), eps))
    if ("mix" %in% names(results_list_pos)) {
      results_pos <- results_pos |> mutate(p_value_mix = pmax(pmin(p_value_mix, 1), eps))
    }
  }
  results_pos
}

# One-call bundle: load neg + pos controls for a dataset and stack them into a
# single results tibble (pair_type in {negative control, positive control}) with
# p-value columns renamed to display names (p_value_SCEPTRE, ...).
load_assoc_bundle <- function(which = c("replogle", "gasperini"), eps = 1e-250) {
  which <- match.arg(which)
  cfg   <- assoc_config[[which]]
  base  <- .get_config_path("LOCAL_BENCHMARKING_DIR")
  base_neg <- file.path(base, "association/neg-control")
  base_pos <- file.path(base, "association/pos-control")

  disc_pairs <- read_rds(file.path(base_neg, "input_data", cfg$dataset_name_neg,
                                   "sceptre/discovery_pairs.rds")) |> as_tibble()

  results_list_neg <- load_associations_neg(file.path(base_neg, "outputs"),
                                             cfg$run_name_neg, cfg$dataset_name_neg)
  results_neg <- combine_and_prepare_results_neg(results_list_neg, disc_pairs, eps = eps)

  results_list_pos <- load_associations_pos(file.path(base_pos, "outputs"),
                                             cfg$run_name_pos, cfg$dataset_name_pos)
  results_pos <- combine_and_prepare_results_pos(results_list_pos, eps = eps)

  results <- rbind(
    results_neg,
    mutate(results_pos, pair_type = "positive control")[names(results_neg)]
  )
  rename_pval_cols(results)
}

# =============================================================================
# §1 EDA — marginal p-value distributions (violins) + pairwise comparison
# =============================================================================
# One figure per dataset: the marginal distribution of -log10 p-values, a violin
# per method, faceted by pair type (positive control | negative control). A quick
# read of each method's overall p-value behaviour before the detailed analysis.
#
#   p_floor  cap p-values (original scale) at this value before -log10, so the
#            -log10 axis tops out at -log10(p_floor). Note load_assoc_bundle()
#            already floors at 1e-250; raise p_floor here to cap tighter and see
#            what reads best (e.g. 1e-50 -> axis caps at 50).
#   log1p_y  put the -log10 p axis on a log1p scale with breaks at 0 and powers of
#            10 (evenly spaced, not the scrunched-up naive log1p breaks). Appends
#            "(log10 scale)" to the y-axis label.
plot_pval_violin <- function(results, dataset_label, p_floor = 1e-250, log1p_y = FALSE) {
  df <- results |>
    dplyr::select(pair_type, dplyr::starts_with("p_value")) |>
    pivot_longer(dplyr::starts_with("p_value"), names_to = "method", values_to = "p") |>
    mutate(
      method    = sub("^p_value_", "", method),
      neglog10p = -log10(pmax(p, p_floor)),
      method    = factor(method, levels = intersect(names(assoc_method_pal), unique(method))),
      pair_type = factor(str_to_sentence(pair_type),
                         levels = c("Positive control", "Negative control"))
    ) |>
    filter(!is.na(neglog10p), is.finite(neglog10p))

  y_lab <- "p-value (-log10 scale)"
  y_scale <- NULL
  if (log1p_y) {
    max_v <- max(df$neglog10p, na.rm = TRUE)
    pw    <- if (max_v >= 1) 10^(0:floor(log10(max_v))) else numeric(0)
    y_scale <- scale_y_continuous(trans = "log1p", breaks = c(0, pw),
                                  labels = scales::label_number(accuracy = 1))
    y_lab <- paste0(y_lab, " (log10 scale)")
  }

  ggplot(df, aes(method, neglog10p, fill = method)) +
    geom_violin(trim = TRUE, scale = "width", colour = "grey30", linewidth = 0.3) +
    facet_wrap(~ pair_type, scales = "free_y") +
    scale_fill_manual(values = assoc_method_pal, guide = "none") +
    y_scale +
    labs(title = paste0(dataset_label, ": marginal p-value distributions"),
         x = "Method", y = y_lab) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# =============================================================================
# §2 EDA — pairwise p-value comparison against SCEPTRE
# =============================================================================
# Faceted scatter of each other method's -log10 p-value vs SCEPTRE's, split by
# pair type. Negative controls can be downsampled so the panel isn't a solid blob.
plot_pval_comparison <- function(results, dataset_label, reference_col = "p_value_SCEPTRE",
                                 downsample_n = 10000, max_log_p = 50,
                                 point_size = 0.8, transparency = 0.5, seed = 42) {
  comparison_cols <- setdiff(grep("^p_value_", names(results), value = TRUE), reference_col)

  df <- results
  if (!is.null(downsample_n)) {
    df_nc <- df |> filter(pair_type == "negative control")
    df_ot <- df |> filter(pair_type != "negative control")
    if (nrow(df_nc) > downsample_n) {
      set.seed(seed)
      df_nc <- df_nc |> sample_n(downsample_n)
    }
    df <- bind_rows(df_nc, df_ot)
  }

  df_long <- df |>
    pivot_longer(all_of(comparison_cols), names_to = "comparison_method",
                 values_to = "comparison_pvalue") |>
    mutate(
      x_plot = pmin(-log10(.data[[reference_col]]), max_log_p),
      y_plot = pmin(-log10(comparison_pvalue), max_log_p),
      comparison_method = str_remove(comparison_method, "^p_value_"),
      pair_type = str_to_sentence(pair_type)
    )

  ggplot(df_long, aes(x_plot, y_plot)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    geom_point(alpha = transparency, size = point_size, colour = "dodgerblue4") +
    facet_wrap(comparison_method ~ pair_type, scales = "free") +
    labs(
      title = paste0(dataset_label, ": comparison of p-values by pair type"),
      x = paste0("-log10(", str_remove(reference_col, "^p_value_"), ")"),
      y = "-log10(comparison p-value)"
    ) +
    theme_thesis()
}

# Hexbin version of plot_pval_comparison: same faceted -log10(p) comparison against
# SCEPTRE, but density is shown as hexbins (viridis, log10 count) instead of points.
# No downsampling is needed -- hexbins encode density directly, so every pair is
# used. The dashed y = x line is drawn on top for reference.
plot_pval_comparison_hex <- function(results, dataset_label, reference_col = "p_value_SCEPTRE",
                                     max_log_p = 50, bins = 50,
                                     title = paste0(dataset_label, ": comparison of p-values by pair type"),
                                     legend = c("right", "inside"), square = TRUE,
                                     neg_clip = NULL) {
  legend <- match.arg(legend)
  comparison_cols <- setdiff(grep("^p_value_", names(results), value = TRUE), reference_col)

  df_long <- results |>
    pivot_longer(all_of(comparison_cols), names_to = "comparison_method",
                 values_to = "comparison_pvalue") |>
    mutate(
      x_plot = pmin(-log10(.data[[reference_col]]), max_log_p),
      y_plot = pmin(-log10(comparison_pvalue), max_log_p),
      comparison_method = str_remove(comparison_method, "^p_value_"),
      pair_type = str_to_sentence(pair_type)
    )

  # Optional per-method cap on -log10(p) for NEGATIVE CONTROLS only: neg_clip is a
  # named vector method -> max -log10 p (equivalently, a p-value floor). Methods not
  # named are left at max_log_p. Applied independently to the reference method (x =
  # SCEPTRE) and to each comparison method (y). Positive controls are untouched.
  if (!is.null(neg_clip)) {
    ref_name <- str_remove(reference_col, "^p_value_")
    is_neg   <- df_long$pair_type == "Negative control"
    y_cap    <- unname(neg_clip[df_long$comparison_method])   # NA where method unspecified
    df_long$y_plot <- ifelse(is_neg & !is.na(y_cap), pmin(df_long$y_plot, y_cap), df_long$y_plot)
    if (ref_name %in% names(neg_clip))
      df_long$x_plot <- ifelse(is_neg, pmin(df_long$x_plot, neg_clip[[ref_name]]), df_long$x_plot)
  }

  # "right": thin vertical colourbar on the side (never overlaps data). "inside":
  # horizontal colourbar in the empty top-RIGHT corner of the top-left facet
  # (FR-Perturb neg control), to free panel width for tight layouts (e.g. the poster).
  if (legend == "inside") {
    fill_guide   <- guide_colourbar(barwidth = unit(3, "lines"), barheight = unit(0.4, "lines"))
    legend_theme <- theme(legend.position = "inside",
                          legend.position.inside = c(0.46, 0.99),
                          legend.justification.inside = c(1, 1),
                          legend.direction = "horizontal",
                          legend.background = element_rect(fill = alpha("white", 0.65), colour = NA),
                          legend.title = element_text(size = rel(0.75)),
                          legend.text  = element_text(size = rel(0.65)))
  } else {
    fill_guide   <- guide_colourbar(barwidth = unit(0.4, "lines"), barheight = unit(5, "lines"))
    legend_theme <- theme(legend.position = "right")
  }

  ggplot(df_long, aes(x_plot, y_plot)) +
    geom_hex(bins = bins) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    scale_fill_viridis_c(trans = "log10", name = "Count",
                         labels = scales::label_number(big.mark = ","),
                         guide = fill_guide) +
    facet_wrap(comparison_method ~ pair_type, scales = "free") +
    labs(
      title = title,
      x = paste0("-log10(", str_remove(reference_col, "^p_value_"), ")"),
      y = "-log10(comparison p-value)"
    ) +
    theme_thesis() +
    # square = TRUE keeps panels square (good for a 1x2 facet set on its own, and for
    # lining the combined 3x2 up); FALSE lets the facets fill the space (so the plot
    # doesn't shrink to a square with whitespace inside a wider/taller cell).
    (if (square) theme(aspect.ratio = 1) else NULL) +
    legend_theme
}

# Both datasets' pairwise-p-value hexbins in one 3x2 figure: Replogle (2 comparison
# methods x 2 pair types = 2x2) tagged (a) on top, Gasperini (1 method x 2 = 1x2)
# tagged (b) below. heights = c(2, 1) plus the square panels line the six up into a
# 3x2 grid. Each block keeps its own count colourbar (per-dataset count ranges).
plot_pval_comparison_hex_combined <- function(results_repl, results_gasp,
                                              title = "Pairwise p-values (vs SCEPTRE)",
                                              neg_clip = NULL) {
  p_r <- plot_pval_comparison_hex(results_repl, "Replogle", title = "Replogle", neg_clip = neg_clip)
  p_g <- plot_pval_comparison_hex(results_gasp, "Gasperini", title = "Gasperini", neg_clip = neg_clip)
  patchwork::wrap_plots(list(p_r, p_g), ncol = 1, heights = c(2, 1)) +
    patchwork::plot_annotation(
      title = title, tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      theme = theme_thesis() + theme(plot.title = element_text(size = 14))) &
    theme(plot.tag = element_text(face = "bold", size = 13))
}

# Histogram of each method's p-values for the pairs where n_nonzero_trt == 0 (no
# perturbed cell expresses the gene, so there is no signal), faceted method (rows)
# x pair type (cols). A well-behaved method should look ~uniform on negative
# controls even here; a pile at 0 or 1 flags degenerate behaviour. y is free per
# panel (negative controls vastly outnumber positives).
plot_pval_hist_zero_trt <- function(results, dataset_label, bins = 30) {
  df <- results |>
    dplyr::filter(n_nonzero_trt == 0) |>
    pivot_longer(dplyr::starts_with("p_value_"), names_to = "method", values_to = "p_value") |>
    dplyr::filter(!is.na(p_value)) |>
    mutate(method    = sub("^p_value_", "", method),
           method    = factor(method, levels = intersect(names(assoc_method_pal), unique(method))),
           pair_type = str_to_sentence(pair_type))
  ggplot(df, aes(p_value, fill = method)) +
    geom_histogram(breaks = seq(0, 1, length.out = bins + 1),
                   colour = "white", linewidth = 0.1) +
    facet_grid(method ~ pair_type, scales = "free_y") +
    scale_fill_manual(values = assoc_method_pal, guide = "none") +
    labs(title = paste0(dataset_label, ": p-values at n_nonzero_trt = 0"),
         x = "p-value", y = "Count") +
    theme_thesis()
}

# =============================================================================
# §2 Calibration — QQ plot (tail) of negative-control p-values by method
# =============================================================================
plot_qq_tail <- function(results, dataset_label, point_size = 1, transparency = 0.9) {
  df <- results |>
    filter(pair_type == "negative control") |>
    dplyr::select(starts_with("p_value")) |>
    pivot_longer(everything(), names_to = "method", values_to = "pvalue") |>
    mutate(method = sub("^p_value_", "", method)) |>
    filter(!is.na(pvalue))

  ggplot(df, aes(y = pvalue, colour = method, group = method)) +
    sceptre:::stat_qq_points(ymin = 1e-8, size = point_size, alpha = transparency) +
    sceptre:::stat_qq_band() +
    scale_x_continuous(trans = sceptre:::revlog_trans(10)) +
    scale_y_continuous(trans = sceptre:::revlog_trans(10)) +
    geom_abline(col = "black", linetype = "dashed") +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    labs(title = paste0(dataset_label, ": miscalibration in negative controls"),
         x = "Expected null p-value", y = "Observed p-value") +
    theme_thesis()
}

# =============================================================================
# §3 Power — ROC separating positive from negative controls
# =============================================================================
plot_assoc_roc <- function(results, dataset_label) {
  if (!setequal(results$pair_type, c("positive control", "negative control"))) {
    stop("'pair_type' contains unexpected values.")
  }
  is_pos   <- results$pair_type == "positive control"
  p_cols   <- grep("^p_value_", names(results), value = TRUE)
  methods  <- sub("^p_value_", "", p_cols)

  roc_list <- lapply(p_cols, function(pc) {
    pROC::roc(response = is_pos, predictor = -log10(results[[pc]]),
              levels = c(FALSE, TRUE), direction = "<")
  })
  aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))

  roc_data <- purrr::pmap_dfr(list(roc_list, methods), function(r, nm) {
    tibble(method = nm, fpr = 1 - r$specificities, tpr = r$sensitivities)
  })
  lev <- intersect(names(assoc_method_pal), methods)   # palette order, present only
  roc_data$method <- factor(roc_data$method, levels = lev)
  auc_labels <- setNames(
    sprintf("%s (AUC = %.3f)", lev, aucs[match(lev, methods)]), lev
  )

  ggplot(roc_data, aes(fpr, tpr, colour = method)) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = assoc_method_pal, labels = auc_labels, name = NULL) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = paste0(dataset_label, ": separating positive and negative controls"),
         x = "False positive rate", y = "Recall") +
    coord_fixed() +
    theme_thesis() +
    theme(legend.position = c(0.62, 0.22),
          legend.background = element_rect(fill = "white", colour = "grey80"))
}

# ROC as above, but two curves per method: solid = all positive controls, dashed =
# after dropping positive controls with n_nonzero_trt == 0 (no perturbed cell
# expresses the gene, so they are undetectable by any method). The gap shows how
# much those "impossible" positives depress each method's ROC. Negative controls
# are unchanged. Method legend labels carry both AUCs (all -> filtered).
plot_assoc_roc_zero_trt <- function(results, dataset_label, title = "ROC without sparse pairs",
                                    zoom = NULL) {
  if (!setequal(results$pair_type, c("positive control", "negative control")))
    stop("'pair_type' contains unexpected values.")
  p_cols  <- grep("^p_value_", names(results), value = TRUE)
  methods <- sub("^p_value_", "", p_cols)
  lev     <- intersect(names(assoc_method_pal), methods)

  # solid = all pairs; dashed = only pairs with n_trt > 0 (the n_trt==0 positives,
  # undetectable by any method, are dropped).
  variants <- c(`all pairs` = "all", `n_trt > 0` = "drop")
  df_of <- function(v) if (v == "all") results else
    dplyr::filter(results, !(pair_type == "positive control" & n_nonzero_trt == 0))

  roc_data <- purrr::imap_dfr(variants, function(v, vlab) {
    df <- df_of(v)
    is_pos <- df$pair_type == "positive control"
    purrr::pmap_dfr(list(p_cols, methods), function(pc, nm) {
      r <- pROC::roc(response = is_pos, predictor = -log10(df[[pc]]),
                     levels = c(FALSE, TRUE), direction = "<")
      tibble(method = nm, variant = vlab, auc = as.numeric(pROC::auc(r)),
             fpr = 1 - r$specificities, tpr = r$sensitivities)
    })
  })
  roc_data$method  <- factor(roc_data$method,  levels = lev)
  roc_data$variant <- factor(roc_data$variant, levels = names(variants))

  aucs <- roc_data |> dplyr::distinct(method, variant, auc) |>
    tidyr::pivot_wider(names_from = variant, values_from = auc)
  method_labels <- setNames(
    sprintf("%s (%.3f -> %.3f)", aucs$method,
            aucs[["all pairs"]], aucs[["n_trt > 0"]]),
    as.character(aucs$method))[lev]

  # zoom = z (0<z<1) crops to the top-left corner [0,z] x [1-z,1] (equal spans keep
  # coord_fixed square) to separate near-perfect, overlapping curves.
  win <- if (is.null(zoom)) list(x = c(0, 1), y = c(0, 1)) else list(x = c(0, zoom), y = c(1 - zoom, 1))

  ggplot(roc_data, aes(fpr, tpr, colour = method, linetype = variant)) +
    geom_line(linewidth = 0.9, alpha = 0.7) +   # alpha: near-perfect curves overlap in the corner
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = assoc_method_pal, labels = method_labels, name = NULL) +
    scale_linetype_manual(values = c(`all pairs` = "solid", `n_trt > 0` = "22"), name = NULL) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = title, x = "False positive rate", y = "Recall") +
    coord_fixed(ratio = 1, xlim = win$x, ylim = win$y) +
    theme_thesis() +
    # long AUC labels -> stack each legend in a single column so it doesn't run wide
    guides(colour = guide_legend(ncol = 1, order = 1),
           linetype = guide_legend(ncol = 1, order = 2))
}

# One method's zero-trt p-value histogram with positive + negative controls STACKED
# in a single panel, coloured by pair type. Helper for plot_assoc_zero_trt_grid.
.hist_zero_trt_stacked <- function(results, method, pal, bins = 30, legend = FALSE) {
  col <- paste0("p_value_", method)
  df  <- results |>
    dplyr::filter(n_nonzero_trt == 0, !is.na(.data[[col]])) |>
    dplyr::mutate(pair_type = str_to_sentence(pair_type), p_value = .data[[col]])
  ggplot(df, aes(p_value, fill = pair_type)) +
    geom_histogram(breaks = seq(0, 1, length.out = bins + 1),
                   position = "stack", colour = "white", linewidth = 0.1) +
    scale_fill_manual(values = pal, name = NULL) +
    labs(title = method, x = "p-value", y = "Count") +
    theme_thesis() +
    theme(legend.position = if (legend) "bottom" else "none")
}

# 2x2 for one dataset (Replogle by default): three per-method stacked zero-trt
# p-value histograms (pos + neg control stacked, coloured by pair type) plus the
# two-curve ROC (plot_assoc_roc_zero_trt), assembled with cowplot. The histograms
# share ONE legend (extracted once) placed below the grid; the ROC keeps its own.
# Pair type is coloured with dataset_colour for now (neg = periwinkle, pos = amber).
plot_assoc_zero_trt_grid <- function(results, dataset_label = "Replogle", bins = 30,
                                     title = NULL, roc_zoom = NULL) {
  methods  <- intersect(names(assoc_method_pal),
                        sub("^p_value_", "", grep("^p_value_", names(results), value = TRUE)))
  pair_pal <- c(`Negative control` = unname(dataset_colour[["Replogle"]]),
                `Positive control` = unname(dataset_colour[["Gasperini"]]))

  hists <- lapply(methods, function(m) .hist_zero_trt_stacked(results, m, pair_pal, bins))
  roc   <- plot_assoc_roc_zero_trt(results, dataset_label, zoom = roc_zoom)

  hist_leg <- cowplot::get_legend(
    .hist_zero_trt_stacked(results, methods[[1]], pair_pal, bins, legend = TRUE) +
      theme(legend.direction = "horizontal"))

  grid <- cowplot::plot_grid(plotlist = c(hists, list(roc)), ncol = 2, labels = "auto")
  body <- cowplot::plot_grid(grid, hist_leg, ncol = 1, rel_heights = c(1, 0.05))
  if (is.null(title)) return(body)
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold", size = 16, x = 0.01, hjust = 0)
  cowplot::plot_grid(title_grob, body, ncol = 1, rel_heights = c(0.045, 1))
}

# =============================================================================
# §4 Power at a false-positive budget — #TP recovered vs #FP allowed
# =============================================================================
plot_fp_budget <- function(results, dataset_label,
                           fp_budgets = c(0, 1, 2, 5, 10, 20, 30, 50, 100, 250, 500),
                           y_var = c("tp", "tpr"),
                           pseudo_log_sigma_x = 1, pseudo_log_sigma_y = 1) {
  y_var     <- match.arg(y_var)
  pos_label <- "positive control"
  neg_label <- "negative control"

  tp_at_fp <- function(df_m) {
    n_pos <- sum(df_m$pair_type == pos_label)
    n_neg <- sum(df_m$pair_type == neg_label)
    steps <- df_m |>
      mutate(is_pos = pair_type == pos_label, is_neg = pair_type == neg_label) |>
      group_by(p_value) |>
      summarise(tp_step = sum(is_pos), fp_step = sum(is_neg), .groups = "drop") |>
      arrange(p_value) |>
      mutate(tp_cum = cumsum(tp_step), fp_cum = cumsum(fp_step))

    idx <- findInterval(fp_budgets, steps$fp_cum)
    fp_actual <- integer(length(fp_budgets)); tp <- integer(length(fp_budgets))
    threshold <- rep(NA_real_, length(fp_budgets))
    ok <- idx > 0
    fp_actual[ok] <- steps$fp_cum[idx[ok]]
    tp[ok]        <- steps$tp_cum[idx[ok]]
    threshold[ok] <- steps$p_value[idx[ok]]
    big <- fp_budgets >= n_neg
    fp_actual[big] <- n_neg; tp[big] <- n_pos; threshold[big] <- 1
    tibble(fp_budget = fp_budgets, fp_actual = fp_actual, threshold = threshold,
           tp = tp, tpr = if (n_pos > 0) tp / n_pos else NA_real_,
           n_pos = n_pos, n_neg = n_neg)
  }

  p_val_cols     <- grep("^p_value_", names(results), value = TRUE)
  p_cols_present <- setNames(sub("^p_value_", "", p_val_cols), p_val_cols)

  long <- results |>
    dplyr::select(pair_type, all_of(names(p_cols_present))) |>
    pivot_longer(all_of(names(p_cols_present)), names_to = "method_col", values_to = "p_value") |>
    mutate(method = recode(method_col, !!!as.list(p_cols_present)),
           p_value = as.numeric(p_value)) |>
    filter(pair_type %in% c(pos_label, neg_label), !is.na(p_value)) |>
    dplyr::select(method, pair_type, p_value)

  plot_df <- long |>
    group_by(method) |>
    group_modify(~ tp_at_fp(.x)) |>
    ungroup() |>
    arrange(method, fp_actual, fp_budget) |>
    group_by(method, fp_actual) |>
    summarise(tp = max(tp), tpr = max(tpr), .groups = "drop") |>
    mutate(method = factor(method, levels = intersect(names(assoc_method_pal), method)))

  max_fp   <- max(plot_df$fp_actual, na.rm = TRUE)
  x_breaks <- c(0, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000,
                2000, 5000, 10000, 20000, 50000, 100000)
  x_breaks <- x_breaks[x_breaks <= max_fp]
  if (!0 %in% x_breaks) x_breaks <- c(0, x_breaks)

  if (y_var == "tp") {
    max_y    <- max(plot_df$tp, na.rm = TRUE)
    y_breaks <- c(0, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000, 1500, 2000, 5000, 10000)
    y_breaks <- y_breaks[y_breaks <= max_y]
    if (!0 %in% y_breaks) y_breaks <- c(0, y_breaks)
    y_lab   <- "#TP (positive controls recovered)"
    y_trans <- scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma_y)
    y_fmt   <- scales::label_number(big.mark = "")
  } else {
    y_breaks <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
    y_lab    <- "TPR (recall on positives)"
    y_trans  <- scales::pseudo_log_trans(base = 10, sigma = 0.01)
    y_fmt    <- scales::label_percent(accuracy = 1)
  }

  ggplot(plot_df, aes(x = fp_actual, y = .data[[y_var]], colour = method)) +
    geom_step(direction = "hv", linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma_x),
                       breaks = x_breaks, labels = scales::label_number(big.mark = "")) +
    scale_y_continuous(trans = y_trans, breaks = y_breaks, labels = y_fmt) +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    labs(title = paste0(dataset_label, ": #TP recovered vs #FP allowed"),
         x = "#FP (negative controls called significant)", y = y_lab) +
    theme_thesis()
}

# =============================================================================
# Combined per-dataset "poster": the four analyses in one tagged 2x2 figure
# =============================================================================
# For ONE dataset, lay out (a) plot_pval_comparison_hex, (b) plot_qq_tail,
# (c) plot_assoc_roc, (d) plot_fp_budget in a 2x2, tagged (a)-(d). Panel (a) is
# itself multi-faceted, so render this LARGE (e.g. fig.width=13, fig.height=11).
# Each panel keeps its own legend (poster style); the dataset is the overall title.
plot_assoc_poster <- function(results, dataset_label, title = dataset_label,
                              square = FALSE, legend = "inside", neg_clip = NULL) {
  a <- plot_pval_comparison_hex(results, dataset_label,
                                title = "Pairwise p-values (vs SCEPTRE)",
                                legend = legend, square = square, neg_clip = neg_clip)
  b <- plot_qq_tail(results, dataset_label) +
    labs(title = "Negative-control calibration")
  c <- plot_assoc_roc(results, dataset_label) +
    labs(title = "Positive vs. negative separation")
  d <- plot_fp_budget(results, dataset_label) +
    labs(title = "Power at a false-positive budget")

  patchwork::wrap_plots(list(a, b, c, d), ncol = 2) +
    patchwork::plot_annotation(
      title = title, tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
      theme = theme_thesis() + theme(plot.title = element_text(size = 16))
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))
}

# =============================================================================
# §5 FDR control — false-discovery proportion across positive prevalence (BH)
# =============================================================================
# Repeatedly draw a mixture of N_total pos/neg-control p-values at prevalence pi,
# BH-adjust at q_target, and record the realized FDP. Median +/- IQR over B draws.
plot_bh_vs_pi <- function(results, dataset_label, N_total, pi_grid,
                          q_target = 0.1, B = 100, p_floor = 1e-250) {
  set.seed(1)
  pos_label <- "positive control"; neg_label <- "negative control"

  p_val_cols     <- grep("^p_value_", names(results), value = TRUE)
  p_cols_present <- setNames(sub("^p_value_", "", p_val_cols), p_val_cols)

  pools <- results |>
    filter(pair_type %in% c(pos_label, neg_label)) |>
    dplyr::select(pair_type, all_of(names(p_cols_present))) |>
    pivot_longer(all_of(names(p_cols_present)), names_to = "method_col", values_to = "p") |>
    filter(!is.na(p)) |>
    mutate(method = recode(method_col, !!!as.list(p_cols_present)),
           p = pmax(as.numeric(p), p_floor),
           is_pos = pair_type == pos_label, is_neg = pair_type == neg_label) |>
    dplyr::select(method, p, is_pos, is_neg)

  bench <- tidyr::crossing(method = unique(pools$method), pi = pi_grid, b = seq_len(B)) |>
    group_by(method, pi, b) |>
    group_modify(~{
      m <- .y$method; pi <- .y$pi
      pos_pool <- pools |> filter(method == m, is_pos) |> pull(p)
      neg_pool <- pools |> filter(method == m, is_neg) |> pull(p)
      N_pos <- floor(pi * N_total); N_neg <- N_total - N_pos
      if (N_pos > length(pos_pool)) {
        stop("Trying to sample ", N_pos, " pos. pairs from ", length(pos_pool),
             " total.\n pi = ", round(pi, 4))
      }
      pos_samp <- sample(pos_pool, size = N_pos, replace = FALSE)
      neg_samp <- sample(neg_pool, size = N_neg, replace = FALSE)
      p_mix <- c(pos_samp, neg_samp)
      lab   <- c(rep(TRUE, N_pos), rep(FALSE, N_neg))
      disc  <- p.adjust(p_mix, method = "BH") <= q_target
      tp <- sum(disc & lab); fp <- sum(disc & !lab); D <- tp + fp
      tibble(discoveries = D, tp = tp, fp = fp,
             fdp = if (D == 0) 0 else fp / D,
             tpr = if (N_pos == 0) NA_real_ else tp / N_pos)
    }) |>
    ungroup()

  bench_sum <- bench |>
    group_by(method, pi) |>
    summarise(fdp_med = mean(fdp, na.rm = TRUE),
              fdp_lo  = quantile(fdp, 0.25, na.rm = TRUE),
              fdp_hi  = quantile(fdp, 0.75, na.rm = TRUE), .groups = "drop") |>
    mutate(method = factor(method, levels = intersect(names(assoc_method_pal), unique(method))))

  ggplot(bench_sum, aes(pi, fdp_med, colour = method, fill = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = fdp_lo, ymax = fdp_hi), alpha = 0.15, colour = NA) +
    geom_hline(yintercept = q_target, linetype = "dashed") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    scale_fill_manual(values = assoc_method_pal, name = "Method") +
    labs(title = paste0(dataset_label, ": FDP for varying positive prevalence"),
         subtitle = paste0("BH q = ", q_target, ", N = ", N_total, ", B = ", B),
         x = "Positive prevalence", y = "FDP") +
    theme_thesis()
}

# =============================================================================
# On-target significance vs. effective sample size
# =============================================================================
# Positive controls only: each method's on-target -log10 p-value against SCEPTRE's
# n_nonzero_trt (the number of perturbed cells with nonzero expression — the
# effective treatment sample size, shared by all methods for a given target/gene).
# Shows whether a method's significance grows with sample size or plateaus.
plot_pval_vs_ntrt <- function(results, dataset_label, size_col = "n_nonzero_trt",
                              smooth_method = "loess") {
  df <- results |>
    filter(pair_type == "positive control") |>
    dplyr::select(dplyr::all_of(size_col), dplyr::starts_with("p_value")) |>
    pivot_longer(dplyr::starts_with("p_value"), names_to = "method", values_to = "p") |>
    mutate(method    = sub("^p_value_", "", method),
           neglog10p = -log10(p),
           method    = factor(method, levels = intersect(names(assoc_method_pal), unique(method)))) |>
    filter(!is.na(neglog10p), is.finite(neglog10p), !is.na(.data[[size_col]]))

  ggplot(df, aes(.data[[size_col]], neglog10p, colour = method)) +
    geom_point(alpha = 0.35, size = 0.8) +
    geom_smooth(se = FALSE, method = smooth_method, linewidth = 0.9) +
    scale_x_log10() +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    labs(title = paste0(dataset_label, ": on-target significance vs. sample size"),
         x = "Effective treatment sample size (SCEPTRE n_nonzero_trt, log scale)",
         y = "p-value (-log10 scale)") +
    theme_thesis()
}

# =============================================================================
# §7 Computational benchmarking — runtime + memory vs number of guides
# =============================================================================
# Per-dataset compute config: the in-memory association run and a glob for the
# pipeline dataset dirs (10 each: ngenes {500,1000} x ntargets). Only ngenes==1000
# points are plotted (x = ntargets); the ngenes=500 dirs are dropped by the filter.
assoc_comp_config <- list(
  replogle = list(
    label           = "Replogle",
    run_names_inmem = "run_repl_max_upto_ng1000nt800",
    scep_pipe_glob  = "replogle-rd7_max_*"
  ),
  gasperini = list(
    label           = "Gasperini",
    run_names_inmem = "run_gasp_t5_upto_ng1000nt800",
    scep_pipe_glob  = "gasperini_t5_*"
  )
)

# In-memory association-method traces (SCEPTRE, Mixscale, FR-Perturb). Parses
# peak_rss -> GB and realtime (Hh Mm Ss) -> seconds.
load_trace_files_in_memory <- function(run_names, base_path_comp) {
  lapply(run_names, function(run_name) {
    read_tsv(file.path(base_path_comp, "outputs", run_name, "trace.txt"),
             show_col_types = FALSE)
  }) |>
    do.call(what = rbind) |>
    tidyr::extract(tag, into = c("ngenes", "ntargets"),
                   regex = "ngenes=([0-9]+)_ntargets=([0-9]+)",
                   remove = FALSE, convert = TRUE) |>
    transmute(method = gsub("_computational", "", process, ignore.case = TRUE) |> tolower(),
              peak_rss, realtime, ngenes, ntargets) |>
    mutate(peak_rss_gb = case_when(
      str_detect(peak_rss, regex("\\bGB\\b", ignore_case = TRUE)) ~ parse_number(peak_rss),
      str_detect(peak_rss, regex("\\bMB\\b", ignore_case = TRUE)) ~ parse_number(peak_rss) / 1024,
      TRUE ~ NA_real_
    )) |>
    mutate(
      h = coalesce(parse_number(str_extract(realtime, "\\d+(?=h)")), 0),
      m = coalesce(parse_number(str_extract(realtime, "\\d+(?=m)")), 0),
      s = coalesce(parse_number(str_extract(realtime, "\\d+(?=s)")), 0),
      realtime_sec = 3600 * h + 60 * m + s
    ) |>
    dplyr::select(-h, -m, -s)
}

# Distributed pipeline traces (run_association_analysis workers). Elapsed time =
# span from first submit to last (drift-corrected) complete; memory = max over
# workers. Time-queued is excluded by adding realtime to submit time.
load_trace_files_scep_pipe <- function(scep_pipe_datasets, scep_pipe_base_fp) {
  parse_mem_to_gb <- function(x) {
    x <- str_trim(x)
    value <- as.numeric(str_extract(x, "^[0-9.]+"))
    unit  <- str_extract(x, "[A-Za-z]+$") |> toupper()
    mult <- case_when(
      unit == "B" ~ 1 / 1024^3, unit == "KB" ~ 1 / 1024^2, unit == "MB" ~ 1 / 1024,
      unit == "GB" ~ 1, unit == "TB" ~ 1024, TRUE ~ NA_real_
    )
    ifelse(is.na(value * mult), 0, value * mult)
  }
  lapply(scep_pipe_datasets, function(dataset) {
    curr_trace <- file.path(scep_pipe_base_fp, dataset, "tracing/trace.tsv") |>
      read_tsv(show_col_types = FALSE) |>
      mutate(start    = as.POSIXct(start,    tz = "America/New_York"),
             complete = as.POSIXct(complete, tz = "America/New_York"),
             submit   = as.POSIXct(submit,   tz = "America/New_York"),
             peak_rss_in_GB = parse_mem_to_gb(peak_rss))
    assoc_data <- curr_trace |>
      filter(grepl("run_analysis_subworkflow_.*:run_association_analysis", process))
    start_time <- min(assoc_data$submit)
    end_time   <- max(assoc_data$submit + (assoc_data$complete - assoc_data$start))
    data.frame(
      dataset = dataset,
      assoc_elapsed_time_in_sec = as.numeric(end_time - start_time, units = "secs"),
      max_peak_rss_in_gb_for_assoc_workers = max(assoc_data$peak_rss_in_GB),
      num_assoc_analysis_workers = nrow(assoc_data)
    )
  }) |>
    do.call(what = rbind) |>
    tidyr::extract(dataset, into = c("ngenes", "ntargets"),
                   regex = "ngenes=([0-9]+)_ntargets=([0-9]+)",
                   remove = FALSE, convert = TRUE) |>
    mutate(method = "scep-pipe")
}

# Bundle: stack in-memory + pipeline traces into one frame
# (method, mem_in_gb, runtime_in_sec, ngenes, ntargets).
load_assoc_comp_bundle <- function(which = c("replogle", "gasperini")) {
  which <- match.arg(which)
  cfg   <- assoc_comp_config[[which]]
  base  <- .get_config_path("LOCAL_BENCHMARKING_DIR")
  base_path_comp    <- file.path(base, "association/computational")
  scep_pipe_base_fp <- file.path(base, "association/computational/outputs/sceptre-pipeline")

  # The pipeline dataset dirs are globbed (10 per dataset) rather than listed.
  scep_pipe_datasets <- basename(Sys.glob(file.path(scep_pipe_base_fp, cfg$scep_pipe_glob)))
  trace_inmem <- load_trace_files_in_memory(cfg$run_names_inmem, base_path_comp)
  trace_pipe  <- load_trace_files_scep_pipe(scep_pipe_datasets, scep_pipe_base_fp)

  rbind(
    trace_inmem |>
      dplyr::select(method, mem_in_gb = peak_rss_gb,
                    runtime_in_sec = realtime_sec, ngenes, ntargets),
    trace_pipe |>
      dplyr::select(method, mem_in_gb = max_peak_rss_in_gb_for_assoc_workers,
                    runtime_in_sec = assoc_elapsed_time_in_sec, ngenes, ntargets)
  )
}

# Filter to one gene count, prettify method, and tag each point's execution mode
# (drives circle vs triangle: in-memory methods vs the distributed pipeline).
.prep_assoc_comp <- function(trace, ngenes_keep = 1000) {
  trace |>
    filter(ngenes == ngenes_keep) |>
    mutate(
      method    = factor(assoc_method_renamer_comp[method], levels = assoc_method_renamer_comp),
      execution = factor(ifelse(method == "pipeline", "Distributed", "In-memory"),
                         levels = c("In-memory", "Distributed"))
    ) |>
    droplevels()
}

# Circle = in-memory method, triangle = distributed pipeline.
comp_shape_scale <- function() {
  scale_shape_manual(values = c(`In-memory` = 16, `Distributed` = 17), name = "Execution")
}

# Fig 6a/6c: runtime vs number of guides (log-log), one dataset per figure.
plot_assoc_runtime <- function(trace, title, ngenes_keep = 1000) {
  fmt_hms <- function(x) { x <- round(x); sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60) }
  all_breaks <- c(30, 60, 120, 300, 600, 1800, 3600, 7200, 14400, 28800)
  d <- .prep_assoc_comp(trace, ngenes_keep)
  current_max <- max(d$runtime_in_sec, na.rm = TRUE)
  ggplot(d, aes(ntargets, runtime_in_sec, colour = method, group = method)) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 3) +
    scale_y_log10(breaks = all_breaks[all_breaks <= current_max * 1.2], labels = fmt_hms) +
    scale_x_log10() +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    comp_shape_scale() +
    labs(x = "Number of targets (log scale)", y = "Runtime (H:M:S, log scale)", title = title) +
    theme_thesis()
}

# Fig 6b/6d: peak memory vs number of guides (log-log), one dataset per figure.
plot_assoc_memory <- function(trace, title, ngenes_keep = 1000) {
  d <- .prep_assoc_comp(trace, ngenes_keep)
  ggplot(d, aes(ntargets, mem_in_gb, colour = method, group = method)) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 3) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                  labels = scales::label_number(suffix = " GB")) +
    scale_x_log10() +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    comp_shape_scale() +
    labs(x = "Number of targets (log scale)", y = "Peak memory (log scale)", title = title) +
    theme_thesis()
}

# Combined faceted view of §8: dataset as facet columns (Replogle / Gasperini),
# runtime and memory stacked as rows. Mirrors the assignment chapter's
# plot_comp_grid_faceted -- two facet_grid(~ dataset) plots, each keeping its own
# y-scale so the exact H:M:S / GB ticks of plot_assoc_runtime/plot_assoc_memory
# survive, stacked with patchwork (aligned columns, one wrapped legend). The two
# datasets arrive as separate traces, so pass both. `title` sets a single overall
# title for the whole figure.
plot_assoc_comp_grid_faceted <- function(trace_repl, trace_gasp, title = NULL,
                                         ngenes_keep = 1000) {
  d <- dplyr::bind_rows(
    .prep_assoc_comp(trace_repl, ngenes_keep) |> dplyr::mutate(dataset = "Replogle"),
    .prep_assoc_comp(trace_gasp, ngenes_keep) |> dplyr::mutate(dataset = "Gasperini")
  ) |>
    dplyr::mutate(dataset = factor(dataset, levels = c("Replogle", "Gasperini")))

  fmt_hms <- function(x) { x <- round(x); sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60) }
  all_breaks <- c(30, 60, 120, 300, 600, 1800, 3600, 7200, 14400, 28800)
  rt_max <- max(d$runtime_in_sec, na.rm = TRUE)

  base_layers <- function(p) p +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 2.6) +
    facet_grid(~ dataset, scales = "free_x") +   # dataset columns; shared y within each metric
    scale_x_log10() +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    comp_shape_scale() +
    theme_thesis() +
    guides(colour = guide_legend(nrow = 2))   # wrap Method legend so it fits

  p_rt <- base_layers(ggplot(d, aes(ntargets, runtime_in_sec, colour = method, group = method))) +
    scale_y_log10(breaks = all_breaks[all_breaks <= rt_max * 1.2], labels = fmt_hms) +
    labs(x = NULL, y = "Runtime (H:M:S, log scale)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  # x shown on bottom plot only

  p_mm <- base_layers(ggplot(d, aes(ntargets, mem_in_gb, colour = method, group = method))) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                  labels = scales::label_number(suffix = " GB")) +
    labs(x = "Number of targets (log scale)", y = "Peak memory (log scale)") +
    theme(strip.text.x = element_blank(),           # dataset strips on the top plot only
          strip.background.x = element_blank())

  grid <- patchwork::wrap_plots(list(p_rt, p_mm), ncol = 1, guides = "collect")
  if (!is.null(title))
    grid <- grid + patchwork::plot_annotation(
      title = title, theme = theme_thesis() + theme(plot.title = element_text(size = 14)))
  grid & theme(legend.position = "bottom")
}

# Like .prep_assoc_comp but keeps ALL gene counts (no ngenes filter); adds ngenes
# as a factor so it can drive linetype.
.prep_assoc_comp_ng <- function(trace) {
  trace |>
    mutate(
      method    = factor(assoc_method_renamer_comp[method], levels = assoc_method_renamer_comp),
      execution = factor(ifelse(method == "pipeline", "Distributed", "In-memory"),
                         levels = c("In-memory", "Distributed")),
      ngenes_f  = factor(ngenes)
    ) |>
    droplevels()
}

# Like plot_assoc_comp_grid_faceted, but two curves per method -- one per gene count
# -- across all 10 datasets (ngenes {500, 1000} x ntargets). colour = method,
# linetype = ngenes (largest = solid), shape = execution. Same stacked runtime /
# memory layout with dataset facet columns and one shared legend.
plot_assoc_comp_grid_faceted_ngenes <- function(trace_repl, trace_gasp, title = NULL) {
  d <- dplyr::bind_rows(
    .prep_assoc_comp_ng(trace_repl) |> dplyr::mutate(dataset = "Replogle"),
    .prep_assoc_comp_ng(trace_gasp) |> dplyr::mutate(dataset = "Gasperini")
  ) |>
    dplyr::mutate(dataset = factor(dataset, levels = c("Replogle", "Gasperini")))

  fmt_hms <- function(x) { x <- round(x); sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60) }
  all_breaks <- c(30, 60, 120, 300, 600, 1800, 3600, 7200, 14400, 28800)
  rt_max <- max(d$runtime_in_sec, na.rm = TRUE)
  ng_lev <- levels(d$ngenes_f)
  ng_lty <- setNames(c("solid", "22", "42", "11")[rank(-as.numeric(ng_lev))], ng_lev)  # largest ngenes -> solid

  base_layers <- function(p) p +
    geom_line(aes(linetype = ngenes_f), linewidth = 0.8, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 2.6) +
    facet_grid(~ dataset, scales = "free_x") +
    scale_x_log10() +
    scale_colour_manual(values = assoc_method_pal, name = "Method") +
    scale_linetype_manual(values = ng_lty, name = "# genes") +
    comp_shape_scale() +
    theme_thesis() +
    guides(colour = guide_legend(nrow = 2))

  p_rt <- base_layers(ggplot(d, aes(ntargets, runtime_in_sec, colour = method,
                                    group = interaction(method, ngenes_f)))) +
    scale_y_log10(breaks = all_breaks[all_breaks <= rt_max * 1.2], labels = fmt_hms) +
    labs(x = NULL, y = "Runtime (H:M:S, log scale)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_mm <- base_layers(ggplot(d, aes(ntargets, mem_in_gb, colour = method,
                                    group = interaction(method, ngenes_f)))) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                  labels = scales::label_number(suffix = " GB")) +
    labs(x = "Number of targets (log scale)", y = "Peak memory (log scale)") +
    theme(strip.text.x = element_blank(), strip.background.x = element_blank())

  grid <- patchwork::wrap_plots(list(p_rt, p_mm), ncol = 1, guides = "collect")
  if (!is.null(title))
    grid <- grid + patchwork::plot_annotation(
      title = title, theme = theme_thesis() + theme(plot.title = element_text(size = 14)))
  grid & theme(legend.position = "bottom")
}
