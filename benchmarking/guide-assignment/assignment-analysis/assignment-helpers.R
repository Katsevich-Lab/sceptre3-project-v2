# =============================================================================
# assignment-helpers.R
# Shared setup + helper functions for the gRNA-assignment benchmarking chapter.
# Source this once at the top of assignment-for-thesis.Rmd; every figure chunk
# then stands alone.
# =============================================================================

library(Matrix)
library(tidyverse)
library(data.table)

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
# Palettes
# -----------------------------------------------------------------------------
# Default ggplot2 discrete hues (hue_pal), matching the defense deck.
method_pal <- c(
  SCEPTRE  = "#F8766D",  # red
  crispat  = "#7CAE00",  # green
  CLEANSER = "#00BFC4",  # teal
  `pertpy (v1.0.5)` = "#C77CFF",  # purple
  pipeline = "#FF0000"   # bright red (compute plots only)
)
sparsity_pal <- c(low = "#56B4E9", high = "#E69F00")

# THE dataset palette, shared across both dissertation chapters (kept identical in
# association-helpers.R): periwinkle for Replogle, amber for Gasperini. Chosen to
# stay clear of the four method hues (red, green, teal, purple) and to read under
# red-green colour-vision deficiency (blue-vs-amber is on the preserved blue-yellow
# axis). Used for colour and (with alpha) for fill.
dataset_colour <- c(Replogle = "#7480CC", Gasperini = "#E3A245")

# Colors for the two per-item statistics used to pick example guides/cells.
# Default ggplot2 2-color hues (red + teal). (Adjust here to restyle every overlay.)
feature_pal <- c(`Nonzero UMIs` = "#F8766D", `Total UMIs` = "#00BFC4")

# -----------------------------------------------------------------------------
# Pretty-name maps
# -----------------------------------------------------------------------------
method_renamer <- c(
  sceptre  = "SCEPTRE",
  crispat  = "crispat",
  cleanser = "CLEANSER",
  pertpy   = "pertpy (v1.0.5)"
)
dataset_renamer <- c(
  ng100 = "100 guides",
  ng200 = "200 guides",
  ng400 = "400 guides",
  ng800 = "800 guides"
)
source_renamer <- c(replogle = "Replogle", gasperini = "Gasperini")

data_sources <- c("replogle-rd7", "gasperini")
methods      <- c("sceptre", "crispat", "cleanser", "pertpy")

# Ambient vs perturbed components in the simulation figures (teal / red).
sim_status_pal <- c(Ambient = "#00BFC4", `Pert.` = "#F8766D")

# Perturbation strength (mu_1 = 971/2^mupow), treated as a discrete factor. Four
# high-contrast colors chosen to NOT collide with the method palette (which owns
# the default ggplot red/green/teal/purple). Names are the mu values so
# scale_colour_manual matches the factor levels.
mu_pal <- c(`121` = "#E69F00", `243` = "#0072B2", `486` = "#CC79A7", `971` = "#000000")

# The two calling-boundary quantities in the Fig 2e boxplots (blue = the noise
# ceiling a method tolerates, orange = the signal floor it requires).
boundary_pal <- c(`Largest UMI left unperturbed`  = "#0072B2",
                  `Smallest UMI called perturbed` = "#E69F00")

# -----------------------------------------------------------------------------
# Print theme (replaces the presentation theme)
# -----------------------------------------------------------------------------
# Tuned for a printed page, not a projector: smaller base size, thin recessive
# grid, on-plot title optional (LaTeX \caption usually carries the description).
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

# =============================================================================
# Loaders: per-method gRNA assignments -> named list (guide -> vector of cells)
# =============================================================================
get_all_guides_per_dataset <- function(dataset_name) {
  curr_path_to_data <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "guide_assignment/input_data", dataset_name, "sceptre/grna_matrix.rds"
  )
  grna_matrix <- read_rds(curr_path_to_data)
  rownames(grna_matrix)
}

make_assn_data_structure <- function(dataset_names, methods) {
  assns_list <- vector("list", length(dataset_names)) |>
    setNames(names(dataset_names))
  for (i in seq_along(dataset_names)) {
    assns_list[[i]] <- vector("list", length(methods)) |>
      setNames(methods)
    curr_guides <- get_all_guides_per_dataset(dataset_name = dataset_names[i])
    guide_to_assns_template <- vector("list", length(curr_guides)) |>
      setNames(curr_guides)
    for (j in seq_along(methods)) {
      assns_list[[i]][[j]] <- guide_to_assns_template
    }
  }
  assns_list
}

get_run_with_dataset_and_method <- function(run_names, dataset_name, method_name, base_data_dir) {
  which_run_has_this_dataset_and_method <- sapply(run_names, function(run_name) {
    curr_run_dir <- file.path(base_data_dir, run_name)
    curr_run_files <- dir(curr_run_dir)
    file_found <- grepl(method_name, curr_run_files, fixed = TRUE) &
      grepl(dataset_name, curr_run_files, fixed = TRUE)
    if (sum(file_found) > 1) {
      stop("multiple files matched for '", dataset_name, "' and ", method_name, " in '", run_name, "'\n")
    }
    any(file_found)
  })
  if (sum(which_run_has_this_dataset_and_method) > 1) {
    stop("multiple runs matched for '", dataset_name, "' and ", method_name, "\n")
  } else if (sum(which_run_has_this_dataset_and_method) == 0) {
    stop("no runs matched for '", dataset_name, "' and ", method_name, "\n")
  }
  run_names[which_run_has_this_dataset_and_method]
}

load_all_assns <- function(run_names, dataset_names, methods) {
  base_data_dir <- file.path(
    .get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/outputs"
  )
  all_assns <- make_assn_data_structure(dataset_names = dataset_names, methods = methods)

  for (dataset_idx in seq_along(dataset_names)) {
    dataset_name <- dataset_names[dataset_idx]
    dataset_name_short <- names(dataset_names)[dataset_idx]

    path_to_cell_names <- file.path(
      .get_config_path("LOCAL_BENCHMARKING_DIR"),
      "guide_assignment/input_data", dataset_name, "sceptre/grna_matrix.rds"
    )
    cell_names_in_order <- read_rds(path_to_cell_names) |> colnames()

    for (method_name in methods) {
      run_with_this_dataset <- get_run_with_dataset_and_method(
        run_names = run_names, dataset_name = dataset_name,
        method_name = method_name, base_data_dir = base_data_dir
      )
      curr_run_dir <- file.path(base_data_dir, run_with_this_dataset)
      curr_files <- dir(curr_run_dir)
      curr_assns_fname <- curr_files[
        grepl(dataset_name, curr_files, fixed = TRUE) &
          grepl(method_name, curr_files, fixed = TRUE)
      ]
      curr_assns <- load_assns(
        path_to_assns = file.path(curr_run_dir, curr_assns_fname),
        method_name = method_name,
        all_guides = names(all_assns[[dataset_name_short]][[method_name]]),
        cell_names_in_order = cell_names_in_order
      )
      if (!all(names(curr_assns) == names(all_assns[[dataset_name_short]][[method_name]]))) {
        stop("mismatch between the guide names returned when ", dataset_name, " and ", method_name, " is loaded.")
      }
      all_assns[[dataset_name_short]][[method_name]] <- curr_assns
    }
  }
  all_assns
}

load_assns <- function(path_to_assns, method_name, all_guides, cell_names_in_order) {
  if (method_name == "sceptre") {
    load_assns_sceptre(path_to_assns = path_to_assns, all_guides = all_guides, cell_names_in_order = cell_names_in_order)
  } else if (method_name == "crispat") {
    load_assns_crispat(path_to_assns = path_to_assns, all_guides = all_guides)
  } else if (method_name == "pertpy") {
    load_assns_pertpy(path_to_assns = path_to_assns, all_guides = all_guides)
  } else if (method_name == "cleanser") {
    load_assns_cleanser(path_to_assns = path_to_assns, all_guides = all_guides, cell_names_in_order = cell_names_in_order)
  } else {
    stop("method '", method_name, "' not recognized.")
  }
}

load_assns_sceptre <- function(path_to_assns, all_guides, cell_names_in_order) {
  scep_assn_mat <- read_rds(path_to_assns)
  guide_to_assns <- vector("list", length(all_guides)) |> setNames(all_guides)
  for (guide in all_guides) {
    guide_to_assns[[guide]] <- cell_names_in_order[which(scep_assn_mat[guide, ])]
  }
  guide_to_assns
}

load_assns_crispat <- function(path_to_assns, all_guides) {
  crisp_assign_raw <- read.csv(path_to_assns)
  unordered_cell_list <- split(crisp_assign_raw$cell_id, crisp_assign_raw$grna_id)
  guide_to_assns <- vector("list", length(all_guides)) |> setNames(all_guides)
  guide_to_assns[names(unordered_cell_list)] <- unordered_cell_list
  guide_to_assns
}

load_assns_pertpy <- function(path_to_assns, all_guides) {
  pert_assign_raw <- read.csv(path_to_assns) |> filter(grna_id != "negative")
  guide_to_assns <- vector("list", length(all_guides)) |> setNames(all_guides)
  dt <- as.data.table(pert_assign_raw)
  spl <- dt[, tstrsplit(grna_id, "\\+", fixed = FALSE)]
  long <- melt(
    cbind(dt[, .(cell_id)], spl),
    id.vars = "cell_id", value.name = "grna", na.rm = TRUE
  )[, .(cell_id, grna)]
  setkey(long, grna, cell_id)
  long <- unique(long)
  grna_cells_list <- long[, .(cells = list(cell_id)), by = grna]
  grna_cells_list <- setNames(grna_cells_list$cells, grna_cells_list$grna)
  guide_to_assns[names(grna_cells_list)] <- grna_cells_list
  guide_to_assns
}

load_assns_cleanser <- function(path_to_assns, all_guides, cell_names_in_order, posterior_thresh = 0.8) {
  clns_assign_raw <- read_csv(path_to_assns, show_col_types = FALSE) |> as.data.frame()
  clns_assign_raw <- clns_assign_raw |>
    filter(posterior >= posterior_thresh) |>
    mutate(grna_id = all_guides[grna_id], cell_name = cell_names_in_order[cell_id])
  guide_to_assns <- vector("list", length(all_guides)) |> setNames(all_guides)
  unordered_cell_list <- split(clns_assign_raw$cell_name, clns_assign_raw$grna_id)
  guide_to_assns[names(unordered_cell_list)] <- unordered_cell_list
  guide_to_assns
}

# =============================================================================
# Set-similarity metrics (each on two vectors of cell ids)
# =============================================================================
precision <- function(set1, set2) {
  denom <- length(set1)
  if (denom == 0) return(NA_real_)
  length(intersect(set1, set2)) / denom
}
recall <- function(set1, set2) {
  denom <- length(set2)
  if (denom == 0) return(NA_real_)
  length(intersect(set1, set2)) / denom
}
F1 <- function(set1, set2) {
  denom <- length(set1) + length(set2)
  if (denom == 0) return(NA_real_)
  2 * length(intersect(set1, set2)) / denom
}
jaccard <- function(set1, set2) {
  u <- length(union(set1, set2))
  if (u == 0) return(NA_real_)
  length(intersect(set1, set2)) / u
}

# =============================================================================
# UMI-count histograms: exact singleton bins 0..k_exact, then doubling-width bins.
# Count axis on log1p so the ambient spike at 0 and the sparse tail are both
# legible. bin_umi_counts() builds ONE shared bin structure (global max) so any
# faceting of the resulting counts is aligned by construction (never
# geom_histogram's per-panel binning).
# =============================================================================
logp1 <- function(x, base = 10) log(x + 1, base = base)

make_mixed_bin_info <- function(max_x, k_exact = 10L) {
  stopifnot(length(max_x) == 1, is.finite(max_x), max_x >= 0)
  upper <- 0:k_exact
  width <- 2L
  while (tail(upper, 1L) < max_x) {
    upper <- c(upper, tail(upper, 1L) + width)
    width <- width * 2L
  }
  lower <- c(0L, head(upper, -1L) + 1L)
  label <- ifelse(lower == upper, as.character(upper), paste0(lower, "-", upper))
  data.frame(
    bin_id = seq_along(label), lower = lower, upper = upper, label = label,
    xmin = seq_along(label) - 0.45, xmax = seq_along(label) + 0.45,
    stringsAsFactors = FALSE
  )
}

bin_with_info <- function(x, bin_info) {
  cut(x, breaks = c(-1, bin_info$upper), labels = bin_info$label, right = TRUE)
}

make_y_breaks_umi <- function(max_count) {
  breaks <- c(0, 1, 5, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5, 5e5)
  breaks[breaks <= max_count]
}

# declutter x labels: keep every `exact_every`-th singleton and every `tail_every`-th tail bin
make_x_axis_spec <- function(bin_info, k_exact, exact_every = 2L, tail_every = 2L) {
  is_exact  <- bin_info$lower == bin_info$upper & bin_info$upper <= k_exact
  exact_idx <- which(is_exact)
  tail_idx  <- which(!is_exact)
  keep <- rep(FALSE, nrow(bin_info))
  if (length(exact_idx) > 0) {
    keep[exact_idx[seq(1, length(exact_idx), by = exact_every)]] <- TRUE
    keep[max(exact_idx)] <- TRUE
  }
  if (length(tail_idx) > 0) {
    keep[tail_idx[seq(1, length(tail_idx), by = tail_every)]] <- TRUE
    keep[max(tail_idx)] <- TRUE
  }
  list(breaks = bin_info$bin_id[keep], labels = bin_info$label[keep])
}

# Core: df with `umi` + grouping columns -> per-group binned counts (with bin
# geometry attached), using one shared bin structure across all groups.
bin_umi_counts <- function(df, group_cols, k_exact = 10L) {
  bin_info <- make_mixed_bin_info(max(df$umi), k_exact = k_exact)
  counts <- df |>
    dplyr::mutate(bin = factor(bin_with_info(umi, bin_info), levels = bin_info$label)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), bin, .drop = FALSE) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(bin = as.character(bin)) |>
    dplyr::left_join(bin_info, by = c("bin" = "label"))
  list(counts = counts, bin_info = bin_info)
}

umi_hist_scales <- function(binned, k_exact = 10L, exact_every = 2L, tail_every = 2L) {
  x_axis   <- make_x_axis_spec(binned$bin_info, k_exact, exact_every, tail_every)
  y_breaks <- make_y_breaks_umi(max(binned$counts$count))
  list(
    scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels,
                       expand = expansion(add = 0.6)),
    scale_y_continuous(breaks = logp1(y_breaks),
                       labels = scales::label_comma()(y_breaks),
                       expand = expansion(mult = c(0, 0.05)))
  )
}

# Space-efficient 2-panel overlay: facet by Quantile (low vs high), and within
# each panel overlay two histograms colored by the statistic (Feature) that
# selected the item. `df` needs columns umi, Quantile, Feature (factors). Shared
# bins => aligned; shared x AND y across the two panels.
plot_overlay_sweep_hist <- function(df, title, y_lab = "Number of cells",
                                    pal = feature_pal, alpha = 0.55,
                                    k_exact = 10L, exact_every = 2L, tail_every = 2L) {
  b <- bin_umi_counts(df, c("Quantile", "Feature"), k_exact)
  ggplot(b$counts) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = logp1(count), fill = Feature),
              alpha = alpha, colour = NA) +
    facet_wrap(~ Quantile) +
    # positional: color follows Feature factor order (not the palette names), so
    # reversing the factor levels in the builders swaps the two colors.
    scale_fill_manual(values = unname(pal)) +
    umi_hist_scales(b, k_exact, exact_every, tail_every) +
    labs(title = title, x = "gRNA UMI count", y = y_lab, fill = "Selected by") +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Two-dataset version of plot_overlay_sweep_hist: builds the overlay df for each
# dataset (via `builder`, i.e. build_guide_overlay_df or build_cell_overlay_df) at
# the quantiles `q`, then lays them out with one dataset per row. Facet labels are
# "<Dataset>: <quantile>" (e.g. "Replogle: 10th percentile"). Always
# length(datasets) rows (default 2: Replogle top, Gasperini bottom), quantiles across
# the columns. One shared bin structure (global max over both datasets) => the x-axis
# is aligned across every panel, just like the single-dataset version aligns quantiles.
# `features` selects which selecting-statistic overlays to draw: "both" (Total UMIs +
# Nonzero UMIs, the default), "nonzero" (Nonzero UMIs only), or "total" (Total UMIs
# only).
plot_overlay_sweep_hist_combined <- function(builder, q, title,
                                             y_lab = "Number of cells",
                                             features = c("both", "nonzero", "total"),
                                             datasets = c(Replogle = "replogle-rd7",
                                                          Gasperini = "gasperini"),
                                             pal = feature_pal, alpha = 0.55,
                                             k_exact = 10L, exact_every = 2L, tail_every = 2L) {
  features <- match.arg(features)
  df <- purrr::imap_dfr(datasets, function(ds, nm) {
    builder(ds, q = q) |> dplyr::mutate(Dataset = nm)
  })

  # Name-keyed palette reproducing the current positional mapping (value[i] -> Feature
  # level[i]), so "both" is unchanged and a single-feature subset keeps its own color.
  feat_levels <- levels(df$Feature)
  named_pal   <- stats::setNames(unname(pal)[seq_along(feat_levels)], feat_levels)

  keep <- switch(features,
                 both    = feat_levels,
                 nonzero = grep("Nonzero", feat_levels, value = TRUE),
                 total   = grep("Total",   feat_levels, value = TRUE))
  df <- df |>
    dplyr::filter(Feature %in% keep) |>
    dplyr::mutate(Feature = droplevels(Feature))

  # Combined facet label, ordered so each dataset is a full row (all quantiles across),
  # datasets stacked in the order of `datasets`.
  q_levels     <- levels(df$Quantile)
  facet_levels <- as.vector(outer(q_levels, names(datasets),
                                   function(qq, dd) paste0(dd, ": ", qq)))
  df <- df |>
    dplyr::mutate(facet = factor(paste0(Dataset, ": ", Quantile), levels = facet_levels))

  b <- bin_umi_counts(df, c("facet", "Feature"), k_exact)
  # With a single selecting-statistic, each panel has one overlay and the only thing
  # distinguishing panels is the dataset (the facet rows), so fill by Dataset
  # (dataset_colour) and drop the now-redundant legend. With "both", keep the
  # per-feature fill + "Selected by" legend.
  single <- features != "both"
  b$counts$Dataset <- factor(sub(":.*$", "", as.character(b$counts$facet)),
                             levels = names(datasets))
  fill_var <- if (single) "Dataset" else "Feature"
  p <- ggplot(b$counts) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = logp1(count),
                  fill = .data[[fill_var]]),
              alpha = alpha, colour = NA) +
    facet_wrap(~ facet, nrow = length(datasets)) +
    umi_hist_scales(b, k_exact, exact_every, tail_every) +
    labs(title = title, x = "gRNA UMI count", y = y_lab) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (single) p + scale_fill_manual(values = dataset_colour, guide = "none")
  else        p + scale_fill_manual(values = named_pal, name = "Selected by")
}

# =============================================================================
# §1 EDA data builders
# =============================================================================

# English ordinal suffix for an integer: 1st, 2nd, 3rd, 4th, ... with the 11th/12th/
# 13th exceptions (and 21st, 22nd, 23rd, ...). Vectorized.
ordinal_suffix <- function(n) {
  n <- as.integer(round(n))
  ifelse(n %% 100 %in% 11:13, "th",
         ifelse(n %% 10 == 1, "st",
                ifelse(n %% 10 == 2, "nd",
                       ifelse(n %% 10 == 3, "rd", "th"))))
}

# Factor labeller for quantile levels, e.g. 0.25 -> "25th percentile", 0.01 -> "1st
# percentile". Ordinal suffix is computed (not hardcoded "th"), so low quantiles read
# correctly.
make_q_lab <- function(q) {
  fmt      <- function(p) { n <- round(p * 100); paste0(n, ordinal_suffix(n), " percentile") }
  levels_  <- fmt(sort(q))
  function(p) factor(fmt(p), levels = levels_)
}

# guide/cell row whose single feature is closest to its p-quantile
pick_by_feature_quantile <- function(feature, p) {
  which.min(abs(feature - quantile(feature, p, names = FALSE, type = 7)))
}

# Fig 1a (guides): 2-panel overlay. Guides are picked SEPARATELY by each statistic
# (grna_nonzero, grna_total) — jointly-extreme guides like "high nonzero, low
# total" don't exist. Each panel (a quantile level) overlays the guide chosen by
# grna_nonzero with the one chosen by grna_total, colored by Feature.
build_guide_overlay_df <- function(data_source, q = c(0.1, 0.9)) {
  path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", data_source, "sceptre-pipeline")
  grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path, "grna.odm"))
  stats <- read_csv(file.path(path, "grna_summary_stats.csv"), show_col_types = FALSE)
  q_lab <- make_q_lab(q)

  feature_cols <- c(`Total UMIs` = "grna_total", `Nonzero UMIs` = "grna_nonzero")
  spec <- expand.grid(p = q, feature = names(feature_cols), stringsAsFactors = FALSE)

  purrr::pmap_dfr(spec, function(p, feature) {
    i <- pick_by_feature_quantile(stats[[feature_cols[[feature]]]], p)
    tibble(
      umi      = grna_odm[stats$grna[i], ],
      Quantile = q_lab(p),
      Feature  = factor(feature, levels = names(feature_cols))
    )
  })
}

# Fig 1b (cells): 2-panel overlay, mirroring build_guide_overlay_df. Cells are
# picked SEPARATELY by each per-cell statistic (grna_n_nonzero, grna_n_umis) from
# the sceptre_object covariate data frame. Each panel (a quantile level) overlays
# the cell chosen by grna_n_nonzero with the one chosen by grna_n_umis, colored by
# Feature. odm access is row-wise (guide-wise), so we loop over every guide and
# pull just the picked-cell entries (the slow part — cache the chunks if it drags).
build_cell_overlay_df <- function(data_source, q = c(0.1, 0.9)) {
  path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", data_source, "sceptre-pipeline")
  scep <- sceptre::read_ondisc_backed_sceptre_object(
    sceptre_object_fp    = file.path(path, "sceptre_object.rds"),
    response_odm_file_fp = file.path(path, "response.odm"),
    grna_odm_file_fp     = file.path(path, "grna.odm")
  )
  cov   <- scep@covariate_data_frame
  q_lab <- make_q_lab(q)

  feature_cols <- c(`Total UMIs` = "grna_n_umis", `Nonzero UMIs` = "grna_n_nonzero")
  spec <- expand.grid(p = q, feature = names(feature_cols), stringsAsFactors = FALSE)
  spec$cell_idx <- purrr::map2_int(spec$feature, spec$p, function(feature, p) {
    pick_by_feature_quantile(cov[[feature_cols[[feature]]]], p)
  })

  grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path, "grna.odm"))
  n_guides <- nrow(grna_odm)
  cols <- matrix(0L, nrow = n_guides, ncol = nrow(spec))  # guides x picked cells
  for (j in seq_len(n_guides)) {
    cols[j, ] <- grna_odm[j, ][spec$cell_idx]
  }

  purrr::map_dfr(seq_len(nrow(spec)), function(k) {
    tibble(
      umi      = cols[, k],
      Quantile = q_lab(spec$p[k]),
      Feature  = factor(spec$feature[k], levels = names(feature_cols))
    )
  })
}

# =============================================================================
# §1 EDA — cell-level gRNA covariates (hexbins)
# =============================================================================
# CLEANSER's per-cell "low-pass" depth: the sum of a cell's gRNA counts that are
# <= 2 (its signal-free ambient scaling factor). Matches datasets.R::kappa_leq2 —
# zero out entries > 2, then colSums.
cleanser_low_pass <- function(grna_mat) {
  s <- grna_mat
  s@x[s@x > 2] <- 0
  Matrix::colSums(s)
}

# Per-cell gRNA covariates for one base dataset ("replogle-rd7" / "gasperini"),
# read from its ondisc-backed sceptre object. grna_n_umis and grna_n_nonzero come
# straight from scep@covariate_data_frame (nothing recomputed); low_pass is the
# only computed quantity — CLEANSER's <=2 depth, accumulated over the odm's guide
# rows so it aligns to the covariate rows by cell position. The odm loop is the
# slow part (~10s Replogle, ~45s Gasperini); cache the chunk that calls this.
build_grna_covariate_df <- function(data_source) {
  path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", data_source, "sceptre-pipeline")
  scep <- sceptre::read_ondisc_backed_sceptre_object(
    sceptre_object_fp    = file.path(path, "sceptre_object.rds"),
    response_odm_file_fp = file.path(path, "response.odm"),
    grna_odm_file_fp     = file.path(path, "grna.odm")
  )
  cov <- scep@covariate_data_frame

  grna_odm <- ondisc::initialize_odm_from_backing_file(file.path(path, "grna.odm"))
  ncells   <- ncol(grna_odm)
  stopifnot(nrow(cov) == ncells)          # covariate rows align to odm cells by position
  low_pass <- numeric(ncells)
  for (j in seq_len(nrow(grna_odm))) {    # accumulate the <=2 depth guide-by-guide
    row <- grna_odm[j, ]                  # dense integer counts, length ncells
    w   <- row > 0L & row <= 2L
    low_pass[w] <- low_pass[w] + row[w]
  }
  tibble::tibble(
    grna_n_umis    = cov$grna_n_umis,
    grna_n_nonzero = cov$grna_n_nonzero,
    low_pass       = low_pass
  )
}

# Single-panel hexbin of CLEANSER's low-pass depth (sum of counts <= 2) vs total
# gRNA UMIs per cell. Uses build_grna_covariate_df() output (needs low_pass, so the
# slow odm loop). Hex fill is the cell count (log10); non-positive dropped on log axes.
plot_low_pass_hex <- function(df, dataset_label, bins = 60, log_axes = TRUE) {
  d <- df
  if (log_axes) d <- dplyr::filter(d, grna_n_umis > 0, low_pass > 0)
  x_lab <- if (log_axes) "Total gRNA UMIs per cell (log scale)" else "Total gRNA UMIs per cell"
  y_lab <- if (log_axes) "CLEANSER low-pass depth (sum of counts ≤ 2, log scale)" else
                         "CLEANSER low-pass depth (sum of counts ≤ 2)"
  p <- ggplot(d, aes(grna_n_umis, low_pass)) +
    geom_hex(bins = bins) +
    scale_fill_viridis_c(trans = "log10", name = "Cell count",
                         labels = scales::label_number(big.mark = ",")) +
    labs(title = paste0(dataset_label, ": CLEANSER low-pass depth vs. total gRNA UMIs"),
         x = x_lab, y = y_lab) +
    theme_thesis()
  if (log_axes) {
    p <- p +
      scale_x_log10(labels = scales::label_number(big.mark = ",")) +
      scale_y_log10(labels = scales::label_number(big.mark = ","))
  }
  p
}

# Per-cell covariates (all four count covariates) for one base dataset, read straight
# from scep@covariate_data_frame — nothing recomputed, no odm loop, so this is fast.
build_cell_covariate_df <- function(data_source) {
  path <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                    "guide_assignment/input_data", data_source, "sceptre-pipeline")
  scep <- sceptre::read_ondisc_backed_sceptre_object(
    sceptre_object_fp    = file.path(path, "sceptre_object.rds"),
    response_odm_file_fp = file.path(path, "response.odm"),
    grna_odm_file_fp     = file.path(path, "grna.odm")
  )
  cov <- scep@covariate_data_frame
  tibble::tibble(
    grna_n_nonzero     = cov$grna_n_nonzero,
    grna_n_umis        = cov$grna_n_umis,
    response_n_nonzero = cov$response_n_nonzero,
    response_n_umis    = cov$response_n_umis
  )
}

# 2x2 hexbin grid of cell-level covariates (one figure per dataset). Each panel is a
# distinct (x, y) covariate pair with its own axis titles (via patchwork, so every
# panel labels both axes). Top row: x = num. expressed gRNAs (y = total gRNA
# expression, num. expressed genes). Bottom row: x = num. expressed genes (y = total
# gene expression, total gRNA expression). Free scales per panel; one shared
# cell-count fill (log10) with a common limit so the legend collapses to one.
# Non-positive values are dropped on the log axes.
plot_covariate_grid <- function(df, dataset_label, bins = 50, log_axes = TRUE) {
  cov_labs <- c(
    grna_n_nonzero     = "Num. expressed gRNAs",
    grna_n_umis        = "Total gRNA expression",
    response_n_nonzero = "Num. expressed genes",
    response_n_umis    = "Total gene expression"
  )
  specs <- list(
    c(x = "grna_n_nonzero",     y = "grna_n_umis"),
    c(x = "grna_n_nonzero",     y = "response_n_nonzero"),
    c(x = "response_n_nonzero", y = "response_n_umis"),
    c(x = "response_n_nonzero", y = "grna_n_umis")
  )

  base_panel <- function(s) {
    d <- tibble::tibble(x = df[[s[["x"]]]], y = df[[s[["y"]]]])
    if (log_axes) d <- dplyr::filter(d, x > 0, y > 0)
    p <- ggplot(d, aes(x, y)) +
      geom_hex(bins = bins) +
      labs(x = cov_labs[[s[["x"]]]], y = cov_labs[[s[["y"]]]]) +
      theme_thesis()
    if (log_axes) {
      p <- p +
        scale_x_log10(labels = scales::label_number(big.mark = ",")) +
        scale_y_log10(labels = scales::label_number(big.mark = ","))
    }
    p
  }
  panels <- lapply(specs, base_panel)

  # Common fill limit across panels -> guides="collect" yields a single legend.
  gmax <- max(vapply(panels, function(p)
    max(ggplot2::ggplot_build(p)$data[[1]]$count, na.rm = TRUE), numeric(1)))
  fill_scale <- scale_fill_viridis_c(trans = "log10", name = "Cell count",
                                     limits = c(1, gmax),
                                     labels = scales::label_number(big.mark = ","))
  panels <- lapply(panels, function(p) p + fill_scale)

  patchwork::wrap_plots(panels, nrow = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(dataset_label, ": cell-level gRNA and gene expressions"),
      theme = theme_thesis()
    ) &
    theme(legend.position = "bottom")
}

# Simpler 1x2 hexbin: num. expressed vs total expression, once for the gRNA
# covariates (left) and once for the response covariates (right). Uses the fast
# build_cell_covariate_df() output (no odm loop). Free scales per panel; one
# shared cell-count fill (log10) with a common limit -> single collected legend.
# Non-positive values dropped on the log axes.
plot_covariate_pair <- function(df, dataset_label, bins = 50, log_axes = TRUE) {
  cov_labs <- c(
    grna_n_nonzero     = "Num. expressed gRNAs",
    grna_n_umis        = "Total gRNA expression",
    response_n_nonzero = "Num. expressed genes",
    response_n_umis    = "Total gene expression"
  )
  specs <- list(
    c(x = "grna_n_nonzero",     y = "grna_n_umis"),
    c(x = "response_n_nonzero", y = "response_n_umis")
  )

  base_panel <- function(s) {
    d <- tibble::tibble(x = df[[s[["x"]]]], y = df[[s[["y"]]]])
    if (log_axes) d <- dplyr::filter(d, x > 0, y > 0)
    p <- ggplot(d, aes(x, y)) +
      geom_hex(bins = bins) +
      labs(x = cov_labs[[s[["x"]]]], y = cov_labs[[s[["y"]]]]) +
      theme_thesis()
    if (log_axes) {
      p <- p +
        scale_x_log10(labels = scales::label_number(big.mark = ",")) +
        scale_y_log10(labels = scales::label_number(big.mark = ","))
    }
    p
  }
  panels <- lapply(specs, base_panel)

  gmax <- max(vapply(panels, function(p)
    max(ggplot2::ggplot_build(p)$data[[1]]$count, na.rm = TRUE), numeric(1)))
  fill_scale <- scale_fill_viridis_c(trans = "log10", name = "Cell count",
                                     limits = c(1, gmax),
                                     labels = scales::label_number(big.mark = ","))
  panels <- lapply(panels, function(p) p + fill_scale)

  patchwork::wrap_plots(panels, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(dataset_label, ": cell-level gRNA and gene expressions"),
      theme = theme_thesis()
    ) &
    theme(legend.position = "bottom")
}

# =============================================================================
# §2 Simulations with ground truth
# =============================================================================

# Load the simulation bundle once: per-method assignments, the ground-truth
# perturbation matrix, and the simulated gRNA count matrix.
load_sim_bundle <- function() {
  run_names_sims     <- "run_sims_4mu_50k_cleanser_params"
  dataset_names_sims <- c(sims = "replogle-rd7_sims_from_clns_50k_4_mu")
  base <- .get_config_path("LOCAL_BENCHMARKING_DIR")
  list(
    sims_assns = load_all_assns(run_names = run_names_sims,
                                dataset_names = dataset_names_sims, methods = methods),
    true_perts = read_rds(file.path(base, "guide_assignment/input_data",
                                    dataset_names_sims[["sims"]], "true_pert_matrix.rds")),
    sim_data   = read_rds(file.path(base, "guide_assignment/input_data",
                                    dataset_names_sims[["sims"]], "sceptre/grna_matrix.rds"))
  )
}

# Fig 2a: long data frame of one or more simulated guides' UMI counts, split by
# ground-truth perturbation status. `guides` is a NAMED character vector whose
# names become the panel labels (e.g. c(`Perturbed Mean = 971` = "..._grna_5")).
build_sim_component_df <- function(sim_data, true_perts, guides) {
  purrr::imap_dfr(guides, function(g, nm) {
    tibble(
      umi         = as.integer(sim_data[g, ]),
      pert_status = factor(ifelse(true_perts[g, ] == 1, "Pert.", "Ambient"),
                           levels = c("Ambient", "Pert.")),
      panel       = factor(nm, levels = names(guides))
    )
  })
}

# Fig 2a: stacked mixed-bin histogram — ambient on the bottom, perturbed stacked
# on top, per UMI-count bin. One shared bin structure (global max) across panels.
plot_sim_component_hist <- function(df, title = NULL, y_lab = "Number of cells",
                                    pal = sim_status_pal, k_exact = 10L,
                                    exact_every = 2L, tail_every = 2L) {
  bin_info <- make_mixed_bin_info(max(df$umi), k_exact = k_exact)
  counts <- df |>
    dplyr::mutate(bin = factor(bin_with_info(umi, bin_info), levels = bin_info$label)) |>
    dplyr::group_by(panel, bin, pert_status, .drop = FALSE) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(panel, bin, pert_status) |>
    dplyr::group_by(panel, bin) |>
    dplyr::mutate(y1 = cumsum(count), y0 = dplyr::lag(y1, default = 0)) |>
    dplyr::ungroup() |>
    dplyr::mutate(bin = as.character(bin)) |>
    dplyr::left_join(bin_info, by = c("bin" = "label"))

  x_axis   <- make_x_axis_spec(bin_info, k_exact, exact_every, tail_every)
  y_breaks <- make_y_breaks_umi(max(counts$y1))

  ggplot(counts) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = logp1(y0), ymax = logp1(y1),
                  fill = pert_status), colour = "black", linewidth = 0.15) +
    facet_wrap(~ panel) +
    scale_fill_manual(values = pal) +
    scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels,
                       expand = expansion(add = 0.6)) +
    scale_y_continuous(breaks = logp1(y_breaks),
                       labels = scales::label_comma()(y_breaks),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = "gRNA UMI count", y = y_lab, fill = NULL) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Fig 2b/2c: per-(guide, method) Jaccard / Precision / Recall against ground truth.
# Ground truth for a guide = truly-perturbed cells that have a nonzero UMI (a
# perturbed cell with zero UMI carries no signal any method could use). Metrics are
# computed in one pass (one intersect per guide-method), scored on 100 guides per
# perturbation level. NA when a denominator is empty (e.g. a method assigns 0).
sim_per_guide_scores <- function(sim, mupow_seq = 0:3, nguides_per = 100) {
  sims_assns <- sim$sims_assns$sims
  true_perts <- sim$true_perts
  sim_data   <- sim$sim_data
  cell_names <- colnames(true_perts)
  purrr::map_dfr(mupow_seq, function(mupow) {
    purrr::map_dfr(seq_len(nguides_per), function(gi) {
      g     <- paste0("params_cleanser_mupow=", mupow, "_grna_", gi)
      truth <- cell_names[true_perts[g, ] == 1 & sim_data[g, ] > 0]
      nt    <- length(truth)
      purrr::map_dfr(methods, function(m) {
        pred  <- sims_assns[[m]][[g]]
        np    <- length(pred)
        inter <- length(intersect(pred, truth))
        uni   <- np + nt - inter
        tibble(
          mupow = mupow, guide_id = gi, method = m,
          Jaccard   = if (uni == 0) NA_real_ else inter / uni,
          Precision = if (np  == 0) NA_real_ else inter / np,
          Recall    = if (nt  == 0) NA_real_ else inter / nt
        )
      })
    })
  })
}

# Collapse per-guide scores to a tidy mean per (perturbation level, method, metric).
sim_avg_metrics <- function(per_guide) {
  per_guide |>
    tidyr::pivot_longer(c(Jaccard, Precision, Recall),
                        names_to = "metric", values_to = "score") |>
    dplyr::group_by(mupow, method, metric) |>
    dplyr::summarise(value = mean(score, na.rm = TRUE), .groups = "drop")
}

# Fig 2b: mean metric vs perturbation strength (mu_1 = 971 / 2^mupow), one line
# per method, faceted by metric.
plot_sim_avg_metrics <- function(avg_df, title = NULL) {
  avg_df |>
    dplyr::mutate(
      method  = factor(method_renamer[method], levels = method_renamer),
      metric  = factor(metric, levels = c("Jaccard", "Precision", "Recall")),
      mu_pert = 971 / 2^mupow
    ) |>
    ggplot(aes(mu_pert, value, colour = method, group = method)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    facet_wrap(~ metric) +
    scale_x_log10(breaks = round(971 / 2^(0:3))) +
    scale_colour_manual(values = method_pal) +
    labs(title = title,
         x = expression("Perturbed-cell NegBin mean" ~ (mu[1])),
         y = "Mean metric value", colour = "Method") +
    theme_thesis()
}

# Fig 2c: distribution of per-guide Jaccard (to ground truth) across guides, one
# panel per perturbation level with the four methods STACKED (overlapping 4
# histograms is unreadable). Panels ordered so mu increases left -> right (mu =
# 971 / 2^mupow, so ascending mu == descending mupow). Shows how per-guide accuracy
# spreads out as the perturbation weakens, and where methods fail (mass at 0).
plot_sim_jaccard_hist <- function(per_guide, mupows = 0:3, bins = 30, title = NULL) {
  levs <- paste0("Perturbed mean = ", round(971 / 2^sort(mupows, decreasing = TRUE)))
  per_guide |>
    dplyr::filter(mupow %in% mupows, !is.na(Jaccard)) |>
    dplyr::mutate(
      method = factor(method_renamer[method], levels = method_renamer),
      panel  = factor(paste0("Perturbed mean = ", round(971 / 2^mupow)), levels = levs)
    ) |>
    ggplot(aes(x = Jaccard, fill = method)) +
    geom_histogram(bins = bins, colour = "black", linewidth = 0.1) +
    facet_wrap(~ panel, nrow = 1) +
    scale_fill_manual(values = method_pal) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = title, x = "Jaccard similarity to ground truth",
         y = "Number of guides", fill = "Method") +
    theme_thesis()
}

# Fig 2d (and §3): number of assignments per guide vs the guide's "tail mean" (the
# mean UMI among counts >= thresh, a proxy for how heavy-tailed the guide is), for
# each method. General over any (grna matrix, per-method assignment list) so it
# serves both the simulations and the real datasets.
build_tail_sensitivity_df <- function(grna_mat, assns_list, thresh = 10) {
  guides <- rownames(grna_mat)
  tail_mean <- sapply(guides, function(g) {
    row  <- grna_mat[g, ]
    tail <- row[row >= thresh]
    if (length(tail) == 0) 0 else mean(tail)
  })
  # Parse perturbation strength from the sim guide names
  # ("params_cleanser_mupow=<j>_grna_<n>"); NA for real-data guides.
  mupow <- suppressWarnings(as.integer(sub(".*mupow=([0-9]+).*", "\\1", guides)))
  mu    <- ifelse(grepl("mupow=", guides), round(971 / 2^mupow), NA_real_)
  names(mu) <- guides
  purrr::map_dfr(names(assns_list), function(m) {
    tibble(
      guide     = guides,
      method    = m,
      mu        = mu[guides],
      tail_mean = tail_mean[guides],
      n_assn    = sapply(assns_list[[m]][guides], length)
    )
  })
}

# Combine both datasets' tail-sensitivity data into one long df (adds a `dataset`
# column), for the overlaid version below. `real` is the load_real_bundle() list.
build_tail_sensitivity_combined_df <- function(real, thresh = 10,
    datasets = c(Replogle = "replogle-rd7", Gasperini = "gasperini")) {
  purrr::imap_dfr(datasets, function(key, label)
    build_tail_sensitivity_df(real$grna_mat[[key]], real$assns[[key]], thresh) |>
      dplyr::mutate(dataset = label))
}

# Overlaid version of plot_tail_sensitivity: both datasets in the SAME panel per
# method (2x2 over the four methods), coloured by dataset (dataset_colour: Replogle
# blue, Gasperini orange). x is log1p too (not just y) because the datasets' tail
# means span very different ranges -- Replogle's tails are far heavier -- so a
# linear x would crush Gasperini into the left edge.
plot_tail_sensitivity_combined <- function(df, title = NULL, thresh = 10,
                                           pal = dataset_colour) {
  df <- dplyr::mutate(df,
    method  = factor(method_renamer[method], levels = method_renamer),
    dataset = factor(dataset, levels = names(pal)))
  ggplot(df, aes(tail_mean, n_assn, colour = dataset)) +
    geom_point(alpha = 0.4, size = 0.9) +
    facet_wrap(~ method) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000, 100000)) +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000, 10000)) +
    scale_colour_manual(values = pal, name = "Dataset") +
    labs(title = title,
         x = paste0("Mean UMI among counts ≥ ", thresh, " (log1p)"),
         y = "Assignments per guide") +
    theme_thesis() +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
}

# Dataset label -> guide_assignment input_data dir key (matches load_real_bundle).
.dataset_dir <- c(Replogle = "replogle-rd7", Gasperini = "gasperini")

# Grid of single-guide UMI histograms, for showcasing guides where a method went
# wrong. `picks` is a list of list(dataset = "Replogle"/"Gasperini", guide = <id>);
# `titles` gives one custom panel title per pick (same length/order; defaults to
# the guide id). Each histogram uses the mixed-bin linear-then-log style over the
# guide's NONZERO counts with its own bins (guides span very different maxima), and
# bars are coloured by dataset (dataset_colour). Arranged with `ncol` (default 2 ->
# 2x2 for four picks). No faceting: each panel is an independent plot.
plot_guide_umi_grid <- function(real, picks, titles = NULL, ncol = 2, k_exact = 10L,
                                pal = dataset_colour) {
  if (is.null(titles)) titles <- vapply(picks, function(p) p$guide, character(1))
  stopifnot(length(titles) == length(picks))
  plots <- purrr::map2(picks, titles, function(p, ttl) {
    row <- real$grna_mat[[.dataset_dir[[p$dataset]]]][p$guide, ]
    nz  <- as.integer(row[row > 0])
    b   <- bin_umi_counts(data.frame(umi = nz, grp = "g"), "grp", k_exact)
    ggplot(b$counts) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = logp1(count)),
                fill = pal[[p$dataset]], colour = NA) +
      umi_hist_scales(b, k_exact) +
      labs(title = ttl, x = "gRNA UMI count (nonzero cells)", y = "Number of cells (log)") +
      theme_thesis() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
            plot.title  = element_text(size = 10))
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

plot_tail_sensitivity <- function(df, title = NULL, thresh = 10,
                                  colour_by = c("method", "mu")) {
  colour_by <- match.arg(colour_by)
  df <- dplyr::mutate(df, method = factor(method_renamer[method],
                                          levels = method_renamer))
  if (colour_by == "mu") {
    # mu ascending (weak -> strong) so the legend orders sensibly; discrete factor.
    df$mu <- factor(df$mu, levels = sort(unique(df$mu)))
    p <- ggplot(df, aes(tail_mean, n_assn, colour = mu)) +
      geom_point(alpha = 0.4, size = 0.9) +
      scale_colour_manual(values = mu_pal, name = "Perturbed mean")
  } else {
    p <- ggplot(df, aes(tail_mean, n_assn, colour = method)) +
      geom_point(alpha = 0.4, size = 0.9) +
      scale_colour_manual(values = method_pal, guide = "none")
  }
  p +
    facet_wrap(~ method) +
    scale_y_continuous(trans = "log1p", breaks = c(0, 1, 10, 100, 1000, 10000)) +
    labs(title = title,
         x = paste0("Mean UMI among counts ≥ ", thresh),
         y = "Assignments per guide") +
    theme_thesis()
}

# Per-guide calling boundary: for each (guide, method), the smallest UMI the method
# calls perturbed (min over its assigned cells) and the largest UMI it leaves
# unperturbed (max over its unassigned cells). Where largest_not > smallest_pert
# the decision boundary is non-monotone in count (leaning on covariates). Returns
# one row per (mupow, method, guide); mu = 971/2^mupow. Feeds both the table and
# the Fig 2e scatter.
sim_threshold_per_guide <- function(sim, mupow_seq = 0:3, nguides_per = 100) {
  sim_data <- sim$sim_data
  assns    <- sim$sims_assns$sims
  cells    <- colnames(sim_data)

  purrr::map_dfr(mupow_seq, function(mupow) {
    purrr::map_dfr(seq_len(nguides_per), function(gi) {
      g   <- paste0("params_cleanser_mupow=", mupow, "_grna_", gi)
      row <- sim_data[g, ]
      purrr::map_dfr(methods, function(m) {
        is_a <- cells %in% assns[[m]][[g]]
        tibble(
          mupow = mupow, mu = round(971 / 2^mupow), guide = gi, method = m,
          smallest_pert = if (any(is_a))  min(row[is_a])  else NA_real_,
          largest_not   = if (any(!is_a)) max(row[!is_a]) else NA_real_
        )
      })
    })
  })
}

# Table: per (method, perturbation level), averaged over the guides at that level —
# the average smallest UMI a method calls perturbed and the average largest UMI it
# leaves unperturbed. A proxy for each method's effective calling threshold and how
# much the two components overlap. Returns a display-ready tibble (mu ascending).
sim_threshold_table <- function(sim, mupow_seq = 0:3, nguides_per = 100) {
  sim_threshold_per_guide(sim, mupow_seq, nguides_per) |>
    dplyr::group_by(mupow, method) |>
    dplyr::summarise(
      avg_smallest_pert = mean(smallest_pert, na.rm = TRUE),
      avg_largest_not   = mean(largest_not,   na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      mu_pert = round(971 / 2^mupow),
      method  = factor(method_renamer[method], levels = method_renamer)
    ) |>
    dplyr::arrange(mu_pert, method) |>
    dplyr::select(`Perturbed mean`               = mu_pert,
                  Method                         = method,
                  `Avg. smallest perturbed UMI`  = avg_smallest_pert,
                  `Avg. largest unperturbed UMI` = avg_largest_not)
}

# Fig 2e: the calling boundary as boxplots. x = method; within each method two
# dodged boxplots — the largest UMI a method leaves unperturbed and the smallest it
# calls perturbed, across the 100 guides. Whiskers span the full range (coef = Inf),
# so each box shows exactly min / Q1 / median / Q3 / max. Faceted by mu (2x2), so a
# panel holds 8 boxplots. When the two boxes overlap, the method's boundary is not a
# clean count cutoff.
plot_sim_threshold_box <- function(per_guide, title = NULL) {
  # Precompute the five-number summary per box and draw with stat = "identity".
  # (geom_boxplot(coef = Inf) breaks when a group has zero IQR: Inf * 0 = NaN.)
  d <- per_guide |>
    tidyr::pivot_longer(c(largest_not, smallest_pert),
                        names_to = "quantity", values_to = "umi") |>
    dplyr::filter(!is.na(umi)) |>   # guides where a method assigned all / no cells
    dplyr::mutate(
      method   = factor(method_renamer[method], levels = method_renamer),
      quantity = factor(quantity,
                        levels = c("largest_not", "smallest_pert"),
                        labels = names(boundary_pal)),
      mu       = factor(mu, levels = sort(unique(mu)))   # ascending, left -> right
    ) |>
    dplyr::group_by(mu, method, quantity) |>
    dplyr::summarise(
      ymin   = min(umi),
      lower  = stats::quantile(umi, 0.25, names = FALSE),
      middle = stats::median(umi),
      upper  = stats::quantile(umi, 0.75, names = FALSE),
      ymax   = max(umi),
      .groups = "drop"
    )
  ggplot(d, aes(x = method, fill = quantity)) +
    geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle,
                     upper = upper, ymax = ymax),
                 stat = "identity", linewidth = 0.3,
                 width = 0.6, position = position_dodge(width = 0.7)) +
    facet_wrap(~ mu, labeller = labeller(
      mu = function(x) paste0("Perturbed mean = ", x))) +
    scale_fill_manual(values = boundary_pal, name = NULL) +
    labs(title = title, x = "Method", y = "gRNA UMI count") +
    theme_thesis()
}

# =============================================================================
# Section 3 — Real-data agreement between methods
# =============================================================================

# Real-data run names holding the assignment outputs (across the guide-count sweep).
run_names_real <- c("run_all_3guides_100k", "run_all_nguides=800_100k",
                    "rerun_pertpy_more_samples")

# Load the real-data assignments + gRNA count matrices for one guide-count config
# (default 800, the richest). Returns per-dataset method -> guide -> cells lists and
# the raw gRNA matrices (for the tail-sensitivity figures). Only the requested
# config is loaded, so this stays fast.
load_real_bundle <- function(nguides = 800L) {
  ds_repl <- c(`replogle-rd7` = sprintf("replogle-rd7_assign_nguides=%d_ncells=100k", nguides))
  ds_gasp <- c(gasperini      = sprintf("gasperini_assign_nguides=%d_ncells=100k", nguides))
  repl_assns <- load_all_assns(run_names_real, ds_repl, methods)[[1]]
  gasp_assns <- load_all_assns(run_names_real, ds_gasp, methods)[[1]]
  base_in <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                       "guide_assignment/input_data")
  grna_repl <- read_rds(file.path(base_in, ds_repl[[1]], "sceptre/grna_matrix.rds"))
  grna_gasp <- read_rds(file.path(base_in, ds_gasp[[1]], "sceptre/grna_matrix.rds"))
  list(
    nguides  = nguides,
    assns    = list(`replogle-rd7` = repl_assns, gasperini = gasp_assns),
    grna_mat = list(`replogle-rd7` = grna_repl,  gasperini = grna_gasp)
  )
}

# Mean per-guide Jaccard between every pair of methods for one dataset. `assns` is
# method -> guide -> cells. Two methods agreeing a guide has no cells counts as a
# perfect match (union 0 -> 1), matching the defense definition.
mean_jaccard_matrix <- function(assns) {
  ms <- names(assns)
  outer_fn <- function(m1, m2) {
    common <- intersect(names(assns[[m1]]), names(assns[[m2]]))
    mean(vapply(common, function(g) {
      a <- assns[[m1]][[g]]; b <- assns[[m2]][[g]]
      uni <- length(union(a, b))
      if (uni == 0) 1 else length(intersect(a, b)) / uni
    }, numeric(1)))
  }
  m <- outer(ms, ms, Vectorize(outer_fn))
  dimnames(m) <- list(ms, ms)
  m
}

# Assemble the split-triangle heatmap frame: upper triangle = `upper` dataset,
# lower triangle = `lower` dataset, diagonal left NA. Methods ordered by `methods`.
build_jaccard_heatmap_df <- function(real, upper = "replogle-rd7",
                                     lower = "gasperini") {
  mat_u <- mean_jaccard_matrix(real$assns[[upper]])
  mat_l <- mean_jaccard_matrix(real$assns[[lower]])
  grid <- expand.grid(Row = methods, Col = methods, stringsAsFactors = FALSE)
  grid$value <- mapply(function(r, c) {
    ri <- match(r, methods); ci <- match(c, methods)
    if (ri < ci)      mat_u[r, c]   # above diagonal -> upper dataset
    else if (ri > ci) mat_l[r, c]   # below diagonal -> lower dataset
    else              NA_real_
  }, grid$Row, grid$Col)
  # Prettify + order: Col left->right in method order, Row reversed so SCEPTRE on top.
  grid$Row <- factor(unname(method_renamer[grid$Row]),
                     levels = rev(unname(method_renamer[methods])))
  grid$Col <- factor(unname(method_renamer[grid$Col]),
                     levels = unname(method_renamer[methods]))
  grid
}

# Fig 3a: split-triangle mean-Jaccard heatmap. Beige -> red fill, cell values
# printed, one dataset per triangle (labelled in the two free diagonal corners).
plot_jaccard_heatmap <- function(df, upper_label = "Replogle",
                                 lower_label = "Gasperini", title = NULL) {
  n <- nlevels(df$Col)
  ggplot(df, aes(Col, Row, fill = value)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value))),
              colour = "black", size = 3) +
    scale_fill_gradient(low = "#F5F5DC", high = "#F8766D",
                        na.value = "white", limits = c(0, 1)) +
    # Dataset tags in the two empty diagonal corners.
    annotate("text", x = 1, y = n, label = paste0(lower_label, "\n(lower)"),
             size = 3, fontface = "italic", colour = "grey25", lineheight = 0.9) +
    annotate("text", x = n, y = 1, label = paste0(upper_label, "\n(upper)"),
             size = 3, fontface = "italic", colour = "grey25", lineheight = 0.9) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_thesis() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

# =============================================================================
# Section 4 — Computational benchmarking (runtime + memory)
# =============================================================================
# Two execution modes are compared: the four assignment methods run in-memory
# (a guide-count sweep), and the SCEPTRE pipeline run distributed over workers
# (a single point per dataset). "pipeline" gets its own color and, in the plots,
# a triangle marker vs the circle used for the in-memory methods.
method_renamer_comp <- c(method_renamer, pipeline = "pipeline")

# In-memory trace files (the four methods across the guide-count sweep). Parses
# peak_rss to GB and realtime (Hh Mm Ss) to seconds.
load_trace_files_in_memory_assignment <- function(run_names, path_to_outputs) {
  lapply(run_names, function(run_name) {
    read_tsv(file.path(path_to_outputs, run_name, "trace.tsv"), show_col_types = FALSE)
  }) |>
    do.call(what = rbind) |>
    tidyr::extract(tag, into = "nguides", regex = "nguides=([0-9]+)",
                   remove = FALSE, convert = TRUE) |>
    mutate(source_data = strsplit(tag, "_") |> sapply(`[`, 1) |>
             gsub(pattern = "-rd7", replacement = "", fixed = TRUE)) |>
    transmute(
      method = gsub("_assign", "", process, ignore.case = TRUE) |> tolower(),
      peak_rss, realtime,
      num_cells = sub(".*ncells=([^_]+).*", "\\1", tag), nguides, source_data
    ) |>
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
    select(-h, -m, -s)
}

# Distributed pipeline trace (assign_grnas workers from the association compute
# runs). Elapsed time = span from first submit to last (drift-corrected) complete;
# memory = max peak_rss over workers.
load_trace_files_pipeline_assignment <- function(assoc_datasets) {
  parse_mem_to_gb <- function(x) {
    x <- str_trim(x)
    value <- as.numeric(str_extract(x, "^[0-9.]+"))
    unit  <- str_extract(x, "[A-Za-z]+$") |> toupper()
    mult <- case_when(
      unit == "B" ~ 1 / 1024^3, unit == "KB" ~ 1 / 1024^2, unit == "MB" ~ 1 / 1024,
      unit == "GB" ~ 1, unit == "TB" ~ 1024, TRUE ~ NA_real_
    )
    ifelse(is.na(value * mult), 0, value * mult)  # unitless 0s -> 0
  }
  lapply(seq_along(assoc_datasets), function(i) {
    curr_trace <- file.path(
      .get_config_path("LOCAL_BENCHMARKING_DIR"),
      "association/computational/outputs/sceptre-pipeline",
      assoc_datasets[[i]], "tracing/trace.tsv"
    ) |>
      read_tsv(show_col_types = FALSE) |>
      mutate(method = "pipeline", source_data = names(assoc_datasets)[i],
             start    = as.POSIXct(start,    tz = "America/New_York"),
             complete = as.POSIXct(complete, tz = "America/New_York"),
             submit   = as.POSIXct(submit,   tz = "America/New_York"),
             peak_rss_in_GB = parse_mem_to_gb(peak_rss)) |>
      filter(process == "assign_grnas")
    start_time <- min(curr_trace$submit)
    end_time   <- max(curr_trace$submit + (curr_trace$complete - curr_trace$start))
    data.frame(
      source_data = names(assoc_datasets)[i], method = "pipeline",
      assign_elapsed_time_in_sec = as.numeric(end_time - start_time, units = "secs"),
      max_peak_rss_in_gb_for_assign_workers = max(curr_trace$peak_rss_in_GB),
      num_assign_analysis_workers = nrow(curr_trace)
    )
  }) |>
    do.call(what = rbind)
}

# Orchestrator: load both modes and stack into one long frame (method, peak_rss_gb,
# realtime_sec, nguides, source_data). The pipeline's guide counts are the full
# (non-downsampled) totals per dataset.
load_comp_bundle <- function() {
  path_to_outputs <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                               "guide_assignment/outputs")
  trace_inmem <- load_trace_files_in_memory_assignment(run_names_real, path_to_outputs)
  assoc_datasets <- list(
    replogle  = "replogle-rd7_comp_ngenes=1000_ntargets=100_gene_thresh=7",
    gasperini = "gasperini_comp_ngenes=1000_ntargets=100_gene_thresh=7"
  )
  trace_pipe <- load_trace_files_pipeline_assignment(assoc_datasets) |>
    left_join(data.frame(source_data = c("replogle", "gasperini"),
                         nguides = c(2666, 13078), stringsAsFactors = FALSE),
              by = "source_data")
  rbind(
    trace_pipe |>
      dplyr::select(method, peak_rss_gb = max_peak_rss_in_gb_for_assign_workers,
                    realtime_sec = assign_elapsed_time_in_sec, nguides, source_data),
    trace_inmem |>
      dplyr::select(method, peak_rss_gb, realtime_sec, nguides, source_data)
  )
}

# Shared prep for the compute plots: filter to a dataset, prettify method, and tag
# each point's execution mode (drives circle vs triangle).
.prep_comp_df <- function(df, source_tag) {
  df |>
    dplyr::filter(source_data == source_tag) |>
    dplyr::mutate(
      method    = factor(method_renamer_comp[method], levels = method_renamer_comp),
      execution = factor(ifelse(method == "pipeline", "Distributed", "In-memory"),
                         levels = c("In-memory", "Distributed"))
    )
}

# Circle = in-memory method, triangle = distributed pipeline.
comp_shape_scale <- function() {
  scale_shape_manual(values = c(`In-memory` = 16, `Distributed` = 17),
                     name = "Execution")
}

# Fig 4a/4c: runtime vs guide count, one dataset per figure.
plot_comp_runtime <- function(df, title, source_tag) {
  fmt_hms <- function(x) {
    x <- round(x)
    sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60)
  }
  all_breaks <- c(30, 60, 120, 300, 600, 1800, 3600, 7200, 14400, 28800)
  d <- .prep_comp_df(df, source_tag)
  current_max <- max(d$realtime_sec, na.rm = TRUE)
  ggplot(d, aes(nguides, realtime_sec, colour = method, group = method)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 3) +
    scale_y_log10(breaks = all_breaks[all_breaks <= current_max * 1.2], labels = fmt_hms) +
    scale_x_log10() +
    scale_colour_manual(values = method_pal, name = "Method") +
    comp_shape_scale() +
    labs(x = "Number of guides (log scale)",
         y = "Runtime (H:M:S, log scale)", title = title) +
    theme_thesis() +
    guides(colour = guide_legend(nrow = 2))   # wrap Method legend so it fits
}

# Fig 4b/4d: peak memory vs guide count, one dataset per figure.
plot_comp_memory <- function(df, title, source_tag) {
  d <- .prep_comp_df(df, source_tag)
  ggplot(d, aes(nguides, peak_rss_gb, colour = method, group = method)) +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 3) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                  labels = scales::label_number(suffix = " GB")) +
    scale_x_log10() +
    scale_colour_manual(values = method_pal, name = "Method") +
    comp_shape_scale() +
    labs(x = "Number of guides (log scale)",
         y = "Peak memory (log scale)", title = title) +
    theme_thesis() +
    guides(colour = guide_legend(nrow = 2))   # wrap Method legend so it fits
}

# Combined 2x2 of §4: rows = dataset (Replogle / Gasperini), columns = metric
# (runtime / memory). Reuses the single-panel functions so styling stays in sync;
# the identical Method + Execution legends collapse to one shared legend at the
# bottom. Concise per-panel titles ("Dataset - Metric"); the x-axis title is kept
# only on the bottom row since all panels share the "number of guides" x.
plot_comp_grid <- function(df, title = NULL) {
  drop_x <- theme(axis.title.x = element_blank())
  p_rt_repl <- plot_comp_runtime(df, "Replogle — Runtime",  "replogle")  + drop_x
  p_mm_repl <- plot_comp_memory (df, "Replogle — Memory",   "replogle")  + drop_x
  p_rt_gasp <- plot_comp_runtime(df, "Gasperini — Runtime", "gasperini")
  p_mm_gasp <- plot_comp_memory (df, "Gasperini — Memory",  "gasperini")

  grid <- patchwork::wrap_plots(
    list(p_rt_repl, p_mm_repl, p_rt_gasp, p_mm_gasp),
    ncol = 2, guides = "collect")   # row-major: dataset rows, metric cols
  if (!is.null(title))
    grid <- grid + patchwork::plot_annotation(
      title = title, theme = theme_thesis() + theme(plot.title = element_text(size = 14)))
  grid & theme(legend.position = "bottom")
}

# Alternative to plot_comp_grid: dataset is a facet (columns) and metric is the
# vertical stack. A single facet_grid(metric ~ dataset) can't give runtime H:M:S
# and memory GB tick labels at once (one y-scale = one label format), so instead
# we build two facet_grid(~ dataset) plots -- runtime and memory -- each with its
# OWN y-scale reusing the exact breaks/labels of plot_comp_runtime/plot_comp_memory,
# then stack them with patchwork. Result: dataset column strips on top, one shared
# legend, and each metric row keeps its real y-axis title and the identical ticks
# (H:M:S, " GB") of the unfaceted figures. patchwork aligns the dataset columns.
plot_comp_grid_faceted <- function(df, title = NULL) {
  d <- dplyr::bind_rows(
    .prep_comp_df(df, "replogle")  |> dplyr::mutate(dataset = "Replogle"),
    .prep_comp_df(df, "gasperini") |> dplyr::mutate(dataset = "Gasperini")
  ) |>
    dplyr::mutate(dataset = factor(dataset, levels = c("Replogle", "Gasperini")))

  fmt_hms <- function(x) {
    x <- round(x)
    sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60)
  }
  all_breaks <- c(30, 60, 120, 300, 600, 1800, 3600, 7200, 14400, 28800)
  rt_max <- max(d$realtime_sec, na.rm = TRUE)

  base_layers <- function(p) p +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = execution), size = 2.6) +
    facet_grid(~ dataset, scales = "free_x") +   # dataset columns; shared y within each metric
    scale_x_log10() +
    scale_colour_manual(values = method_pal, name = "Method") +
    comp_shape_scale() +
    theme_thesis() +
    guides(colour = guide_legend(nrow = 2))   # wrap Method legend so it fits

  p_rt <- base_layers(ggplot(d, aes(nguides, realtime_sec, colour = method, group = method))) +
    scale_y_log10(breaks = all_breaks[all_breaks <= rt_max * 1.2], labels = fmt_hms) +
    labs(x = NULL, y = "Runtime (H:M:S, log scale)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  # x shown on bottom plot only

  p_mm <- base_layers(ggplot(d, aes(nguides, peak_rss_gb, colour = method, group = method))) +
    scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                  labels = scales::label_number(suffix = " GB")) +
    labs(x = "Number of guides (log scale)", y = "Peak memory (log scale)") +
    theme(strip.text.x = element_blank(),           # dataset strips on the top plot only
          strip.background.x = element_blank())

  grid <- patchwork::wrap_plots(list(p_rt, p_mm), ncol = 1, guides = "collect")
  if (!is.null(title))
    grid <- grid + patchwork::plot_annotation(
      title = title, theme = theme_thesis() + theme(plot.title = element_text(size = 14)))
  grid & theme(legend.position = "bottom")
}

# =============================================================================
# §3 Poisson bug-fix simulations: precision / recall / jaccard vs ground truth
# =============================================================================
# Adapted from benchmarking/guide-assignment/sceptre-nb/bugfix-sims-metrics.R.
# Two assignment methods (buggy vs fixed Poisson variance) are scored against the
# simulation ground truth on three simulation datasets. Guide names encode the
# noise tail (part 2, "NP*") and the perturbation level (part 3, "P*"); each
# (dataset, noise-tail) pair is one simulation "row", and mu_pert differs per
# dataset by orders of magnitude, so each row keeps its own log2(mu) x-axis.

# The simulation datasets, each with its P-label -> mu_pert map (copied from
# additive-model.Rmd / bugfix-sims-metrics.R). The dispersion sweep
# (sims_sum_1np_3p_disp, "Extreme overdisp.") is intentionally excluded.
sim_bugfix_datasets <- list(
  sims_sum_2np_3p = tibble::tibble(
    p = c("P1", "P2", "P3"), mu_pert = c(500, 750, 1000)
  ),
  sims_sum_repeat_old = tibble::tibble(
    p = c("Psmall", "Pmed", "Plarge"), mu_pert = c(970 / 8, 970 / 4, 970)
  )
)

# Method labels: poisbug is the SCEPTRE assignment (uses the chapter-wide SCEPTRE
# red); poisfix is the fixed-variance Poisson comparison. threshglmpois1000_* are
# the thresholded (>=1000) variants. NOTE: the "(thresh. 1000)" names are
# placeholders — rename here if you call them something else.
sim_method_renamer <- c(
  glmpois_poisbug           = "SCEPTRE",
  glmpois_poisfix           = "Fixed Poisson",
  threshglmpois1000_poisbug = "SCEPTRE (thresh. 1000)",
  threshglmpois1000_poisfix = "Fixed Poisson (thresh. 1000)"
)
# Warm pair = SCEPTRE variants, cool pair = fixed-Poisson variants. Names/order
# here set the legend order and the color + shape mapping.
sim_method_pal <- c(
  SCEPTRE                        = unname(method_pal[["SCEPTRE"]]),  # red
  `SCEPTRE (thresh. 1000)`       = "#E69F00",                        # orange
  `Fixed Poisson`                = "#6C96D6",                        # soft blue (dataset-softness)
  `Fixed Poisson (thresh. 1000)` = "#56B4E9"                         # light blue
)
sim_method_shapes <- c(
  SCEPTRE                        = 16,   # circle
  `SCEPTRE (thresh. 1000)`       = 17,   # triangle
  `Fixed Poisson`                = 15,   # square
  `Fixed Poisson (thresh. 1000)` = 18    # diamond
)

# Row labels, one per (dataset, noise-tail) simulation. Newlines are intentional
# (two-line facet strips read better). Keyed by "<dataset>|<np>". Definition order
# sets row order. "burst" = endogenous burst component (NPhighvar bigger than
# NPlowvar; NPpois has none); "overdisp." = perturbation dispersion theta_pert.
sim_row_renamer <- c(
  "sims_sum_2np_3p|NPhighvar"  = "Medium overdisp.\nlarge burst",
  "sims_sum_2np_3p|NPlowvar"   = "Medium overdisp.\nsmall burst",
  "sims_sum_repeat_old|NPpois" = "Large overdisp.\nno burst"
)

# --- metric helpers (inlined from nb-bench_v3-helpers.R) ----------------------
.sim_set_metrics <- function(pred, truth) {
  pred <- unique(pred); truth <- unique(truth)
  n_inter <- length(intersect(pred, truth))
  n_pred  <- length(pred); n_truth <- length(truth)
  n_union <- length(union(pred, truth))
  tibble::tibble(
    precision = if (n_pred  == 0) NA_real_ else n_inter / n_pred,
    recall    = if (n_truth == 0) NA_real_ else n_inter / n_truth,
    jaccard   = if (n_union == 0) NA_real_ else n_inter / n_union
  )
}
.sim_part <- function(x, k)
  vapply(strsplit(x, "_", fixed = TRUE),
         function(z) if (length(z) >= k) z[[k]] else NA_character_, character(1))
.sim_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
.sim_sd   <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)

# Load one simulation dataset's run: read each assignment_matrix_*_<dataset>.rds,
# turn each guide's assigned cells into integer indices, and build ground truth
# (perturbed AND observed nonzero) from true_pert_matrix.rds.
.sim_load_run <- function(dataset_name) {
  ga_dir      <- file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment")
  outputs_dir <- file.path(ga_dir, "outputs", "run_bugfix_sims")
  input_dir   <- file.path(ga_dir, "input_data", dataset_name)

  grna_matrix <- readRDS(file.path(input_dir, "sceptre", "grna_matrix.rds"))
  true_perts  <- readRDS(file.path(input_dir, "true_pert_matrix.rds"))
  guide_names <- rownames(grna_matrix)

  fnames <- dir(outputs_dir)
  assn_fnames <- fnames[grepl("^assignment_matrix", fnames) &
                          grepl(dataset_name, fnames, fixed = TRUE)]

  all_assns <- list()
  for (fname in assn_fnames) {
    method_name <- fname |>
      sub(pattern = "^assignment_matrix_script_", replacement = "") |>
      sub(pattern = paste0("_", dataset_name, ".rds"), replacement = "", fixed = TRUE)
    mat <- readRDS(file.path(outputs_dir, fname))
    all_assns[[method_name]] <- lapply(guide_names,
      function(g) as.integer(which(mat[g, ]))) |> setNames(guide_names)
  }
  true_assns <- lapply(guide_names, function(g)
    as.integer(which(true_perts[g, ] == 1 & grna_matrix[g, ] > 0))) |>
    setNames(guide_names)

  list(grna_matrix = grna_matrix, all_assns = all_assns, true_assns = true_assns)
}

# Summarized per-(method, np, mu_pert, metric) mean +/- CI for one dataset.
.sim_summarize_dataset <- function(dataset_name, p_to_mean) {
  run <- .sim_load_run(dataset_name)
  guide_ids <- rownames(run$grna_matrix)
  guide_meta <- tibble::tibble(
    grna_id = guide_ids,
    np      = .sim_part(guide_ids, 2),
    p       = .sim_part(guide_ids, 3)
  ) |> dplyr::left_join(p_to_mean, by = "p")

  purrr::imap_dfr(run$all_assns, function(method_assns, method_name) {
    purrr::map_dfr(guide_ids, function(g)
      .sim_set_metrics(method_assns[[g]], run$true_assns[[g]]) |>
        dplyr::mutate(grna_id = g, method = method_name, .before = 1))
  }) |>
    dplyr::left_join(guide_meta, by = "grna_id") |>
    tidyr::pivot_longer(c(precision, recall, jaccard),
                        names_to = "metric", values_to = "value") |>
    dplyr::group_by(method, np, p, mu_pert, metric) |>
    dplyr::summarize(mean = .sim_mean(value), sd = .sim_sd(value),
                     n = sum(!is.na(value)), .groups = "drop") |>
    dplyr::mutate(dataset = dataset_name)
}

# Load + score every simulation dataset and return one long tibble ready to plot:
# columns dataset, np, method (relabeled), mu_pert, metric, mean, ymin, ymax,
# sim_label (dataset x noise-tail, the row key). Cache the chunk that calls this.
load_sim_bugfix_metrics <- function(datasets = sim_bugfix_datasets) {
  purrr::imap_dfr(datasets, function(p_to_mean, ds)
    .sim_summarize_dataset(ds, p_to_mean)) |>
    dplyr::mutate(
      ci     = 2 * sd / sqrt(n),
      ymin   = pmax(0, mean - ci),
      ymax   = pmin(1, mean + ci),
      method = dplyr::recode(method, !!!sim_method_renamer),
      sim_label = factor(unname(sim_row_renamer[paste(dataset, np, sep = "|")]),
                         levels = unname(sim_row_renamer))
    )
}

# One row per simulation (dataset x noise-tail); columns = the metrics in
# `metric_order` (subset/reorder to drop e.g. jaccard), lines = the methods in
# `methods` (subset to compare fewer). Each row is its own ggplot (own log2(mu)
# x-axis, since mu_pert differs by orders of magnitude across datasets) stacked
# with patchwork; single method legend.
plot_sim_bugfix_metrics <- function(df, metric_order = c("jaccard", "precision", "recall"),
                                    methods = names(sim_method_pal),
                                    title = "Assignment accuracy vs. ground truth across simulations",
                                    title_size = 14, y_lim_min = 0) {
  metric_labs <- c(jaccard = "Jaccard", precision = "Precision", recall = "Recall")
  if (!is.factor(df$sim_label))
    df$sim_label <- factor(df$sim_label, levels = unique(df$sim_label))
  method_levels <- names(sim_method_pal)[names(sim_method_pal) %in% methods]  # keep palette order
  df <- df |>
    dplyr::filter(metric %in% metric_order, method %in% methods) |>
    dplyr::mutate(
      metric    = factor(metric, levels = metric_order, labels = metric_labs[metric_order]),
      sim_label = droplevels(sim_label),
      method    = factor(method, levels = method_levels)   # legend + palette order
    )
  sims <- levels(df$sim_label)

  row_plot <- function(s, is_top, is_bottom) {
    d       <- dplyr::filter(df, sim_label == s)
    x_breaks <- sort(unique(d$mu_pert))
    # Label as 10^k only when the breaks round to *distinct* powers of 10 (the
    # overdispersion row, where mu ~ 10^5..10^7); otherwise plain comma numbers so
    # rows like 500/750/1000 don't all collapse to "10^3".
    pw <- round(log10(x_breaks))
    x_labels <- if (!anyDuplicated(pw))
      function(b) parse(text = sprintf("10^%d", as.integer(round(log10(b)))))
    else scales::label_number(big.mark = ",")
    p <- ggplot(d, aes(mu_pert, mean, colour = method, shape = method, group = method)) +
      geom_line(linewidth = 0.5) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.04, alpha = 0.7) +
      geom_point(size = 1.8) +
      facet_grid(sim_label ~ metric) +
      scale_x_continuous(trans = "log2", breaks = x_breaks, labels = x_labels) +
      scale_y_continuous(limits = c(y_lim_min, 1)) +   # y_lim_min: 0, or raise to zoom in
      scale_colour_manual(values = sim_method_pal, name = "Method") +
      scale_shape_manual(values = sim_method_shapes, name = "Method") +
      labs(x = if (is_bottom) expression(mu[pert]) else NULL, y = NULL) +
      theme_thesis()
    if (!is_top)    p <- p + theme(strip.text.x = element_blank())  # metric labels once
    p
  }

  rows <- purrr::imap(sims, function(s, i)
    row_plot(s, is_top = i == 1, is_bottom = i == length(sims)))

  patchwork::wrap_plots(rows, ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme_thesis() + theme(plot.title = element_text(size = title_size))
    ) &
    theme(legend.position = "bottom")
}

# Transposed layout of plot_sim_bugfix_metrics: simulation (dataset x noise-tail)
# as *columns*, metric as *rows*. Each column is its own ggplot (own log2(mu)
# x-axis) so simulations with wildly different mu_pert ranges keep independent
# scales; columns share the 0-1 y-axis so duplicate y ticks are dropped on the
# non-leftmost columns. Good when only a couple of metrics are shown (fewer,
# wider rows read better than many stacked rows).
# y_floor: lower edge of the displayed (uncompressed) y-window. All metrics here
# sit well above 0, so 0-1 wastes half the panel; we crop to [y_floor, 1] but keep
# an honest broken-axis cue -- a real 0 tick at the bottom and a squiggle marking
# the compressed [0, y_floor] band. "auto" = nearest 0.1 below the smallest lower
# CI; a number overrides; NA/0 disables the break (plain 0-1 axis).
plot_sim_bugfix_metrics_t <- function(df, metric_order = c("jaccard", "precision", "recall"),
                                      methods = names(sim_method_pal),
                                      title = "Assignment accuracy vs. ground truth across simulations",
                                      title_size = 14, y_floor = "auto") {
  metric_labs <- c(jaccard = "Jaccard", precision = "Precision", recall = "Recall")
  if (!is.factor(df$sim_label))
    df$sim_label <- factor(df$sim_label, levels = unique(df$sim_label))
  method_levels <- names(sim_method_pal)[names(sim_method_pal) %in% methods]  # keep palette order
  df <- df |>
    dplyr::filter(metric %in% metric_order, method %in% methods) |>
    dplyr::mutate(
      metric    = factor(metric, levels = metric_order, labels = metric_labs[metric_order]),
      sim_label = droplevels(sim_label),
      method    = factor(method, levels = method_levels)   # legend + palette order
    )
  sims <- levels(df$sim_label)

  # Resolve the break floor and build a transform that compresses [0, y_lo] into a
  # thin band (`gap` tall in axis units) while leaving [y_lo, 1] unchanged.
  y_lo <- if (identical(y_floor, "auto"))
    floor(min(df$ymin, na.rm = TRUE) * 10) / 10 else y_floor
  do_break <- length(y_lo) == 1 && is.finite(y_lo) && y_lo > 0.1
  gap <- 0.05
  if (do_break) {
    fwd <- function(y) ifelse(y >= y_lo, y, y_lo - gap + (y / y_lo) * gap)
    inv <- function(z) ifelse(z >= y_lo, z, (z - (y_lo - gap)) / gap * y_lo)
    break_trans <- scales::trans_new("ybreak", fwd, inv, domain = c(0, 1))
    y_breaks    <- c(0, seq(y_lo, 1, by = 0.1))
  }

  col_plot <- function(s, is_left, is_right) {
    d        <- dplyr::filter(df, sim_label == s)
    x_breaks <- sort(unique(d$mu_pert))
    # 10^k labels only when breaks round to distinct powers of 10 (see comment in
    # plot_sim_bugfix_metrics); else plain comma numbers.
    pw <- round(log10(x_breaks))
    x_labels <- if (!anyDuplicated(pw))
      function(b) parse(text = sprintf("10^%d", as.integer(round(log10(b)))))
    else scales::label_number(big.mark = ",")
    p <- ggplot(d, aes(mu_pert, mean, colour = method, shape = method, group = method)) +
      geom_line(linewidth = 0.5) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.04, alpha = 0.7) +
      geom_point(size = 1.8) +
      facet_grid(metric ~ sim_label) +
      scale_x_continuous(trans = "log2", breaks = x_breaks, labels = x_labels) +
      scale_colour_manual(values = sim_method_pal, name = "Method") +
      scale_shape_manual(values = sim_method_shapes, name = "Method") +
      labs(x = expression(mu[pert]), y = NULL) +
      theme_thesis()
    if (do_break) {
      # Left panel edge in data space (x is log2 with the default 5% expansion), so
      # the squiggle lands exactly on the y-axis spine -- x = -Inf can't be used
      # here because log2(-Inf) is NaN and the glyph would be dropped.
      lx     <- log2(range(x_breaks)); x_edge <- 2^(lx[1] - 0.05 * diff(lx))
      p <- p +
        scale_y_continuous(trans = break_trans, breaks = y_breaks, minor_breaks = NULL,
                           limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
        coord_cartesian(clip = "off") +
        # squiggle sits mid-band (fwd(y_lo/2) = y_lo - gap/2); text size is in mm so
        # it stays the same width across columns despite differing log2(mu) x-scales.
        annotate("text", x = x_edge, y = y_lo / 2, label = "≈",
                 angle = 90, size = 3, colour = "grey45")
    } else {
      p <- p + scale_y_continuous(limits = c(0, 1))
    }
    if (!is_right) p <- p + theme(strip.text.y = element_blank())  # metric labels once, on right
    if (!is_left)  p <- p + theme(axis.text.y = element_blank(),    # y ticks once, on left
                                  axis.ticks.y = element_blank())
    p
  }

  cols <- purrr::imap(sims, function(s, i)
    col_plot(s, is_left = i == 1, is_right = i == length(sims)))

  patchwork::wrap_plots(cols, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = title,
      theme = theme_thesis() + theme(plot.title = element_text(size = title_size))
    ) &
    theme(legend.position = "bottom")
}
