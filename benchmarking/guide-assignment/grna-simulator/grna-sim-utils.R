#!/usr/bin/env Rscript
# Shared helpers for gRNA simulation work. Simulator-agnostic: nothing here
# depends on which generative model produced the counts.
#
# Sourced by the simulation scripts in this directory and by
# grna-count-modeling/scripts/{sim_validate_hist,sim_poc_paired}.R.


#' Side-by-side UMI histograms, one real guide against one simulated guide.
#'
#' Shared log1p bins across both panels; simulated bars are stacked
#' perturbed/non-perturbed. The y-axis is log10, so bars are drawn from 1 rather
#' than 0 -- bar tops are still the true bin counts.
#'
#' @param umis_real  UMI counts for a real guide, across cells.
#' @param umis_sim   UMI counts for a simulated guide, across cells.
#' @param is_pert    Logical, same length as `umis_sim`: ground-truth perturbation.
#' @param bins       Number of shared bins.
#' @param title      Plot title.
#' @param x_breaks   UMI values to label on the (log1p) x-axis.
plot_umi_histogram_real_vs_sim <- function(
    umis_real,
    umis_sim,
    is_pert,
    bins = 80,
    title = NULL,
    x_breaks = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000,
                 2000, 5000, 10000, 20000, 50000, 100000)
) {
  stopifnot(all(umis_real >= 0, na.rm = TRUE))
  stopifnot(all(umis_sim >= 0, na.rm = TRUE))
  stopifnot(length(umis_sim) == length(is_pert))

  umis_real <- umis_real[is.finite(umis_real)]

  keep_sim <- is.finite(umis_sim) & !is.na(is_pert)
  umis_sim <- umis_sim[keep_sim]
  is_pert <- is_pert[keep_sim]

  x_real <- log1p(umis_real)
  x_sim <- log1p(umis_sim)
  x_all <- c(x_real, x_sim)

  # Shared x-axis bins across real and simulated data.
  h_all <- hist(x_all, breaks = bins, plot = FALSE)
  breaks <- h_all$breaks

  make_bin_df <- function(x, group, panel) {
    bin_id <- cut(
      x,
      breaks = breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    )

    df <- data.frame(
      bin_id = bin_id,
      group = group
    )

    df <- df[!is.na(df$bin_id), , drop = FALSE]

    tab <- as.data.frame(table(df$bin_id, df$group))
    names(tab) <- c("bin_id", "group", "count")

    tab$bin_id <- as.integer(as.character(tab$bin_id))
    tab$count <- as.integer(tab$count)
    tab <- tab[tab$count > 0, , drop = FALSE]

    tab$xmin <- breaks[tab$bin_id]
    tab$xmax <- breaks[tab$bin_id + 1L]
    tab$panel <- panel

    tab
  }

  df_real <- make_bin_df(
    x = x_real,
    group = "Real",
    panel = "Real"
  )

  df_sim <- make_bin_df(
    x = x_sim,
    group = ifelse(is_pert, "Perturbed", "Non-perturbed"),
    panel = "Simulated"
  )

  # Stack simulated bars manually so total bar height is the true bin count.
  # Real panel has only one group.
  df <- rbind(df_real, df_sim)

  df$group <- factor(
    df$group,
    levels = c("Real", "Non-perturbed", "Perturbed")
  )

  df <- df[order(df$panel, df$bin_id, df$group), , drop = FALSE]

  df$ymin <- 0
  df$ymax <- 0

  split_ids <- split(seq_len(nrow(df)), paste(df$panel, df$bin_id, sep = "___"))

  for (idx in split_ids) {
    counts <- df$count[idx]
    cum_counts <- cumsum(counts)

    df$ymin[idx] <- c(0, head(cum_counts, -1))
    df$ymax[idx] <- cum_counts
  }

  # Log y-scale cannot display 0, so bars start visually at 1.
  # The top of each stacked bar is still the true bin count.
  df$ymin_plot <- pmax(df$ymin, 1)
  df$ymax_plot <- df$ymax

  df <- df[df$ymax_plot > df$ymin_plot, , drop = FALSE]

  # Keep only breaks inside the observed range.
  max_count <- max(c(umis_real, umis_sim), na.rm = TRUE)
  x_breaks <- x_breaks[x_breaks <= max_count]
  if (!0 %in% x_breaks) {
    x_breaks <- c(0, x_breaks)
  }

  facet_scales <- if (length(umis_real) == length(umis_sim)) {
    "fixed"
  } else {
    "free_y"
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin_plot,
        ymax = ymax_plot,
        fill = group
      )
    ) +
    ggplot2::facet_wrap(~ panel, nrow = 1, scales = facet_scales) +
    ggplot2::scale_x_continuous(
      breaks = log1p(x_breaks),
      labels = scales::label_number()(x_breaks),
      limits = range(breaks),
      name = "UMI count"
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number(),
      name = "Number of cells"
    ) +
    ggplot2::labs(
      title = title,
      fill = NULL
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 7,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )
}
