
######################################################################
###                     1.  STATISTICAL RESULTS                    ###
######################################################################

# 1a. functions for loading statistical results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load_associations_neg <- function(base_fp, run_name, dataset_name) {
  curr_fp <- file.path(base_fp, run_name)
  scep = read.csv(paste0(curr_fp,  "/association_neg_control_sceptre_", dataset_name, ".csv")) 
  
  if(file.exists(mixscale_fp <- paste0(curr_fp, "/association_neg_control_mixscale_", dataset_name, ".csv"))) {
    mix = read.csv(mixscale_fp) |> dplyr::select(-gene_ID)  # redundant with `gene`
  } else {
    mix = NULL
  }
  
  # FR-Perturb: p-values and logFC (base e, not base 2, so we convert) 
  frpert_pvals_raw <- read.table(paste0(curr_fp, "/frperturb/frperturb_results_", dataset_name, "_pvals.txt"))
  
  # gotta melt 
  frpert <- frpert_pvals_raw %>%
    mutate(gene = rownames(.)) |>
    pivot_longer(cols = !gene, values_to = "p_value_frpert", names_to = "target") |>
    mutate(target = str_replace_all(target, "\\.", "-"))
  
  out = list(scep=scep, mix=mix, frpert = frpert)
  out[!sapply(out, is.null)]
}

combine_and_prepare_results_neg <- function(results_list_neg, disc_pairs, eps = NULL) {
  results_neg <- disc_pairs |>
    transmute(
      target = grna_target, gene = response_id
    ) |>
    left_join(results_list_neg$scep |> dplyr::select(target = grna_target, gene = response_id, n_nonzero_trt, n_nonzero_cntrl, p_value_scep = p_value),
              by = c("target", "gene")
    )  |>
    left_join(results_list_neg$frpert |> dplyr::select(target, gene, p_value_frpert),
              by = c("target", "gene")
    ) |>
    mutate(pair_type = "negative control")
  
  if(is.numeric(eps)) {
    results_neg <- results_neg |>
      mutate(
        p_value_scep = pmax(pmin(p_value_scep, 1), eps),
        p_value_frpert = pmax(pmin(p_value_frpert, 1), eps)
      )
  }
  if("mix" %in% names(results_list_neg)) {
    results_neg <- results_neg |>
      left_join(results_list_neg$mix |> dplyr::select(target = grna_target, gene = response_id, p_value_mix = p_weight),
                by = c("target", "gene")
      )
    if(is.numeric(eps)) {
      results_neg <- results_neg |>
        mutate(
          p_value_mix = pmax(pmin(p_value_mix, 1), eps)
        )
    }
  }
  results_neg
}

load_associations_pos <- function(base_fp, run_name, dataset_name) {
  base_fp = file.path(base_fp, run_name)
  scep = read_csv(paste0(base_fp, "/association_on_target_sceptre_", dataset_name, ".csv"), show_col_types = FALSE) |>
    dplyr::select(target = grna_target, gene = response_id, n_nonzero_trt, n_nonzero_cntrl, p_value_scep = p_value)
  if(file.exists(mixscale_fp <- paste0(base_fp, "/association_on_target_mixscale_", dataset_name, ".csv"))) {
    mix = suppressMessages(read_csv(mixscale_fp, show_col_types = FALSE)) |>
      dplyr::select(
        target = ...1,
        gene = gene_ID, 
        p_value_mix = p_weight
      )
  } else {
    mix = NULL
  }
  
  # FR-Perturb: p-values and logFC (base e, not base 2, so we convert) 
  frpert_pvals_raw <- read.table(paste0(base_fp, "/frperturb/frperturb_results_", dataset_name, "_pvals.txt"))
  frpert <- tibble(
    target = rownames(frpert_pvals_raw),
    gene = target,
    p_value_frpert = sapply(rownames(frpert_pvals_raw), function(target) frpert_pvals_raw[target,target])
  )
  
  out = list(scep=scep, mix=mix, frpert = frpert)
  out[!sapply(out, is.null)]
}

combine_and_prepare_results_pos <- function(results_list_pos, eps = NULL) {
  results_pos <- results_list_pos$scep  |>
    left_join(results_list_pos$frpert, by = c("target", "gene"))
  if("mix" %in% names(results_list_pos)) {
    results_pos <- results_pos |>
      left_join(results_list_pos$mix, by = c("target", "gene"))
  }
  
  if(is.numeric(eps)) {
    results_pos <- results_pos |>
      mutate(
        p_value_scep = pmax(pmin(p_value_scep, 1), eps),
        p_value_frpert = pmax(pmin(p_value_frpert, 1), eps)
      )
    
    if("mix" %in% names(results_list_pos)) {
      results_pos <- results_pos |>
        mutate(
          p_value_mix = pmax(pmin(p_value_mix, 1), eps)
        )
    }
  }
  results_pos
}

# 1b. functions for displaying dataset properties ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summarize_statistical_dataset <- function(dataset_name, dataset_type) {
  
  if(! dataset_type %in% c("neg", "pos")) {
    stop("`dataset_type` must be 'pos' or 'neg'.")
  }
  
  curr_scep_fp = paste0(
    .get_config_path("LOCAL_BENCHMARKING_DIR"),
    "/association/", dataset_type, "-control/input_data/",
    dataset_name, "/sceptre"
  )
  
  
  response_matrix = read_rds(file.path(curr_scep_fp, "response_matrix.rds"))
  cat("Dataset:", dataset_name, "\n")
  cat("  #genes =", nrow(response_matrix), "; #cells =", ncol(response_matrix), "\n")
  
  
  grna_target_df = read_csv(file.path(curr_scep_fp, "grna_target_data_frame.csv"), show_col_types = FALSE)
  targets_table = table(grna_target_df$grna_target)
  cat("  #guides =", nrow(grna_target_df), "; #targets =", length(targets_table), "; avg # guides / target:", mean(targets_table) |> round(2), "\n")
  
  grna_matrix = read_rds(file.path(curr_scep_fp, "grna_matrix.rds"))
  num_cells_per_guide = Matrix::rowSums(grna_matrix)
  num_guides_per_cell = Matrix::colSums(grna_matrix)
  num_cells_per_target = tibble(
    num_cells_per_guide = num_cells_per_guide,
    guide = rownames(grna_matrix)
  ) |>
    left_join(
      grna_target_df, by = join_by(guide == grna_id)
    ) |>
    group_by(grna_target) |>
    summarize(num_cells_per_target = sum(num_cells_per_guide)) |>
    pull(num_cells_per_target)
  
  cat("  #cells per guide: median =", median(num_cells_per_guide), "; min =", min(num_cells_per_guide), "; max =", max(num_cells_per_guide), "\n")
  
  cat("  #cells per target: avg =", median(num_cells_per_target), "; min =", min(num_cells_per_target), "; max =", max(num_cells_per_target), "\n")
  
  cat("  #guides per cell: moi (avg) =", mean(num_guides_per_cell), "; min =", min(num_guides_per_cell), "; max =", max(num_guides_per_cell), "\n")
  
  if(file.exists(disc_path <- file.path(curr_scep_fp, "discovery_pairs.rds"))) {
    curr_disc_pairs = read_rds(disc_path)
    cat("  #discovery pairs:", nrow(curr_disc_pairs), "\n")
  }
}

# 1c. functions for plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

revlog_trans <- function (base = exp(1)) 
{
  trans <- function(x) {
    -log(x, base)
  }
  inv <- function(x) {
    base^(-x)
  }
  scales::trans_new(name = paste("revlog-", base, sep = ""), 
                    transform = trans, inverse = inv, breaks = scales::log_breaks(base = base), 
                    domain = c(1e-250, Inf))
}

#' `plot_pvalues()` creates QQ plots and optionally a volcano plot from a vector of p-values.
#' This function can be used independently of the sceptre pipeline when you have p-values
#' from external analyses.
#'
#' @param p_values A numeric vector of p-values (required)
#' @param log_2_fold_change A numeric vector of log-2 fold changes (optional). If provided,
#'   a volcano plot will be included. Must be the same length as p_values if provided.
#' @param alpha Significance threshold for multiple testing (default 0.05)
#' @param multiple_testing_method Method for p-value adjustment (default "BH" for Benjamini-Hochberg)
#' @param plot_title A string to use as the overall plot title (default NULL, no title)
#' @param x_limits For volcano plot: numeric vector of length 2 giving x-axis limits (default c(-1.5, 1.5))
#' @param point_size Size of points in plots (default 0.55)
#' @param transparency Transparency of points (default 0.8)
#' @param return_indiv_plots If FALSE (default), returns combined cowplot; if TRUE, returns list of individual plots
#' @param verbose If TRUE, prints diagnostic information (default FALSE)
#'
#' @return A cowplot object (if return_indiv_plots=FALSE) or list of ggplot objects (if return_indiv_plots=TRUE)
#' @export
#' @examples
#' # With just p-values (3 plots: 2 QQ plots + text summary)
#' pvals <- runif(1000)
#' plot_pvalues(pvals)
#'
#' # With p-values and fold changes (4 plots: 2 QQ plots + volcano + text summary)
#' pvals <- runif(1000)
#' fc <- rnorm(1000, mean = 0, sd = 0.5)
#' plot_pvalues(pvals, log_2_fold_change = fc, plot_title = "My Analysis Results")
plot_pvalues <- function(p_values,
                         is_significant = NULL,
                         log_2_fold_change = NULL,
                         alpha = 0.1,
                         multiple_testing_method = "BH",
                         plot_title = NULL,
                         x_limits = c(-1.5, 1.5),
                         point_size = 0.55,
                         transparency = 0.8,
                         return_indiv_plots = FALSE,
                         verbose = FALSE,
                         qc_name = NULL) {
  
  # Input validation
  if (!is.numeric(p_values) || any(p_values < 0, na.rm = TRUE) || any(p_values > 1, na.rm = TRUE)) {
    stop("p_values must be a numeric vector with values between 0 and 1")
  }
  
  if (!is.null(log_2_fold_change)) {
    if (!is.numeric(log_2_fold_change)) {
      stop("log_2_fold_change must be a numeric vector")
    }
    if (length(log_2_fold_change) != length(p_values)) {
      stop("log_2_fold_change must be the same length as p_values")
    }
  }
  
  # Remove NAs
  valid_idx <- !is.na(p_values)
  if (!is.null(log_2_fold_change)) {
    valid_idx <- valid_idx & !is.na(log_2_fold_change)
  }
  
  p_values <- p_values[valid_idx]
  if (!is.null(log_2_fold_change)) {
    log_2_fold_change <- log_2_fold_change[valid_idx]
  }
  
  # Create data frame
  if (is.null(log_2_fold_change)) {
    result_df <- data.frame(p_value = p_values)
  } else {
    result_df <- data.frame(
      p_value = p_values,
      log_2_fold_change = log_2_fold_change
    )
  }
  
  # Calculate significance
  if(is.null(is_significant)) {
    # used to be result_df$significant but now just its own thing
    significant <- stats::p.adjust(result_df$p_value, method = multiple_testing_method) < alpha
  } else {
    significant <- is_significant
  }
  
  # Diagnostic output if requested
  if (verbose) {
    cat("=== Diagnostic Information ===\n")
    cat("Number of p-values:", nrow(result_df), "\n")
    cat("Range of raw p-values:",
        min(result_df$p_value, na.rm = TRUE), "to",
        max(result_df$p_value, na.rm = TRUE), "\n")
    
    adjusted_pvals <- stats::p.adjust(result_df$p_value, method = multiple_testing_method)
    cat("Range of adjusted p-values:",
        min(adjusted_pvals, na.rm = TRUE), "to",
        max(adjusted_pvals, na.rm = TRUE), "\n")
    
    cat("Alpha threshold:", alpha, "\n")
    cat("Number significant (raw p <", alpha, "):",
        sum(result_df$p_value < alpha, na.rm = TRUE), "\n")
    cat("Number significant (adjusted p <", alpha, "):",
        sum(significant, na.rm = TRUE), "\n")
    cat("Any NAs in significant column?:", any(is.na(significant)), "\n")
    cat("==============================\n\n")
  }
  
  # Get theme (you'll need this helper from sceptre)
  my_theme <- sceptre:::get_my_theme()
  
  # Create QQ plot (bulk)
  p1 <- ggplot2::ggplot(
    data = result_df,
    mapping = ggplot2::aes(y = p_value)
  ) +
    sceptre:::stat_qq_points(
      ymin = 1e-8, size = point_size,
      col = "dodgerblue3",
      alpha = transparency
    ) +
    sceptre:::stat_qq_band() +
    ggplot2::scale_x_reverse() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (bulk)") +
    my_theme
  
  # Create QQ plot (tail)
  p2 <- ggplot2::ggplot(
    data = result_df,
    mapping = ggplot2::aes(y = p_value)
  ) +
    sceptre:::stat_qq_points(
      ymin = 1e-8, size = point_size,
      col = "dodgerblue3",
      alpha = transparency
    ) +
    sceptre:::stat_qq_band() +
    ggplot2::scale_x_continuous(trans = sceptre:::revlog_trans(10)) +
    ggplot2::scale_y_continuous(trans = sceptre:::revlog_trans(10)) +
    ggplot2::labs(x = "Expected null p-value", y = "Observed p-value") +
    ggplot2::geom_abline(col = "black") +
    ggplot2::ggtitle("QQ plot (tail)") +
    my_theme
  
  # Create volcano plot if fold changes provided
  if (!is.null(log_2_fold_change)) {
    discovery_set <- result_df |>
      dplyr::filter(significant & !is.na(significant))
    p_thresh <- if (nrow(discovery_set) >= 1L) max(discovery_set$p_value) else NA
    
    p3 <- sceptre:::make_volcano_plot(
      discovery_result = result_df,
      p_thresh = p_thresh,
      transparency = transparency,
      point_size = point_size,
      x_limits = x_limits
    )
  }
  
  # Create text summary
  n_rejections <- sum(significant, na.rm = TRUE)
  n_tests <- nrow(result_df)
  str <- paste0(
    "Num. signif\n(at alpha ",
    signif(alpha, 1), "):\n",
    n_rejections, " of ", n_tests
  )
  if(!is.null(qc_name)) {
    str <- paste0("QC: ", qc_name, "\n", str)
  }
  
  p_text <- ggplot2::ggplot() +
    ggplot2::annotate(geom = "text", label = str, x = 1.1, y = 1.2) +
    ggplot2::theme_void() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = "white")) +
    ggplot2::xlim(c(0, 2)) +
    ggplot2::ylim(c(0, 2))
  
  # Combine plots
  if (return_indiv_plots) {
    if (!is.null(log_2_fold_change)) {
      out <- list(p1, p2, p3, p_text)
    } else {
      out <- list(p1, p2, p_text)
    }
  } else {
    if (!is.null(log_2_fold_change)) {
      plot_grid_obj <- cowplot::plot_grid(
        p1, p2, p3, p_text,
        labels = c("a", "b", "c", ""),
        rel_heights = c(0.55, 0.45),
        nrow = 2
      )
    } else {
      plot_grid_obj <- cowplot::plot_grid(
        p1, p2, p_text,
        labels = c("a", "b", ""),
        rel_widths = c(1, 1, 0.7),
        nrow = 1
      )
    }
    
    # Add title if provided
    if (!is.null(plot_title)) {
      title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(
          plot_title,
          fontface = 'bold',
          size = 14,
          x = 0.5,
          hjust = 0.5
        )
      
      out <- cowplot::plot_grid(
        title_grob,
        plot_grid_obj,
        ncol = 1,
        rel_heights = c(0.05, 1)
      )
    } else {
      out <- plot_grid_obj
    }
  }
  
  return(out)
}

make_qq_plots <- function(results, qc_col=NULL) {
  if(is.null(qc_col)) {
    results_ = results
  } else {
    results_ = results[results[,qc_col,drop=TRUE], ]
  }

  pvals = results_ |>
    dplyr::filter(pair_type == "negative control") |>
    dplyr::select(starts_with("p_value"))
  plot_names = gsub("p_value_", "", names(pvals))
  
  lapply(plot_names, function(name) {
    plot_pvalues(dplyr::select(pvals, ends_with(name))[,1,drop=T], qc_name=qc_col)
  }) |>
    set_names(plot_names)
}

make_roc_plot <- function(results, qc_col=NULL, dataset_name=NULL) {
  if(!setequal(results$pair_type, c("positive control", "negative control"))) {
    stop("'pair_type' contains unexpected values.")
  }
  
  results_clean <- results |>
    mutate(
      is_positive = pair_type == "positive control"
    )
  if(!is.null(qc_col)) {
    results_clean <- results_clean %>%
      dplyr::filter(.data[[qc_col]])
  }
  
  p_val_cols = names(results_clean)[grepl("^p_value_", names(results_clean))]
  
  roc_list <- lapply(p_val_cols, function(pval_col) {
    roc(
      response = results_clean$is_positive,
      predictor = -log10(results_clean[,pval_col,drop=TRUE]),
      levels = c(FALSE, TRUE),
      direction = "<"
    )
  })
  
  aucs <- lapply(roc_list, auc)
  
  roc_data = Map(
    function(roc_, auc_, name_) {
      tibble(
        method = name_,
        fpr = 1 - roc_$specificities,
        tpr = roc_$sensitivities,
        auc = as.numeric(auc_)  # Convert to numeric
      )
    }, 
    roc_list, aucs, gsub("p_value_", "", p_val_cols)
  ) |>
    do.call(what=rbind) |>
    mutate(
      method_label = sprintf("%s (AUC = %.3f)", method, auc)
    )
  
  title_ <- "Separation of pos. and neg."
  if(!is.null(dataset_name)) {
    title_ <- paste0(dataset_name, ": ", title_)
  }
  if(!is.null(qc_col)) {
    title_ <- paste0(title_, " (QC =", qc_col, ")")
  }
  
  # Plot ROC curves
  p_roc <- ggplot(roc_data, aes(x = fpr, y = tpr, color = method_label)) +
    geom_line(linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c("#E41A1C", "#4DAF4A", "#377EB8"),
      name = NULL
    ) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
      title = title_,
      x = "False Positive Rate",
      y = "Recall"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = c(0.65, 0.25),
      legend.background = element_rect(fill = "white", color = "gray80"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    coord_fixed()
  p_roc
}




make_trt_zero_plots <- function(results, qc_col=NULL, cols_to_drop=NULL, dataset_name=NULL) {
  results_ = results
  if(!is.null(qc_col)) {
    results_ = results[results[,qc_col,drop=TRUE],]
    title_ = paste0("All pairs with n_nonzero_trt = 0 (QC = ", qc_col, ")")
  } else {
    title_ = paste0("All pairs with n_nonzero_trt = 0 (no QC)")
  }
  if(!is.null(cols_to_drop)) {
    results_ = dplyr::select(results_, -any_of(cols_to_drop))
  }
  if(!is.null(dataset_name)) {
    title_ <- paste0(dataset_name, ": ", title_)
  }
  
  results_ |>
    filter(n_nonzero_trt == 0) |> 
    # mutate(
    #   cntrl_is_pos = n_nonzero_cntrl > 0
    # ) |>
    dplyr::select(starts_with("p_value"), pair_type = pair_type) |>
    pivot_longer(cols = !pair_type, names_to = "method", values_to = "p-value") |>
    ggplot(aes(x = `p-value`, fill = method)) +
    geom_histogram(bins=50) +
    facet_wrap(method ~ pair_type, scales="free_y")+
    theme_bw() +
    labs(y ="count", title = title_)
}



make_fp_budget_plot <- function(
    results,
    qc_col = NULL,
    fp_budgets = c(0, 1, 2, 3, 4, 5, 10, 15, 20, 25, 30),
    y_var = c("tp", "tpr"),
    pseudo_log_sigma_x = 1,
    pseudo_log_sigma_y = 1,
    dataset_name=NULL
) {
  
  y_var <- match.arg(y_var)
  
  pair_col   <- "pair_type"
  pos_label  <- "positive control"
  neg_label  <- "negative control"
  
  tp_at_fp <- function(df_m) {
    n_pos <- sum(df_m$pair_type == pos_label)
    n_neg <- sum(df_m$pair_type == neg_label)
    
    steps <- df_m %>%
      dplyr::mutate(
        is_pos = pair_type == pos_label,
        is_neg = pair_type == neg_label
      ) %>%
      dplyr::group_by(p_value) %>%
      dplyr::summarise(
        tp_step = sum(is_pos),
        fp_step = sum(is_neg),
        .groups = "drop"
      ) %>%
      dplyr::arrange(p_value) %>%
      dplyr::mutate(
        tp_cum = cumsum(tp_step),
        fp_cum = cumsum(fp_step)
      )
    
    idx <- findInterval(fp_budgets, steps$fp_cum)
    
    fp_actual <- integer(length(fp_budgets))
    tp        <- integer(length(fp_budgets))
    threshold <- rep(NA_real_, length(fp_budgets))
    
    ok <- idx > 0
    fp_actual[ok] <- steps$fp_cum[idx[ok]]
    tp[ok]        <- steps$tp_cum[idx[ok]]
    threshold[ok] <- steps$p_value[idx[ok]]
    
    big <- fp_budgets >= n_neg
    fp_actual[big] <- n_neg
    tp[big]        <- n_pos
    threshold[big] <- 1
    
    tibble::tibble(
      fp_budget = fp_budgets,
      fp_actual = fp_actual,
      threshold = threshold,
      tp        = tp,
      tpr       = if (n_pos > 0) tp / n_pos else NA_real_,
      n_pos     = n_pos,
      n_neg     = n_neg
    )
  }
  
  if (!is.null(qc_col)) {
    results_in <- results %>%
      dplyr::filter(.data[[qc_col]])
  } else {
    results_in <- results
  }
  
  p_val_cols <- names(results_in)[grepl("^p_value_", names(results_in))]
  p_cols_present <- gsub("^p_value_", "", p_val_cols)
  names(p_cols_present) <- p_val_cols
  
  stopifnot(length(p_cols_present) > 0, pair_col %in% names(results_in))
  
  long <- results_in %>%
    dplyr::select(dplyr::all_of(c(pair_col, names(p_cols_present)))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(p_cols_present)),
      names_to = "method_col",
      values_to = "p_value"
    ) %>%
    dplyr::mutate(
      method    = dplyr::recode(method_col, !!!as.list(p_cols_present)),
      pair_type = .data[[pair_col]],
      p_value   = as.numeric(p_value)
    ) %>%
    dplyr::filter(pair_type %in% c(pos_label, neg_label), !is.na(p_value)) %>%
    dplyr::select(method, pair_type, p_value)
  
  tp_fp_budget_tbl <- long %>%
    dplyr::group_by(method) %>%
    dplyr::group_modify(~ tp_at_fp(.x)) %>%
    dplyr::ungroup()
  
  plot_df <- tp_fp_budget_tbl %>%
    dplyr::arrange(method, fp_actual, fp_budget) %>%
    dplyr::group_by(method, fp_actual) %>%
    dplyr::summarise(
      tp  = max(tp),
      tpr = max(tpr),
      .groups = "drop"
    )
  
  title_ <- if (y_var == "tp") {
    "#TP recovered vs #FP allowed"
  } else {
    "TPR vs #FP allowed"
  }
  if(!is.null(dataset_name)) {
    title_ <- paste0(dataset_name, ": ", title_)
  }
  
  if (!is.null(qc_col)) {
    title_ <- paste0(title_, " (QC = ", qc_col, ")")
  }
  
  max_fp <- max(plot_df$fp_actual, na.rm = TRUE)
  
  x_breaks <- c(0, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500, 1000,
                2000, 5000, 10000, 20000, 50000, 100000)
  x_breaks <- x_breaks[x_breaks <= max_fp]
  if (!0 %in% x_breaks) x_breaks <- c(0, x_breaks)
  
  if (y_var == "tp") {
    max_y <- max(plot_df$tp, na.rm = TRUE)
    y_breaks <- c(0, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 500,
                  1000, 1500, 2000, 5000, 10000)
    y_breaks <- y_breaks[y_breaks <= max_y]
    if (!0 %in% y_breaks) y_breaks <- c(0, y_breaks)
    y_lab <- "#TP (pos. controls recovered)"
    y_trans <- scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma_y)
  } else {
    y_breaks <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05,
                  0.1, 0.2, 0.5, 1)
    y_lab <- "TPR (recall on positives)"
    y_trans <- scales::pseudo_log_trans(base = 10, sigma = 0.01)
  }
  
  p_tp_fp <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = fp_actual, y = .data[[y_var]], color = method)
  ) +
    ggplot2::geom_step(direction = "hv", linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = pseudo_log_sigma_x),
      breaks = x_breaks,
      labels = scales::label_number(big.mark = "")
    ) +
    ggplot2::scale_y_continuous(
      trans = y_trans,
      breaks = y_breaks,
      labels = if (y_var == "tp") scales::label_number(big.mark = "") else scales::label_percent(accuracy = 1)
    ) +
    ggplot2::labs(
      title = title_,
      x = "#FP (neg. controls called significant)",
      y = y_lab,
      color = "Method"
    ) +
    ggplot2::theme_bw(base_size = 14)
  
  p_tp_fp
}



make_bh_vs_pi_plot <- function(
    results,
    N_total,
    pi_grid,
    qc_col = NULL,
    q_target = 0.1,
    B        = 100,
    p_floor  = 1e-250,
    dataset_name=NULL
) {
  
  set.seed(1)
  pair_col   <- "pair_type"
  pos_label  <- "positive control"
  neg_label  <- "negative control"
  
  p_val_cols = names(results)[grepl("^p_value_", names(results))]
  p_cols_present <- gsub("^p_value_", "", p_val_cols)
  names(p_cols_present) <- p_val_cols
  
  stopifnot(length(p_cols_present) > 0, pair_col %in% names(results))
  
  results_ <- if(is.null(qc_col)) results else filter(results, .data[[qc_col]])
  pools <- results_ %>%
    dplyr::filter(
      pair_type %in% c("positive control", "negative control")
    ) %>%
    select(pair_type, all_of(names(p_cols_present))) %>%
    pivot_longer(all_of(names(p_cols_present)),
                 names_to = "method_col",
                 values_to = "p") %>%
    filter(!is.na(p)) %>%
    mutate(
      method = recode(method_col, !!!as.list(p_cols_present)),
      p      = pmax(as.numeric(p), p_floor),
      is_pos = pair_type == "positive control",
      is_neg = pair_type == "negative control"
    ) %>%
    select(method, p, is_pos, is_neg)
  
  bench <- tidyr::crossing(method = unique(pools$method), pi = pi_grid, b = seq_len(B)) %>%
    group_by(method, pi, b) %>%
    group_modify(~{
      m  <- .y$method
      pi <- .y$pi
      
      pos_pool <- pools %>% filter(method == m, is_pos) %>% pull(p)
      neg_pool <- pools %>% filter(method == m, is_neg) %>% pull(p)
      
      N_pos <- floor(pi * N_total)
      N_neg <- N_total - N_pos
      
      # sample (use replacement only if needed)
      if(N_pos > length(pos_pool)) {
        stop("Trying to sample ", N_pos, "pos. pairs from ", length(pos_pool), " total.\n pi = ", round(pi, 4))
      }
      pos_samp <- sample(pos_pool, size = N_pos, replace = FALSE)#(N_pos > length(pos_pool)))
      neg_samp <- sample(neg_pool, size = N_neg, replace = FALSE)#(N_neg > length(neg_pool)))
      
      p_mix <- c(pos_samp, neg_samp)
      lab   <- c(rep(TRUE, N_pos), rep(FALSE, N_neg))  # TRUE = positive control
      
      p_adj <- p.adjust(p_mix, method = "BH")
      disc  <- p_adj <= q_target
      
      tp <- sum(disc &  lab)
      fp <- sum(disc & !lab)
      D  <- tp + fp
      
      tibble(
        discoveries = D,
        tp = tp,
        fp = fp,
        fdp = ifelse(D == 0, 0, fp / D),
        tpr = ifelse(N_pos == 0, NA_real_, tp / N_pos)
      )
    }) %>%
    ungroup()
  
  
  bench_sum <- bench %>%
    group_by(method, pi) %>%
    summarise(
      fdp_med = mean(fdp, na.rm = TRUE),
      fdp_lo  = quantile(fdp, 0.25, na.rm = TRUE),
      fdp_hi  = quantile(fdp, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(method = fct_inorder(method))
  
  title_ <- paste0("BH q = ", q_target, " (N=", N_total, ", B=", B)
  if(!is.null(qc_col)) {
    title_ <- paste0(title_, ", QC=", qc_col)
  }
  title_ <- paste0(title_, ")")
  if(!is.null(dataset_name)) {
    title_ <- paste0(dataset_name, ": ", title_)
  }
  
  p_fdp_vs_pi_q01 <- ggplot(bench_sum, aes(x = pi, y = fdp_med, color = method, fill = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = fdp_lo, ymax = fdp_hi), alpha = 0.15, color = NA) +
    geom_hline(yintercept = q_target, linetype = "dashed") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,NA)) +
    labs(
      x = expression(pi~"(downsample / prevalence)"),
      y = "FDP",
      color = "Method",
      fill  = "Method",
      title = title_
    ) +
    theme_bw(base_size = 14)
  p_fdp_vs_pi_q01
}


######################################################################
###                   2.  COMPUTATIONAL RESULTS                    ###
######################################################################

load_trace_files_in_memory <- function(run_names, base_path_comp) {
  lapply(run_names, function(run_name) {
    read_tsv(
      file.path(base_path_comp, "outputs", run_name, "trace.txt"),
      show_col_types = FALSE
    )
  }) |>
    do.call(what = rbind)  |>
    tidyr::extract(
      tag,
      into = c("ngenes", "ntargets"),
      regex = "ngenes=([0-9]+)_ntargets=([0-9]+)",
      remove = FALSE,
      convert = TRUE
    ) |>
    transmute(
      method = gsub("_computational", "", process, ignore.case = TRUE) |> tolower(),
      peak_rss, realtime,
      num_cells =  sub(".*ncells=([^_]+).*", "\\1", tag), ngenes, ntargets
    ) |>
    mutate(
      peak_rss_gb = case_when(
        str_detect(peak_rss, regex("\\bGB\\b", ignore_case = TRUE)) ~ parse_number(peak_rss),
        str_detect(peak_rss, regex("\\bMB\\b", ignore_case = TRUE)) ~ parse_number(peak_rss) / 1024,
        TRUE ~ NA_real_
      )
    ) |>
    mutate(
      h = parse_number(str_extract(realtime, "\\d+(?=h)")),
      m = parse_number(str_extract(realtime, "\\d+(?=m)")),
      s = parse_number(str_extract(realtime, "\\d+(?=s)")),
      h = coalesce(h, 0),
      m = coalesce(m, 0),
      s = coalesce(s, 0),
      realtime_sec = 3600*h + 60*m + s,
      realtime_min = realtime_sec / 60
    ) %>%
    select(-h, -m, -s) 
}



load_trace_files_scep_pipe <- function(scep_pipe_datasets, scep_pipe_base_fp) {
  parse_mem_to_gb <- function(x) {
    x <- str_trim(x)
    
    value <- as.numeric(str_extract(x, "^[0-9.]+"))
    unit  <- str_extract(x, "[A-Za-z]+$") |> toupper()
    
    mult <- case_when(
      unit == "B"  ~ 1 / 1024^3,
      unit == "KB" ~ 1 / 1024^2,
      unit == "MB" ~ 1 / 1024,
      unit == "GB" ~ 1,
      unit == "TB" ~ 1024,
      TRUE ~ NA_real_
    )
    
    # hacky wrapper in case any unitless 0's got in there, which become NA
    ifelse(is.na(value * mult), 0, value * mult)
  }
  
  lapply(
    scep_pipe_datasets, function(dataset) {
      curr_path = file.path(
        scep_pipe_base_fp, dataset, "tracing/trace.tsv"
      )
      curr_trace = read_tsv(curr_path, show_col_types = FALSE) |>
        mutate(
          start = as.POSIXct(start, tz = "America/New_York"),
          complete = as.POSIXct(complete, tz = "America/New_York"),
          submit = as.POSIXct(submit, tz = "America/New_York"),
          
          peak_rss_in_GB = parse_mem_to_gb(peak_rss)
        )
      
      # i am treating the time it was submitted as the start, and then adding the
      # actual runtime to this, to avoid counting time spent queued
      # i also am ignoring prepare_association_analysis, and the dummy calib and power checks.
      # these are only ~20sec total, and have the potential to have queue waits
      # and the ~10s from the two dummy checks is equal to the ~10s from prepare_association_analysis,
      # so i think i can call it even.
      
      assoc_data = curr_trace |>
        filter(grepl("run_analysis_subworkflow_.*:run_association_analysis", process))
      
      start_time <- min(assoc_data$submit)
      
      corrected_endtimes = assoc_data$submit + (assoc_data$complete - assoc_data$start)
      end_time = max(corrected_endtimes)
      
      # get max memory for any worker
      peak_rss_in_gb_for_a_worker = curr_trace |>
        filter(grepl("run_analysis_subworkflow_.*:run_association_analysis", process)) |>
        pull(peak_rss_in_GB) |>
        max()
      
      data.frame(
        dataset = dataset,
        assoc_elapsed_time_in_sec = as.numeric(end_time - start_time, units = "secs"),
        max_peak_rss_in_gb_for_assoc_workers = peak_rss_in_gb_for_a_worker,
        num_assoc_analysis_workers = sum(grepl("run_analysis_subworkflow_.*:run_association_analysis", curr_trace$process))
      )
    }
  ) |>
    do.call(what = rbind) |>
    tidyr::extract(
      dataset,
      into = c("ngenes", "ntargets"),
      regex = "ngenes=([0-9]+)_ntargets=([0-9]+)",
      remove = FALSE,
      convert = TRUE
    ) |>
    mutate(
      method = "scep-pipe",
      num_cells =  sub(".*ncells=([^_]+).*", "\\1", dataset)
    ) 
}




plot_runtime <- function(trace_file, title) {
  fmt_hms <- function(x) {
    x <- round(x)
    sprintf("%d:%02d:%02d", x %/% 3600, (x %% 3600) %/% 60, x %% 60)
  }
  
  
  time_breaks <- c(
    60,      # 1 min
    120,     # 2 min
    300,     # 5 min
    600,     # 10 min
    1800,    # 30 min
    3600,    # 1 hr
    7200,    # 2 hr
    14400    # 4 hr
  )
  
  ggplot(
    trace_file,
    aes(
      x = num_cells,
      y = runtime_in_sec,
      color = method,
      group = interaction(method, npairs)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(
      ~ npairs,
      labeller = labeller(npairs = function(x) paste0("num. pairs = ", x, "k"))
    ) +
    scale_y_log10(
      breaks = time_breaks[time_breaks <= max(trace_file$runtime_in_sec)],
      labels = fmt_hms
    ) +
    labs(
      x = "Number of cells",
      y = "Runtime (log scale)",
      color = "Method",
      linetype = NULL,
      title = title
    ) +
    theme_bw()
}

plot_memory <- function(trace_file, title) {
  ggplot(
    trace_file,
    aes(
      x = num_cells,
      y = mem_in_gb,
      color = method,
      group = interaction(method, npairs)
    )
  ) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(
      ~ npairs,
      labeller = labeller(npairs = function(x) paste0("num. pairs = ", x, "k"))
    ) +
    scale_y_log10(
      breaks = c(0.5, 1, 2, 4, 8, 16, 32, 64, 128),
      labels = scales::label_number()
    ) +
    labs(
      x = "Number of cells",
      y = "Peak memory (GB, log2 scale)",
      color = "Method",
      title = title,
      linetype = NULL
    ) +
    theme_bw()
}