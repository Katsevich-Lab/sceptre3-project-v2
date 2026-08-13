#!/usr/bin/env Rscript

# The guidebender simulations parameterize infection by the pre-selection
# Poisson rate lambda, then draw a zero-truncated Poisson for infected cells and
# place `hurdle_prob = 0.1` of recovered cells in a zero-guide component.
# Report the resulting mathematical expectation on the post-selection scale
# while retaining lambda on a secondary axis for paper traceability.
fishash_preselection_moi <- c(0.1, 0.3, 0.5, 1, 2, 3, 5, 10)
fishash_hurdle_prob <- 0.1

fishash_mean_guides_per_recovered_cell <- function(lambda) {
  (1 - fishash_hurdle_prob) * lambda / (-expm1(-lambda))
}

fishash_recovered_guide_labels <- function(lambda) {
  guide_load <- fishash_mean_guides_per_recovered_cell(lambda)
  ifelse(
    guide_load < 1.2,
    formatC(guide_load, format = "f", digits = 2),
    formatC(guide_load, format = "f", digits = 1)
  )
}

scale_x_fishash_guide_load <- function(..., secondary = TRUE) {
  ggplot2::scale_x_log10(
    breaks = fishash_preselection_moi,
    labels = fishash_recovered_guide_labels(fishash_preselection_moi),
    sec.axis = if (secondary) {
      ggplot2::dup_axis(
        name = expression("Pre-selection MOI " * lambda * " (paper parameter)"),
        breaks = fishash_preselection_moi,
        labels = format(fishash_preselection_moi, trim = TRUE)
      )
    } else {
      ggplot2::waiver()
    },
    ...
  )
}
