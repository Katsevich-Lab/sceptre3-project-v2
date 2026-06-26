# script_smoothed_valley: smoothed discrete PMF on the log(count+1) scale, cut
# at the empirical valley (min smoothed PMF) between the two dominant modes:
#   t = argmin_{m1 < k < m2} p_hat(k).  Assign count >= t. Abstains when the
# distribution is not genuinely bimodal.
source(file.path(bin_dir, "script", "lib", "threshold_methods.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  out <- assign_by_threshold(grna_matrix, smoothed_valley_threshold)
  list(assignment_matrix = out$assignment_matrix, extras = out$diagnostics)
}
