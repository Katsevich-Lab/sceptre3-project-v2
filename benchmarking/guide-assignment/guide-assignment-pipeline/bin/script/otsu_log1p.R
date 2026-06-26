# script_otsu_log1p: Otsu's between-class-variance threshold on log(1 + count),
# applied per gRNA. Nonparametric, no covariate adjustment. Assign count >= t.
source(file.path(bin_dir, "script", "lib", "threshold_methods.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  out <- assign_by_threshold(grna_matrix, otsu_threshold_log1p)
  list(assignment_matrix = out$assignment_matrix, extras = out$diagnostics)
}
