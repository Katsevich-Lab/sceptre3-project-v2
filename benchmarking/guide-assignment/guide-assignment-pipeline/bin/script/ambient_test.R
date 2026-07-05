# script_ambient_test: assign guide g to cell c when its count significantly
# exceeds the expected ambient count (cell library size x guide ambient share),
# at BH-FDR level q. Adaptive per-(cell, guide) threshold; one knob q (ambient
# FDR). Contingency-table / hypergeometric family (geomux/fishash); Poisson form.
source(file.path(bin_dir, "script", "lib", "threshold_methods.R"))

assign_grnas_script <- function(response_matrix, grna_matrix, grna_target_df,
                                extra_covariates, formula, moi, cpus) {
  out <- ambient_test_assign(grna_matrix, q = 0.05, model = "hypergeometric")
  list(assignment_matrix = out$assignment_matrix, extras = NULL)
}
