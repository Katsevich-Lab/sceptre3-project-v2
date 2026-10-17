# Fishash+ positive-control comparisons

Retained drivers for the October 8, 2026 analyses. Data and generated results
live outside Git. These scripts use the lab's `.get_config_path()` helper from
the user's R startup profile and paths in `~/.research_config`; do not run with
`Rscript --vanilla`. The recorded analyses used SCEPTRE 0.99.0. R package
dependencies include Matrix, sceptre, ondisc, ggplot2 and sparseMatrixStats
(see each driver's imports).

## DC-TAP K562 high MOI

The shared `datasets.R` loader resolves inputs under
`LOCAL_SCEPTRE3_DATA_DIR/perturbseq-survey/dctap_k562_highmoi_GSE303901/sceptre/`
on both Mac and Betty. The input directory must contain `response_matrix.rds`,
`grna_matrix_aligned.rds` and `grna_target_data_frame.csv`.

From the repository root:

```sh
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/compare_dctap_highmoi.R
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/plot_dctap_comparison.R
```

Results go under
`LOCAL_BENCHMARKING_DIR/association/pos-control/dctap-highmoi-mixture-fishashplus-20261008/`.
The comparison fits mixture and R Poisson Fishash+ assignments, then tests three
promoter unions and 12 individual guides (six promoter and six enhancer).
`fishash_repro_lib.R` supplies the Fishash+ implementation.

## Gasperini

Inputs are under `LOCAL_BENCHMARKING_DIR/guide_assignment/input_data/gasperini/`.
The positive-control comparison reads full-library Fishash+ calls from
`guide_assignment/outputs/gasperini/fishashplus-poisson-20261008/` beneath that
same configured benchmarking root. If those calls are absent,
`run_full_fishashplus.R <dataset-directory> <new-output-directory>` generates
them after validating the RDS, ODM and Matrix Market input representations.

For the mixture arm, `align_mixture_and_plot_overlap.R` aligns the saved calls in
`guide_assignment/outputs/gasperini/sceptre-pipeline/grna_assignment_matrix.rds`
to the current full input. These are saved pipeline assignments, rather than a
freshly fitted mixture model. The script verifies every call against nonzero
input counts and saves `mixture_assignment_matrix_aligned.rds` beside the
Fishash+ calls.

```sh
# Run in a shell with the lab config loaded.
source ~/.research_config
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/align_mixture_and_plot_overlap.R
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/compare_gasperini_positive_controls.R \
  "${LOCAL_BENCHMARKING_DIR}association/pos-control/fishashplus-comparison-20261008" \
  full threshold5,mixture,fishashplus
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/summarize_gasperini_positive_controls.R
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/summarize_three_positive_control_methods.R
Rscript benchmarking/guide-assignment/grna-count-modeling/scripts/plot_positive_control_violins.R mixture
```

The summary and plotting drivers use that same configured output directory.
The Gasperini analysis tests 377 self-targeting promoter-positive-control pairs;
it does not test the globin enhancer controls discussed separately.

## Reruns and interpretation

The comparison drivers reuse existing assignments/results by filename. Use a
separate configured benchmarking root (or the Gasperini output-directory argument)
for a fresh run when changing inputs, package versions or analysis settings.
Summary/plotting scripts write their derived tables and figures in the output
directory. Original October 8 script snapshots in the data directories are
historical copies, not the portable repository drivers.

These are positive-control sensitivity comparisons; neither establishes
negative-control calibration or assignment precision. The guides are not
independent biological replicates, and extreme p-values can hit the numerical
floor or use empirical fallback.
