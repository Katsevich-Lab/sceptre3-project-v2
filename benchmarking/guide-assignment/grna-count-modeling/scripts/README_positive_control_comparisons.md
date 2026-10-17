# Fishash+ positive-control comparisons

Retained drivers for the October 8, 2026 analyses. Data and generated results
live outside Git. These scripts use the lab's `.get_config_path()` helper from
the user's R startup profile and paths in `~/.research_config`; do not run with
`Rscript --vanilla`. The recorded analyses used SCEPTRE 0.99.0. R package
dependencies include Matrix, sceptre, ondisc, ggplot2 and sparseMatrixStats
(see each driver's imports).

## The 17-dataset survey

For the broader dataset landscape and comparison of Fishash, Fishash+ and
CLEANSER, see [Guide assignment across datasets: how datasets vary, how methods
compare](../method_comparison_2.qmd) ([rendered HTML](../method_comparison_2.html)).
That writeup surveys 17 real gRNA-count datasets; it does not run a positive-control
association benchmark across all 17. The completed comparisons documented below
cover Gasperini and high-MOI DC-TAP K562.

The 10 survey dataset directories are available on Betty under
`LOCAL_SCEPTRE3_DATA_DIR/perturbseq-survey/`, currently
`/vast/projects/ekatsevi/wharton/data/projects/sceptre3/perturbseq-survey/`.
Gasperini and Replogle inputs are in the separate configured benchmarking tree;
barnyard and McCutcheon are also separate from this survey directory. Not every
dataset has a prepared expression matrix or curated positive-control pair list;
see [expression-data status](../DE_DATASETS_STATUS.md) for preparation caveats.

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

## Additional positive-control pairs to investigate

Other datasets in the survey may provide useful perturbation-to-gene positive
controls, including effects on genes other than the perturbation's own target.
The following candidates came from checking published experiments and source
annotations. They are **not additional completed analyses in these drivers**.

| Dataset | Candidate pairs or reference set | Evidence and next step |
|---|---|---|
| Gasperini | Globin LCR HS1–HS4 perturbations versus globin genes, including HBE1, HBG1 and HBG2 | The original [GSE120861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861) `all_deg_results.at_scale` table has a `positive_ctrl` category containing both enhancer and promoter guide groups. Curate the intended enhancer–gene links: not every gene tested against a control perturbation is a positive pair. These are separate from the 377 `selfTSS` controls analyzed above. |
| DC-TAP K562 low MOI | e-GATA1 → GATA1; e-HDAC6 → HDAC6; MYC-e2 → MYC | The local guide annotations contain the same three enhancer controls, with two guides per enhancer, as the high-MOI dataset. Extend the matched analysis to low MOI. The high-MOI individual-guide tests already include these six enhancer guides. |
| McCutcheon CD8 T cells | DNMT1 knockdown → decreased CCR7; FOXO1 knockdown → decreased CCR7 | [McCutcheon et al., Fig. 2 and Supplementary Fig. 7](https://pmc.ncbi.nlm.nih.gov/articles/PMC10703699/) report confirmation in individual-guide experiments. These provide non-self trans-response candidates; distinguish the protein-based validation readout from the RNA endpoint being benchmarked. |
| Zhu CD4 T cells | Published arrayed validation set of nine regulators and 12 interactions involving IL10/IL21; additional trans effects compared with a separate arrayed knockout RNA-seq study | [Zhu et al., Fig. 2 and validation methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC13614707/) report two guides per target, bulk RNA-seq and protein measurements. All 12 RNA effects agreed in direction; 10 were significant at FDR <10%. Extract the exact supported pairs and match stimulation condition and donor/replicate coverage before use. Our prepared expression copy covers one replicate. |

For the McCutcheon study, separate BATF3 overexpression experiments also support
BATF3 → IL7R regulation, but overexpression is a different intervention from our
local CRISPRi sample and should not automatically define its expected response.
A549/MethNet is another regulatory-element screen, but we have not established
an independently validated positive-control subset for it.

Outside the 17-dataset collection, [Alda-Catalinas et al. (2024), primary T-cell
crisprQTL](https://pmc.ncbi.nlm.nih.gov/articles/PMC10835965/), is a possible extension:
three CD2 LCR elements → CD2 and 28 enhancer–gene pairs previously identified in
Gasperini's K562 screen. The latter require a cell-type-transfer caveat (17 of
26 evaluable pairs reproduced in primary CD4 T cells). This dataset is not part
of the survey directory transferred to Betty.

Here, independent validation means separate experimental measurements, which may
be within the same paper; it need not mean a different laboratory. A supported
perturbation–response pair need not represent direct molecular regulation.

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
