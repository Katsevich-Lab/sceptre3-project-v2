# ambient-ceiling / method-comparison suite

This session added an analysis suite comparing **fishash**, **fishash+** (= `depth_fix`), and **CLEANSER**
on the ambient-noise model each uses, across all real gRNA-count datasets. It culminates in the
collaborator writeup **`method_comparison.qmd` → `.html`**.

All scripts run **from the `grna-count-modeling/` folder root** (e.g. `Rscript scripts/ambient_fit_cache.R`).

## Shared library — `scripts/datasets.R`

Single source of truth for the suite (source it after `scripts/contingency_method.R`):

- `dataset_paths()` — named list of the **17 real datasets'** `grna_matrix.rds` paths.
- `load_grna_matrix(ds_or_path)` — load a guide×cell matrix (by name or path) as a numeric CsparseMatrix.
- `kappa_leq2(gm)` — CLEANSER-style per-cell ambient depth `L_i` (sum of ≤2 counts).
- `fishash_assign(gm)` / `fishashplus_assign(gm)` — `contingency_assign` with the finalized config
  (observed vs ambient cell margin; `tail="hyper", fdr="GS"`).
- `dctap_sceptre_dir("highmoi"|"lowmoi")` — the DC-TAP sceptre dir (aligned matrix + GEX + target df).
- `load_cleanser_assignment(gm, ds)` — read `results/ambient_ceiling/cleanser_calls/<ds>.csv` into an
  assignment matrix aligned to `gm`.

## Writeup pipeline (run order → feeds `method_comparison.qmd`)

1. `ambient_fit_cache.R` → `results/ambient_ceiling/fit_cache/` (fishash+ rank-1 fit params + assignment
   per dataset; **gitignored**, regenerable). Prerequisite for the diagnostics below.
2. `writeup_compute.R` → `results/ambient_ceiling/writeup/*.csv` (fishash+fishash+ assignments, depths,
   ceilings, Jaccard) **then** `depth_used.R` (overwrites `depths_sampled.csv`/`depth_summary.csv` with
   the **depth each method actually uses**: raw library / rank-1 / ≤2). ⚠️ run `depth_used.R` *after*
   `writeup_compute.R`; both need `fit_cache/` present.
3. CLEANSER assignments: `writeup_cleanser_export.R` → `run_cleanser_batch.sh` → `cleanser_calls/<ds>.csv`
   (**gitignored**, regenerable via `.venv_cleanser`; MCMC is slow, so only the ≤253-guide datasets).
4. `num_ambient_perguide.R {highmoi,lowmoi}` → the per-guide `num_ambient_dctap_*.png` embedded in the writeup.
5. `quarto render method_comparison.qmd`.

## Exploratory / diagnostic (not embedded in the writeup)

`dctap_poscontrol_compare.R` (the causal recall-vs-knockdown validation on GATA1/HDAC6/MYC positive
controls — the most decisive real-data check), `dctap_discordance.R`, `dctap_fishash_vs_fishashplus_ambient.R`,
`cleanser_vs_fishash_depth.R`, `ambient_depth_hist.R`, `ambient_ceiling_survey.R` (ambient-ceiling
histograms; compute superseded by `ambient_fit_cache.R`), `init_sensitivity.R`, `model_vs_calls.R`,
`replogle_top_ambient_guide_fit.R`. `depth_mix_panel.R` is a hand-sourced tool
(`source(); depth_mix_panel(ds, guide)`).

## Known remaining tech debt (not yet consolidated)

Duplicated helpers *outside* this suite still await extraction (deferred to avoid silent numeric drift
in the committed writeups they feed — refactor one, re-run, diff before the next):
- **barnyard QC loader**: `barnyard_io.R` is the intended source of truth but `load_basic()`/`barnyard_qc()`
  are re-transcribed inline in `barnyard_ambient_figures.R`, `barnyard_five_cells.R`,
  `barnyard_marginal_fit.R`, `barnyard_benchmark_fishashplus.R`.
- **`fit_rank1_denoised()`**: ~4 near-verbatim copies (`ambient_validation.R` is the cleanest).
- **Replogle Dm/Dp loader**: shared by `collab_dose_response_aggregate.R`, `collab_eno1_violin.R`,
  `collab_eno1_gapfit.R`, `method_comparison_lowermode.R`.
