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
- `replogle_rd7_pipeline_dir()` / `load_replogle_rd7_de()` — the replogle-rd7 sceptre-pipeline dir, and
  the `(mc, resp, so)` triple (count matrix + response ODM + fitted `sceptre_object`) its DE-backed
  figures share. Requires `ondisc` + `sceptre` attached by the caller.

## Other shared helper libraries

- `scripts/ambient_lib.R` — `fit_rank1_denoised()` (rank-1 Poisson denoised ambient field, masked EM) and
  `detect_gap()` (clean ambient/signal gap detector). Sourced by the ambient-Poisson ladder scripts
  (`ambient_validation.R`, `global_ambient_gap.R`, `weak_integration_prevalence.R`).
- `scripts/barnyard_io.R` — canonical Liu-2025 barnyard loaders: `load_barnyard()` (purity-gated cohort +
  ambient mask), `load_barnyard_basic()` (GEX-QC-only, caller thresholds purity), `.barnyard_qc()`.
  Sourced by the §1 barnyard figure scripts.

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

## method_comparison_2 figures (run order → feeds `method_comparison_2.qmd`)

`method_comparison_2.qmd` is the prose-narrative companion to `method_comparison.qmd`; it has **no code
chunks** — it `![]()`-references four static figures. The generators are the committed `scripts/mc2_*.R`
below (run from the folder root, after the Writeup-pipeline prerequisites above are present: `fit_cache/`,
`writeup/*.csv`, `writeup/assign/*.rds`, CLEANSER calls). All outputs land in `results/ambient_ceiling/`.

1. `mc2_cleanser_cs_calls.sh` — regenerate the **chemistry-correct** CLEANSER calls for the CROP-seq
   datasets (`--cs`, Poisson) into `cleanser_calls_cs/` (gitignored). Needed for Figs 3–4 to match the
   committed versions; currently covers `gastric_organoid` + `a549`. ⚠️ `run_cleanser_batch.sh` hardcodes
   `--dc`; the mc2 scripts **prefer** `cleanser_calls_cs/` over `cleanser_calls/` per dataset, so with only
   `--dc` calls present the CROP-seq CLEANSER ceilings come out inflated (the doc's Notes flag this).
2. `mc2_dataset_landscape.R` → `dataset_landscape.csv` (per-dataset metrics + authoritative capture
   modality — the **shared table** the other three read) + `dataset_landscape.png` (Fig 1) +
   `dataset_landscape_cross.png`. Reads `writeup/{depth_summary,ceilings,assign_counts,jaccard_ff,depths_sampled}.csv`.
3. `mc2_variability_strips.R` → `dataset_variability_strips.png` (Fig 2). Reads `dataset_landscape.csv`
   (modality) + `writeup/{depths_sampled,ceilings}.csv`. Needs `patchwork`.
4. `mc2_guide_hist.R` → `guide_hist_methods.png` (Fig 3) + `method_ceilings_fig.csv`. Reads the raw
   matrices (`dataset_paths()`), the cached `writeup/assign/*.rds`, and the CLEANSER calls.
5. `mc2_concordance.R` → `method_concordance.png` (Fig 4) + `method_concordance_summary.csv`. Same inputs.
6. `quarto render method_comparison_2.qmd`.

Companion tool (not embedded in the doc): `mc2_guide_discord.R <dataset>` reproduces the per-guide
discordant-vs-concordant case-study panels (`a549_discordant_hist.png`, `cd8_discordant_hist.png`, …) used
in the CLEANSER-vs-fishash+ investigation. `dataset_landscape_cross.png` is a byproduct of step 2.

## Exploratory / diagnostic (not embedded in the writeup)

`dctap_poscontrol_compare.R` (the causal recall-vs-knockdown validation on GATA1/HDAC6/MYC positive
controls — the most decisive real-data check), `dctap_discordance.R`, `dctap_fishash_vs_fishashplus_ambient.R`,
`cleanser_vs_fishash_depth.R`, `ambient_depth_hist.R`, `ambient_ceiling_survey.R` (ambient-ceiling
histograms; compute superseded by `ambient_fit_cache.R`), `init_sensitivity.R`, `model_vs_calls.R`,
`replogle_top_ambient_guide_fit.R`. `depth_mix_panel.R` is a hand-sourced tool
(`source(); depth_mix_panel(ds, guide)`).

## Consolidation status

The previously-duplicated helpers have been extracted to the shared libraries above (each refactor
verified behavior-preserving by re-running and byte-diffing the committed outputs it feeds):
- `fit_rank1_denoised()` / `detect_gap()` → `ambient_lib.R`.
- barnyard `load_basic()` / `barnyard_qc()` → `barnyard_io.R` (`load_barnyard_basic()` / `.barnyard_qc()`).
- the replogle-rd7 `Dm`/`Dp` load boilerplate → `datasets.R` (`load_grna_matrix("replogle-rd7")`,
  `load_replogle_rd7_de()`, `replogle_rd7_pipeline_dir()`); adopted across the collab/ENO1/RAMAC scripts.

Documented exceptions — filesystem-scanning drivers that build their own dataset list by walking
`~/data/external/perturbseq-survey` (so a fixed registry does not cleanly apply); they keep a single
literal for the one non-survey path (replogle-rd7):
- `profile_datasets.R` — landscape profiler; `repl` is used structurally as a directory.
- `run_ambient_validation.R` — scans the survey dir, then appends replogle-rd7 (`REPL`).
Scripts with a *fixed* dataset list whose members are all in the registry (e.g.
`weak_integration_prevalence.R`) do source their paths from `dataset_paths()`.
