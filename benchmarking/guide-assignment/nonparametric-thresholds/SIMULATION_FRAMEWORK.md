# Guide-assignment simulation framework

A barnyard-anchored, triangulated simulation framework for benchmarking gRNA
guide-assignment methods. Built on top of the `nonparametric-thresholds/`
investigation. **Canonical writeup: `simulation_framework_report.qmd` → `.html`.**

## The two data-generating models

**Model A — semi-synthetic (the ground-truth anchor).** Real barnyard ambient
(the wrong-species submatrix of a QC'd Liu-2025 sample, provably ambient) +
*injected* synthetic signal. Ambient shape is real by construction; only the
signal is synthetic, so per-(guide,cell) truth is exact. Script: `sim_modelA.R`.
Grid: 4 barnyard samples × purity {0.90 = authors' QC verbatim, 0.99 = strict} ×
3 mu_pert × 3 MOI {low=1, high=5, vhigh=15} = 72 datasets. **Barnyard is owned by
Model A** (doublet-heavy direct-capture is hard to characterize parametrically).

**Model B — mechanistic, as DATA-DERIVED REGIMES (one per real dataset).** An
owned CellBender-style generative model (`build_modelB()` in `sim_modelB.R`:
Dirichlet infection + ambient pools, decoupled cell/droplet depth, rank-1
Poisson/NB ambient, per-guide signal). Rather than hand-pick settings, we
**characterize all real datasets and fit each guide**, then build one regime per
dataset. NOT fishash's package → contingency-family wins are not home-field.

## Model B: characterize → per-guide fit → simulate (the current design)

1. **Characterize** (`sim_characterize.R` → `real_characterization.csv`,
   `real_perguide.rds`): for every real dataset, fit each guide's bimodal UMI
   structure (smoothed-valley). Per guide save: **mu** (mean of above-valley
   "signal" cells), **theta** (NB dispersion, method-of-moments), **signal_frac**
   (fraction above valley), right_mode, separation. Dataset-level: MOI (ambient
   test), ambient "ambiguous middle" (% of below-valley counts ≥3), etc.
2. **Build regimes** (`sim_regimes.R` → `regimes.csv`, `manifest_modelB.csv`,
   datasets `B__regime__<name>`): for each dataset, sample G=150 real guides and
   **simulate each sim guide from a real guide's fit** — its (mu, theta,
   signal_frac). Minimal d/w-sigma so signal ≈ NB(mu, theta) directly. Shared
   rank-1 ambient pool calibrated (kappa_bar, phi_noise, alpha_amb) so the
   simulated below-valley distribution matches the real one.
3. **Validate** (`sim_validate_hist.R` → `regime_histograms*.png`): real-vs-sim
   per-guide UMI histograms, sim guide matched to real guide by signal level.

**Why per-guide fit:** earlier versions used a dataset-level *median* (then a
dataset *pool* draw), so no sim guide matched any real guide and histograms were
off (esp. wide-signal datasets). Fitting (mu, theta, frac) per guide reproduces
the real per-guide histograms — proven in `results/sim_framework/poc_perguide_fit.png`
(dctap mu=365/size=0.9 and replogle mu=1266/size=1.9 both reproduce the wide
signal mode). `build_modelB()` accepts per-guide vectors for `mu_pert`,
`theta_sig`, and `pert_rate` (Bernoulli per-guide perturbation at the real
signal_frac — fixes the many-guide over-density: signal_frac = MOI/n_guides, and
with 150 guides vs real thousands you can match per-cell MOI OR per-guide density,
not both; we match density, MOI/masking lives in Model A).

## Method panel (11 methods)

`sim_lib.R::run_methods()` (R, in-process): `thresh3/10/20` (fixed cutoffs),
`otsu`, `valley` (per-guide thresholds), `ambient` (our ambient test), `fishash`,
`depth_fix` (ours = fishash + ambient-depth cell margin, `contingency_method.R`),
`sceptre` (current mixture, `sceptre_assign()`). Out-of-process: `geomux`
(`.venv_geomux_v5`, `run_geomux_sims.py`), `crispat` Poisson-Gaussian
(`.venv_crispat`, `run_crispat_sims.py`, idempotent/shardable, plotting patched out).
Scoring (`sim_lib.R`): per-guide precision/recall/Jaccard + pooled FDR vs
"integrated AND observed" truth. `realism_stats()` is the method-agnostic gate.

## Pipeline (run order, from folder root)

```
Rscript scripts/sim_modelA.R           # Model A (72 datasets)
Rscript scripts/sim_characterize.R     # characterize 17 real -> real_perguide.rds
Rscript scripts/sim_regimes.R          # Model B = 13 data-derived regimes (barnyard dropped)
Rscript scripts/sim_realism_gate.R     # realism gate vs real + lab sum-process sim
Rscript scripts/sim_run_methods.R      # R panel on all datasets -> scores_R.csv + ext_counts.mtx
.venv_geomux_v5/bin/python scripts/run_geomux_sims.py
for s in 0..5: .venv_crispat/bin/python scripts/run_crispat_sims.py $s 6   # parallel shards
Rscript scripts/sim_score_external.R   # + geomux/crispat -> scores_all.csv
Rscript scripts/sim_analysis.R         # figures + tables
Rscript scripts/sim_validate_hist.R ; scripts/sim_plot_real_distributions.R ; sim_plot_within_regime.R
quarto render simulation_framework_report.qmd --to html
```

Outputs in `results/sim_framework/` (datasets/ gitignored, large). Key figures:
`fig_regime_sweep.png`, `fig_regime_heatmap.png`, `regime_histograms_all.png`,
`real_distributions.png`, `real_characterization.png`, `fig_fdr.png`, `fig_purity.png`.

## Key findings (stable across iterations)

- **depth_fix wins** (mean Jaccard ~0.85, FDR controlled) > crispat > thresh3 >
  fishash > sceptre > ambient > … > geomux (collapses on direct-capture).
- **The cutoff trap**: a fixed cutoff ties the principled methods on easy/low-
  ambient data but its FDR explodes (~0.3–0.5) on high-ambient regimes; the
  verdict only emerges because we span the real hardness (separation 1.9–6.8,
  ambient %≥3 0–34%, MOI 1–39).
- **depth_fix beats fishash more as MOI rises** (Jaccard gain +0.02→+0.12 at
  MOI 1→15) — co-occurrence masking; shown on Model A (real ambient) AND Model B.
- **Corrections to the older report**: the lab sum-process sim's ambient is NOT
  "unrealistically lumpy" (Fano ≈2.8 ≈ real CROP-seq 2.4); its limits are
  structural + coverage. Real ambient differs by chemistry; barnyard direct-
  capture's heavy tail is doublets surviving GEX-purity QC (collapses at 0.99).

## CURRENT STATE (2026-06-26) — IMPORTANT for continuation

- **Last fully-run result** (scores_all.csv, report HTML): the *within-dataset-
  variety* version — Model A (72) + 17 dataset-pool-draw Model B regimes (89
  datasets), depth_fix #1 at 0.849. THIS is what the committed scores/report show.
- **In progress, NOT yet run**: the **per-guide-fit** refinement (this doc's
  design). Scripts `sim_modelB.R`, `sim_characterize.R`, `sim_regimes.R`,
  `sim_validate_hist.R` are EDITED for per-guide fit + barnyard-dropped (13
  regimes), but **the rebuild + re-run + re-render has NOT happened yet** (user
  paused to commit). `real_perguide.rds` on disk is still the OLD (right_mode,
  signal_frac) format until `sim_characterize.R` is re-run.

### PENDING next steps (resume here)
1. `Rscript scripts/sim_characterize.R` (saves new per-guide mu/theta/frac).
2. `rm -rf results/sim_framework/datasets/B__regime__*` then `Rscript scripts/sim_regimes.R` (13 regimes).
3. `Rscript scripts/sim_validate_hist.R` → review `regime_histograms_all.png` (should now match per-guide).
4. Re-run panel: delete stale `B__regime__*/{geomux,crispat}_calls.csv`; `sim_run_methods.R`; geomux; crispat shards; `sim_score_external.R`; `sim_analysis.R`.
5. Update the report's regime count (17→13; barnyard now Model-A-only) and re-render.
6. Deferred (by user): downstream SCEPTRE DE test (calibration/power) — the decisive metric.

## Data / envs
- Barnyard repro inputs: `external/repro_work/*_grna_counts.mtx` (+ `_meta.csv`, `_guides.csv`).
- Survey/Gasperini/Replogle: `~/data/external/perturbseq-survey/`, `~/data/projects/sceptre3/benchmarking/guide_assignment/input_data/`.
- Methods: `.venv_geomux_v5` (geomux 0.5.5, correct), `.venv_crispat` (crispat from jdeu1023 fork per `modules/crispat/environment.yml`), R `fishash` 0.99.0, `sceptre` 0.99.0.
- See `CLAUDE.md` for the older nonparametric-thresholds investigation; memory `grna-assignment-sim-framework`.
