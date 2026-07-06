# gRNA-assignment method investigation — get up to speed

This folder investigates **simple, principled methods for gRNA guide assignment** in single-cell
CRISPR screens, benchmarked against the published state of the art.

**Bottom line:** adopt fishash's contingency core + **decontaminate the cell margin** (use the
ambient depth, not raw library size). Method = **`scripts/contingency_method.R`** ("depth_fix", also
called **"fishash+"** in the newer docs — same thing); with `cell_margin="observed"` it reproduces
fishash exactly, with `cell_margin="ambient"` it is our method. See memory
`grna-assignment-recommendation`. Persistent notes: `grna-assignment-investigation`,
`grna-assignment-sim-framework`, `grna-assignment-data-locations`, `cleanser-scoring-artifact`.

> **Folder note:** this folder was renamed from `nonparametric-thresholds/` → `grna-count-modeling/`
> (commit `b84574f`). Scripts run from the folder root and use relative `source("scripts/...")`.
> The shared dataset registry + loaders are in **`scripts/datasets.R`** (see `scripts/README_ambient_ceiling.md`).

## Documents (read in this order)

| Doc | Role |
|---|---|
| **`literature_review.qmd` → `.html`** | **Conceptual synthesis / statistical story.** The methods landscape and the unifying account: assignment = a per-entry test whose null (rank-1 ambient rate) is contaminated by its own signal; the 2×2 table has three background cells, and each method cleans a different subset (fishash cleans the guide side, CLEANSER the cell side, depth_fix both). NEWEST + sharpest framing. |
| **`sceptre3_assignment_report.qmd` → `.html`** | **Our method + real-data evidence + recommendation + spec.** Barnyard ground truth, CLEANSER scoring-bug fix, FDR-control experiments, the depth_fix recommendation, what we rejected. (Its conceptual "framework / failure-modes" sections are an older, rougher version of `literature_review`'s; its "comprehensive-sim" section is superseded by the sim framework below.) |
| **`simulation_framework_report.qmd` → `.html`** | **Simulation apparatus + sim benchmark evidence.** Model A (semi-synthetic) + Model B (mechanistic regimes), realism gate, 11-method panel. Companion working doc: `SIMULATION_FRAMEWORK.md` (design + DE-pushthrough foundation + PENDING re-run steps). |
| **`method_comparison.qmd` → `.html`** | **NEWEST flagship (collaborator-facing).** fishash vs fishash+ vs CLEANSER across all real datasets: per-cell ambient depth *used* (raw library / rank-1 / ≤2), per-guide "number of cells called ambient", assignment agreement (Jaccard), recall. Built by the **ambient-ceiling suite** (`scripts/README_ambient_ceiling.md`, shared lib `scripts/datasets.R`). |
| **ambient-Poisson ladder** (`replogle_ambient_poisson` → `global_ambient_poisson` → `doublet_overdispersion` → `canonical_model` → `grna-count-modeling.qmd`), each `.qmd`→`.html` | The "is the ambient noise Poisson?" evidence chain: single-dataset (RAMAC) → ~760-guide/6-dataset generalization → direct-capture overdispersion is doublet-driven → the parametric generative model (Poisson ambient + gated NB signal) + EM → collaborator roll-up. All live, cross-linked. |
| **`DE_DATASETS_STATUS.md`** | Live inventory for the downstream SCEPTRE-DE pushthrough (the decisive test). 9 DE-ready datasets. |

Older writeups (`report.qmd/pdf`, `methods_review.qmd/html`, `SUMMARY.md`) and the old
comprehensive-sim script pipeline are in **`archive/`** — provenance only.

**Project context:** continuation of the **sceptre3 project's ongoing gRNA-assignment benchmarking**
(led by Louis Deutsch; team incl. Tim Barry, Eugene Katsevich). See **`sceptre3-context/CONTEXT.md`**
(the #sceptre3 Slack discussion, the lab's methods/datasets/metrics, the long-standing SCEPTRE
mixture-model "bug", the Bioconductor timeline) and `sceptre3-context/reports/`. Our sims here ARE the
lab's sims; our baselines (sceptre/crispat/threshGLM) ARE the lab's methods.

**Why this investigation exists:** reimplementing crispat was floated as a fallback, but the lab is
**not** content with that. **SCEPTRE3 is the lab's next flagship paper**, so the goal is a
*legitimate, principled guide-assignment method of our own* — competitive with or better than the
published state of the art. The contingency/ambient-table direction is that method.

## The bottom line (current understanding)

- **Recommended method: depth_fix** (`scripts/contingency_method.R`). Take fishash's noise-conditioned
  hypergeometric test and supply the focal cell's exposure ("draws") from the **denoised ambient
  depth** (the rank-1 noise column fishash already computes) instead of the raw library size
  `N_{:,c}`. It closes the one gap fishash leaves open.
- **The clean statistical story** (see `literature_review`): every method tests whether `N_{g,c}`
  exceeds the ambient-noise rate `a_g · d_c`; that rate must be estimated from a matrix containing
  the signal you are detecting. On the 2×2 table, signal leaks into the three *background* cells:
  - cleaning the two **other-cells** cells removes a Simpson's-paradox **false-positive** inflation
    → this is **fishash's** correction (a *precision* fix, bites in low-SNR/dense regimes);
  - cleaning the **focal cell's own-signal** cell removes co-occurring-guide **self-masking** → this
    is the **cell-margin** fix (a *recall* fix, bites at high MOI). **CLEANSER already does this**
    (its `L` from ≤2-count guides); **fishash does not** — verified from its Eq 7 (draws = raw
    `N_{:,c}`) and from our reproduction (`cell_margin="observed"` ≡ fishash).
  - **depth_fix cleans all three** — the empty corner of the map; fishash + CLEANSER each do one.
- **The bare ambient-proportion test** (`ambient_test_assign()` in
  `../guide-assignment-pipeline/bin/script/lib/threshold_methods.R`) is geomux's core (Fisher-
  hypergeometric + BH) with neither correction — a fast high-signal default, **not** the flagship
  (it over-rejects on the per-guide FDR metric; its barnyard "win" was a forgiving per-cell metric).
- **Verdict (as last recorded; sim numbers provisional pending the per-guide-fit re-run):**
  - Real barnyard (Table-2 metric): depth_fix ties fishash (0.937 vs 0.942, no regression) incl.
    direct capture. (`results/barnyard_exact_corrected.csv`)
  - Sims (per-guide F1, MOI 1/3/5): depth_fix 0.96/0.93/0.89 vs fishash 0.96/0.89/0.82 — recovers
    the high-MOI recall fishash loses to co-occurring guides.
  - The Simpson fix is a pure precision gain (Fig-7 repro: FDR 0.20→0.003, recall unchanged).

## ⚠️ GOTCHAS (read before re-running)

1. **geomux version matters.** `.venv_geomux_v5/` = geomux **0.5.5 (CORRECT** — adaptive log-odds
   threshold, matches the paper). `.venv_geomux/` = 0.2.10 (OLD, fixed threshold — wrong barnyard
   numbers). Use **v5** for barnyard; the panel/PR geomux runs used the old venv with an explicit
   `lor_threshold=0`, so version is moot there.
2. **Two scoring metrics in play.** Barnyard uses the Fishash Table-2 metric (per cell: ≥1
   correct-species & 0 wrong-species guide). Sims use per-guide precision/recall/Jaccard vs the
   integrated-AND-observed truth.
3. **Two path conventions.** Older figure scripts use a `HERE`/`GA = HERE/..` relative-to-own-location
   idiom; the **ambient-ceiling suite runs from the folder root** and uses `source("scripts/...")` +
   `scripts/datasets.R` (`dataset_paths()`) for all data paths. Run those from `grna-count-modeling/`.
4. **Finalized ambient config: `model="hypergeometric", n_iter=1`** (iteration didn't help the bare
   test; depth_fix does iterate the fishash mask schedule — that's separate).

## Directory map

```
grna-count-modeling/
├── CLAUDE.md                          # you are here — orientation index to the docs above
├── literature_review.{qmd,html}       # conceptual synthesis (READ FIRST)
├── sceptre3_assignment_report.{qmd,html}   # our method + real-data evidence + recommendation
├── simulation_framework_report.{qmd,html}  # sim apparatus + sim evidence
├── method_comparison.{qmd,html}       # NEWEST: fishash vs fishash+ vs CLEANSER (ambient-ceiling suite)
├── {replogle_ambient_poisson,global_ambient_poisson,doublet_overdispersion,canonical_model,grna-count-modeling}.{qmd,html}  # ambient-Poisson ladder
├── SIMULATION_FRAMEWORK.md            # sim framework design + DE foundation + PENDING re-run steps
├── DE_DATASETS_STATUS.md              # DE-pushthrough dataset inventory
├── scripts/                           # LIVE scripts. contingency_method.R (the method);
│                                      #   datasets.R (shared registry+loaders); sim_*.R (framework);
│                                      #   ambient_*/dctap_*/depth_*/model_vs_calls/num_ambient_* (ambient-
│                                      #   ceiling suite — see scripts/README_ambient_ceiling.md);
│                                      #   barnyard_*/collab_*/replogle_* figure scripts; simpson_paradox_*
├── results/                           # figures + CSVs the reports read
│   ├── ambient_ceiling/               # method_comparison outputs (fit_cache/, cleanser_{calls,batch}/,
│   │                                  #   writeup/assign/ all gitignored+regenerable)
│   ├── sim_framework/                 # sim-framework outputs (datasets/ + de/ gitignored)
│   ├── benchmark_update/              # depth_fix / build-vs-adopt figures for the assignment report
│   ├── {ambient_intuition,global_ambient_poisson,replogle_ambient_poisson,collaborator_writeup,method_decision}/  # per-writeup outputs
│   ├── barnyard_cohort_export/        # our QC'd barnyard cohort export (gitignored)
│   └── _archive/                      # superseded outputs incl. comprehensive_sim_legacy/ (gitignored)
├── external/repro_work/               # EXACT barnyard Table-2 reproduction (canonical barnyard)
├── literature/                        # method PDFs (cleanser, crispat, fishash, geomux, sceptre)
├── sceptre3-context/                  # #sceptre3 Slack + the lab's report files
├── archive/                           # provenance only — old report.qmd/pdf, methods_review,
│                                      #   SUMMARY.md, old_comprehensive_sim_pipeline/ (dead scripts)
└── .venv_geomux_v5/ (correct), .venv_geomux/ (old), .venv_crispat/   # gitignored
```

## Pipelines (re-run order)

**Simulation framework + DE (current):** the authoritative run order is in **`SIMULATION_FRAMEWORK.md`**
(§Pipeline). In short: `sim_modelA.R` → `sim_characterize.R` → `sim_regimes.R` → `sim_realism_gate.R`
→ `sim_run_methods.R` → geomux/crispat runners → `sim_score_external.R` → `sim_analysis.R` → render
`simulation_framework_report.qmd`. **NOTE:** the per-guide-fit refinement is edited-but-not-re-run
(see SIMULATION_FRAMEWORK.md §PENDING); committed sim scores are the prior 89-dataset version.

**Assignment-report figures + barnyard (current):** figure scripts in `scripts/` feed
`sceptre3_assignment_report.qmd` — `final_barnyard_figure.R`, `comprehensive_depthfix.R`,
`build_vs_adopt.R`, `ambient_survey_{compute,plots}.R`, `ambient_intuition_{compute,plots}.R`,
`simpson_paradox_{repro,plot}.R`. Barnyard EXACT reproduction lives in `external/repro_work/`
(`build_inputs.py`, `run_fishash_local.R`, geomux via `.venv_geomux_v5`, `score_table2.R` +
`score_ours_table2.R` → `results/barnyard_exact_corrected.csv`).

**Method-comparison / ambient-ceiling suite (current):** feeds `method_comparison.qmd`. Run order
(all from the folder root; full detail in `scripts/README_ambient_ceiling.md`):
`ambient_fit_cache.R` (→ gitignored `results/ambient_ceiling/fit_cache/`) → `writeup_compute.R` →
`depth_used.R` (⚠️ **after** writeup_compute — it overwrites `depths_sampled/depth_summary.csv` with the
depth each method *uses*) → `writeup_cleanser_export.R` → `run_cleanser_batch.sh` (CLEANSER via
`.venv_cleanser` → `cleanser_calls/`) → `num_ambient_perguide.R {highmoi,lowmoi}` → render
`method_comparison.qmd`. Shared dataset registry/loaders: `scripts/datasets.R`.

**Old comprehensive-sim pipeline (ARCHIVED):** the pre-framework apparatus
(`characterize_grna.R`, `build_comprehensive_sim.R`, `comprehensive_sim_eval.R`,
`bench_export_and_ours.R`, `score_comprehensive_panel.R`, `pr_curves.R`, `combine_pr.R`, …) is in
`archive/old_comprehensive_sim_pipeline/`. Superseded by the `sim_*.R` framework; kept for provenance.

## Data locations
- Benchmark inputs/outputs: `~/data/projects/sceptre3/benchmarking/guide_assignment/{input_data,outputs}/`
  (HPC pull: `hpcc pull SCEPTRE3 <subdir>`).
- Barnyard raw (Liu 2025): `~/data/external/liu-2025-cleanser/GSE272457/`.
- Recent (2024-26) survey datasets: `~/data/external/perturbseq-survey/`.
- Envs: `.venv_geomux_v5` (geomux 0.5.5), `.venv_crispat` (jdeu1023 fork), R `fishash` 0.99.0,
  `sceptre` 0.99.0.

## Reference: fishash's Simpson's-paradox correction (understood + reproduced)

The plain one-sided Fisher test on full (signal+noise) margins can call a false positive even when
guide g is truly absent from cell c and g,c are noise-independent: signal pooled into the background
cells of the 2×2 flips the odds ratio > 1 (Simpson's paradox; fishash Prop 1, Fig 2). fishash's fix:
replace the *other-cells* counts with **noise-only** counts → adjusted odds ratio R\* (Eq 7); the noise
counts are latent, estimated by **rank-1 Poisson matrix completion** of the assigned-masked matrix
(`impute_masked_counts`, Eq 8–9), iterated with a mask schedule (Alg 1: fresh mask for i≤3, monotone-AND
after). Multiple testing = **Guo–Sarkar 2020 block** correction (`padj_method="GS"`, cells as blocks).
Reproduced (paper Fig 7) by `scripts/simpson_paradox_repro.R` + `simpson_paradox_plot.R` using the
authors' simulator + real package (`fishash(refit=0)` vs `refit=10`) → `results/simpson_paradox_repro.{csv,png}`.
The real package (v0.99.0) feeds every reported fishash number; **reimplementing it is not necessary.**
Full account in `literature_review`.

## Open threads
- **The decisive test: downstream SCEPTRE DE.** All current metrics are upstream (assignment vs truth);
  the arbiter for FDR-vs-recall trade-offs is DE calibration (neg-control type-I error) + power under
  each assignment. Foundation built (`scripts/sim_de_*.R`, `DE_DATASETS_STATUS.md`); cluster-scale run
  pending. See `SIMULATION_FRAMEWORK.md` §Downstream DE.
- **depth_fix vs a correctly-scored CLEANSER, head to head.** The single most informative comparison
  for our contribution, since CLEANSER shares the cell-margin fix — not yet run directly (see
  `cleanser-scoring-artifact`).
- **Per-guide-fit sim re-run** (SIMULATION_FRAMEWORK.md §PENDING) — scripts edited, not yet re-run.
- **Literature review — DONE.** `literature_review.qmd` reconciles cleanser/crispat/fishash/geomux/
  sceptre against our findings and places each on the map.
