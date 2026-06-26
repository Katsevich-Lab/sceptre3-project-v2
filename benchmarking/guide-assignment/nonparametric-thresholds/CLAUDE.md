# gRNA-assignment method investigation — get up to speed

This folder investigates **simple, principled methods for gRNA guide assignment** in
single-cell CRISPR screens, benchmarked against the published state of the art. Read this
first; the full writeup is **`report.qmd` → `report.pdf`**. Persistent notes are in project
memory (`grna-assignment-investigation`, `grna-assignment-data-locations`).

**Project context:** this is a continuation of the **sceptre3 project's ongoing gRNA-assignment
benchmarking** (led by Louis Deutsch; team incl. Tim Barry, Eugene Katsevich). See
**`sceptre3-context/CONTEXT.md`** for the surrounding context — the #sceptre3 Slack discussion,
the lab's methods/datasets/metrics, the long-standing SCEPTRE mixture-model "bug", and the
Bioconductor timeline — and `sceptre3-context/reports/` for the project's current assignment
report files. Our sims here ARE the lab's sims; our baselines (sceptre/crispat/threshGLM) ARE the
lab's methods.

**Why this investigation exists:** reimplementing crispat was floated as a fallback, but the lab
is **not** content with that. **SCEPTRE3 is the lab's next flagship paper**, so the goal is a
*legitimate, principled guide-assignment method of our own* — competitive with or better than the
published state of the art — not a reimplementation. This investigation is the search for that
method (the ambient/contingency-table direction is the leading candidate).

## The bottom line (current understanding)

- **Recommended method: the ambient-proportion test.** Assign a (cell, guide) pair when its
  UMI count significantly exceeds the ambient null `lambda = n_c * pi_g` (cell library size ×
  guide's UMI share) at BH-FDR `q`. One knob (`q`). Code: `ambient_test_assign()` in
  `../guide-assignment-pipeline/bin/script/lib/threshold_methods.R` (the methods library — also
  holds `otsu_threshold_log1p`, `smoothed_valley_threshold`, `assign_by_threshold`). Pipeline
  variant: `../guide-assignment-pipeline/bin/script/ambient_test.R` (method `script_ambient_test`).
- **This method IS geomux's statistical core** (Fisher-hypergeometric + BH), minus geomux's
  adaptive log-odds threshold. **fishash** adds a Simpson's-paradox noise correction we did not
  reproduce by hand.
- **Verdict (rigorous, same-cohort):**
  - Barnyard (real, high-signal): ours mean acc **0.951** > fishash 0.942 > otsu/valley ~0.914 >
    geomux 0.796. Ours wins on direct capture, ties on CROP-seq. (`results/barnyard_exact_corrected.csv`)
  - Hard sim regimes (PR curves): **only fishash reaches precision ≥ 0.90**; ours == geomux on a
    lower-precision/higher-recall frontier. In easy regimes all tie. (`results/PR_curves.png`)
  - Net: simple ambient test is best on clean/high-signal data + recall-favoring metrics;
    fishash's correction wins on precision in hard, low-signal regimes (CD4/iPSC/Gasperini-like).
- Simple per-guide thresholds (Otsu, smoothed valley) excel only when modes are cleanly
  separated; they ignore the per-cell library margin, which is what the ambient/contingency
  methods exploit.

## ⚠️ GOTCHAS (read before re-running)

1. **geomux version matters.** `.venv_geomux_v5/` = geomux **0.5.5 (CORRECT** — adaptive
   log-odds threshold, matches the paper). `.venv_geomux/` = 0.2.10 (OLD, fixed threshold —
   gave wrong barnyard numbers). Use **v5** for barnyard reproduction. The panel/PR geomux runs
   used the old venv but with `lor_threshold=0` (explicit), so version is moot there.
2. **Two scoring metrics in play.** Barnyard uses the Fishash Table-2 metric (per cell: ≥1
   correct-species & 0 wrong-species guide). Sims use per-guide precision/recall/Jaccard vs
   `true_pert_matrix.rds` (perturbed AND observed).
3. **Scripts use relative paths** computed from their own location (`HERE`, and
   `GA = HERE/..` = the `guide-assignment/` dir). Keep canonical scripts at this folder's top
   level or fix the path logic if you move them.
4. **Finalized ambient config: `model="hypergeometric", n_iter=1`** (iteration didn't help).

## Directory map

```
nonparametric-thresholds/
├── CLAUDE.md                     # you are here — the only doc; report.qmd/pdf is the writeup
├── report.qmd  → report.pdf      # THE writeup / single source of truth (quarto render report.qmd --to pdf)
├── scripts/                      # all analysis scripts (.R + .py) — see pipeline order below
├── results/                      # current figures + CSVs (report reads 6 png + 4 csv from here)
│   ├── comprehensive_bench/      # exported subset for the panel/PR (regenerable; gitignored)
│   ├── barnyard_cohort_export/   # our QC'd barnyard cohort export
│   └── _archive/                 # superseded/duplicate outputs (gitignored)
├── external/                     # cloned upstream repos + EXACT barnyard reproduction
│   ├── fishash/ , fishash_analysis/   # authors' repos (gitignored)
│   └── repro_work/               # EXACT Table-2 reproduction + score_ours_table2.R (canonical barnyard)
├── literature/                   # method PDFs to read (cleanser, crispat, fishash, geomux, sceptre)
├── sceptre3-context/             # project context: #sceptre3 Slack (gRNA assignment) + the lab's report files
├── archive/                      # superseded scripts + stale SUMMARY.md — provenance only
└── .venv_geomux_v5/ (correct), .venv_geomux/ (old)   # gitignored
```

**Path conventions (after the reorg):** R scripts in `scripts/` start with a shim
`HERE <- dirname(HERE)` so `HERE` resolves to this folder root and all `results/` and
`../guide-assignment-pipeline/...` paths work — run them from anywhere
(`Rscript scripts/foo.R`). The `.py` geomux runners use **CWD-relative** `results/...` paths,
so **run them from this folder root** (`python3 scripts/run_geomux_*.py`). `report.qmd` reads
`results/...` relative to itself, so render it from this folder root. SUMMARY.md was superseded
by the report and is archived.

## Canonical pipeline (re-run order)

All scripts below live in `scripts/`. Run from this folder root, e.g. `Rscript scripts/characterize_grna.R`.


1. **Characterize real data** → parameter ranges: `characterize_grna.R`
   (reads `~/data/external/perturbseq-survey/*/grna_matrix.rds` + Gasperini/Replogle refs).
2. **Separation axes**: `realism_separation.R`, `survey_separation.R`.
3. **Build sims**: `gasperini_sim_make.R`, `build_comprehensive_sim.R` (→
   `input_data/sims_comprehensive/`), `verify_comprehensive.R`.
4. **Sim benchmark vs baselines**: `comprehensive_sim_eval.R`.
5. **Full panel** (ambient/otsu/valley + real geomux/fishash): `bench_export_and_ours.R` →
   `run_geomux_bench.py` (in `.venv_geomux`) → `score_comprehensive_panel.R`.
6. **PR curves**: `pr_curves.R` → `run_geomux_pr.py` → `combine_pr.R`.
7. **Barnyard (EXACT)**: in `external/repro_work/` — `build_inputs.py`, `run_fishash_local.R`,
   geomux via `.venv_geomux_v5`, `score_table2.R` (theirs) + `score_ours_table2.R` (adds ours).
   Combined table written by an inline step → `results/barnyard_exact_corrected.csv`.
8. **Diagnostics**: `explore_grnas.R` (per-guide log-log histograms with thresholds).
9. **Report**: `report.qmd`.

## Data locations
- Benchmark inputs/outputs: `~/data/projects/sceptre3/benchmarking/guide_assignment/{input_data,outputs}/`
  (HPC pull: `hpcc pull SCEPTRE3 <subdir>`). Comprehensive sim at `input_data/sims_comprehensive/`.
- Barnyard raw (Liu 2025): `~/data/external/liu-2025-cleanser/GSE272457/`.
- 10 recent (2024-26) survey datasets: `~/data/external/perturbseq-survey/`.

## Open threads

### ⚑ KNOWN GAP / CRITICAL-PATH PROBLEM: reproduce fishash's Simpson's-paradox correction
We were **unable to reproduce fishash's Simpson's-paradox correction** from the paper text alone.
This matters because that correction is *exactly* what lets fishash reach high precision in the
hard, low-signal regime (PR curves, `results/PR_curves.png`) — the regime our plain ambient test
(= geomux core) cannot crack. **For a flagship method we must close this gap**: a method that only
matches geomux's core but not fishash's correction is not state-of-the-art where it counts.

What we tried and why it failed: our hand version (`archive/fishash_faithful.R`) implemented the
iterative noise re-estimation (Algorithm 1) with a *marginal* noise profile + plain BH. It behaved
like generic iteration — helped one artificial-hard sim (`repeat_old`) but *hurt* CROP-seq and
Gasperini (the opposite of real fishash). The piece we did not capture is the **block-dependent
multiple-testing correction** (their "GS" method, Appendix B) and the exact noise-conditioned 2×2
odds-ratio.

To address (next session): read **fishash Appendix B** (`literature/fishash.pdf`) AND the actual
**fishash R source** (`external/fishash/`, the `padj_method="GS"` path + the `refit` noise loop),
implement it faithfully, and test whether a simple/transparent version recovers the precision
ceiling on the comprehensive sim's hard groups. This is the key to a legitimate flagship method.

### Other
- **Read the `literature/` PDFs** (cleanser, crispat(+supp), fishash, geomux, sceptre) and
  reconcile each against these findings (place on simple↔contingency↔mixture spectrum). Not done.
- Run remaining baselines (cleanser/sceptre, more of crispat) on the comprehensive sim via
  the Nextflow pipeline for a complete head-to-head.
