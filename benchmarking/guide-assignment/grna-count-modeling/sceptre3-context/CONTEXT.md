# How this investigation fits into the sceptre3 project

This `nonparametric-thresholds/` investigation is **not standalone** — it is a continuation of
the **sceptre3 project's ongoing gRNA guide-assignment benchmarking**, which has been running in
the lab (Katsevich Lab) and discussed in the **#sceptre3 Slack channel** (C05MA6XN5LZ) from
2023-08 through 2026-06. This folder captures that surrounding context so a future session has it.

Sources captured here:
- `sceptre3_slack_grna_assignment.md` — 596 assignment-relevant Slack messages (of 2542),
  2023-08 → 2026-06. `sceptre3_slack_raw.json` (full dump), `sceptre3_slack_links.txt` (130 links).
- `reports/` — the project's current gRNA-assignment report files (source + compiled), copied
  from the repo (the files Slack points to). See list below.

## The project's benchmarking framing (manuscript/benchmarking.md, competitors.md)

SCEPTRE3 benchmarks competitor methods across **data import → guide assignment → DE testing**.
For **guide assignment** the competitors are **CLEANSER** (ambient-vs-native mixture) and
**crispat** (11 methods incl. GLM/mixture/threshold; the lab uses its Poisson-Gaussian as
representative). **pertpy is excluded** — its guide assignment just reimplements crispat's
Poisson-Gaussian, so there is nothing to compare separately. Datasets: **simulated + Gasperini (~200K) +
Replogle RPE1/K562 essential-wide (~600K) + Replogle K562 genome-wide (~1.2M)**. Evaluation:
accuracy on simulated, crispat's metrics on real, RAM/runtime. (These are exactly the methods,
datasets, and metrics this investigation works with.)

## Current state of the guide-assignment effort (from recent Slack, June 2026)

**Louis Deutsch** has led the guide-assignment work (for his dissertation). Key findings he reported:
- The decision is essentially driven by the count `y_i` alone: `y_i >> 20` is almost surely
  perturbed, `y_i << 20` almost surely not; **the only hard regime is `y ~ 10`**.
- Covariates/offsets are *subtle and weak*: `grna_n_umis`/`grna_n_nonzero` have a **circularity**
  (they include the cell's own `y_i`); GLM fits are dominated by the heavy tail; on Replogle
  `grna_n_nonzero` can even get a *negative* coefficient. Crispat ignores covariates; that helps.
- **The SCEPTRE `method="mixture"` has a long-standing bug**: it computes `exp(o_i)+gamma` instead
  of `exp(o_i+gamma)` (an *additive* shift). It accidentally works well by focusing on the
  decision boundary; it has been hard to beat. Louis's method families: offsets
  `{glmpois (=sceptre), glmpoisgrnafix, threshglmpois1000}` × mixtures `{poisbug (=sceptre bug),
  poisthresh10, pois0nb}`, plus crispat. **His simulation settings are exactly the sims used here**
  (`sims_sum_2np_3p`, `sims_sum_repeat_old`, `sims_sum_1np_3p_disp`:
  `Pois(small·L) + Ber(endog)·Pois(larger·L) + Ber(pert)·NB(mu·L, theta)`).
- **2026-06-25 meeting:** no fully principled method without weird failure modes / hard tuning had
  been found *yet*. A **pure-R crispat (Poisson-Gaussian) reimplementation** was floated as a
  fallback — but **the lab is NOT satisfied with merely shipping crispat.** SCEPTRE3 is the lab's
  next **flagship paper**, so the bar is a *legitimate, defensible guide-assignment method of our
  own*, not a reimplementation of someone else's. **That is the whole reason this
  `nonparametric-thresholds/` investigation exists** — to find/justify a real method (simple,
  principled, and competitive with or better than the published state of the art).
- **Bioconductor angle (Eugene Katsevich):** submit `sceptre` with the current mixture documented
  honestly as *legacy* behavior (maybe alias `method="legacy_mixture"`); add+validate a principled
  Poisson-Gaussian (crispat-style) option in **July–Aug 2026**, settle the high-MOI story before
  the first Bioc release (~mid–late September target).

## Report files copied into `reports/` (source + compiled)

| File | What it is |
|---|---|
| `assignment-writeup.Rmd/.pdf` | The main guide-assignment writeup (real-data Jaccard/precision/recall comparisons of sceptre/crispat/cleanser; the Table-2-style barnyard + sim analyses). |
| `assignment-analysis.Rmd/.pdf` | Earlier assignment analysis (inter-method agreement on real data). |
| `nb-bench_v3.Rmd/.pdf` | Louis's NB-benchmarking: the mu_pert-sweep plots over his sim settings (the framework this investigation's `comprehensive_sim_eval.R` extends). |
| `additive-model.Rmd/.pdf` | Analysis of the additive (bugged-sceptre) mixture form. |
| `manuscript-benchmarking.md`, `manuscript-competitors.md` | Manuscript sections framing the benchmark. |

## How this investigation extends the above

This folder adds, on top of the project's existing work:
1. **Two simple nonparametric methods** (Otsu on log1p, smoothed valley) and a **principled
   ambient-proportion test** (= geomux's core; `ambient_test_assign()` in the pipeline lib).
2. A **rigorous, exact reproduction** of the Fishash (Kamm 2026) Table-2 geomux/fishash barnyard
   numbers, and a same-cohort head-to-head incl. our methods (`external/repro_work/`).
3. A **characterization of 10 recent (2024-26) datasets** → realistic gRNA parameter ranges, and a
   **comprehensive simulation** spanning them (filling the gap that the lab's Replogle-calibrated
   sims left at the low-separation / high-overdispersion end — exactly Louis's "hard regime").
4. **Precision-recall curves** showing fishash's Simpson's-paradox correction wins on precision in
   the hard regime, while the simple ambient test is best on clean/high-signal data.

See `../CLAUDE.md` (orientation) and `../report.pdf` (full writeup) for this investigation's results.
The natural tie-in to the project decision: the **ambient/contingency-table test (geomux/fishash
family)** is the principled direction the team was looking for; fishash specifically addresses the
hard `y~10` regime Louis identified.

**Critical open problem for the flagship method:** we could **not** reproduce fishash's
Simpson's-paradox correction by hand — and that correction is precisely what wins the hard regime
(our plain test = geomux's core cannot reach its precision there). Owning/improving that correction
(transparently, ideally more simply) is the key to a legitimate flagship method that beats the SOTA.
See the "KNOWN GAP" section of `../CLAUDE.md`.
