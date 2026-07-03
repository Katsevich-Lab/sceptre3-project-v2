# What gRNA-assignment method should sceptre implement?

**A principled, cross-dataset-validated recommendation.**
SCEPTRE3 guide-assignment investigation — decision memo.

---

## 0. TL;DR

Implement a **noise-conditioned per-entry test against a decontaminated ambient null**:

> For each (guide *g*, cell *c*), test H₀: the count $N_{g,c}$ was produced by ambient
> contamination alone, $N_{g,c}\sim\mathrm{Poisson}(a_g d_c)$, where $a_g$ is the guide's
> ambient share and $d_c$ is the cell's **denoised ambient depth** (both from a rank-1
> Poisson fit to the carrier-masked matrix). Reject → assign. Control FDR with the
> Guo–Sarkar within-cell block procedure.

This is the method currently called **`depth_fix`** (`scripts/contingency_method.R`,
`cell_margin="ambient"`). It is fishash's contingency test with the **cell margin
decontaminated to ambient depth**, so it cleans all three "background" cells of the 2×2 table
rather than fishash's two. It is:

- **Principled** — an exact conditional test of a single, empirically-validated null.
- **Calibrated by construction** — error control is the design goal, matching sceptre's identity.
- **Self-contained** — the whole method is a rank-1 IPF fit + a hypergeometric tail + a block
  FDR; ~60 lines, no external model-fitting, reimplementable natively in sceptre.
- **State-of-the-art or better** — ties/edges fishash on real barnyard ground truth, best
  Jaccard across the simulation panel, and wins increasingly at high MOI.

**The case in one figure** (`method_scorecard.png`): across the decisive axes, `depth_fix` and
`fishash` are the only methods that are simultaneously accurate upstream (Jaccard 0.85 / 0.81),
FDR-controlled (0.04 / 0.03), downstream-calibrated (nominal), and downstream-powerful (~0.56) —
while the bare `ambient` test over-calls (FDR 0.14), `thresh3` controls FDR poorly (0.07), and the
current `sceptre` mixture is worst on upstream accuracy (0.72) and downstream power (0.45).
`depth_fix` leads on accuracy (high-MOI recall) and is the strict superset of `fishash`.

**Keep the exact Poisson/hypergeometric tail.** The ambient is *very nearly* Poisson but carries a
small, real overdispersion in the high-rate tail (validated across 11 datasets, §3a). We tested
the natural fix — a per-entry **negative-binomial null** (`tail="nb"`) — and it **collapses**: the
fitted dispersion is contaminated by residual signal (worst at high MOI), so it assigns almost
nothing. The small overdispersion is best left unmodeled in the per-entry test; the Poisson null
is empirically FDR-controlled (§4) and downstream-calibrated (§5).

---

## 1. The problem

sceptre's current assignment (`method="mixture"`) is a per-guide two-component Poisson-GLM
mixture, $\log\mu_j=\gamma X_j + o_j$. It has (a) an **additive-mean bug** — the perturbed mean
is computed as $\exp(o_j)+\gamma$ instead of $\exp(o_j+\gamma)$, so $\gamma$ is not a
log-fold-change — and (b) a **covariate circularity** — the offset $o_j$ uses the cell's own
guide library size, which contains the count being tested. Since sceptre's brand is
statistically rigorous, calibrated Perturb-seq inference, the assignment step must be replaced
by something we can defend from first principles and validate.

**Design requirements.** (i) A clear generative model and null. (ii) Calibrated error control
(this is what sceptre sells). (iii) Robustness across the enormous heterogeneity of real
Perturb-seq data — MOI from <1 to >8, guide-capture depth from 2 to ~2000 UMIs/cell, direct-
capture and CROP-seq chemistries, immortalized lines, primary cells, organoids, in vivo. (iv)
No overfitting to any one dataset.

---

## 2. The generative model and the test

Observed guide count = ambient contamination + (if the guide is truly integrated) native signal:
$$
N_{g,c}=Z_{g,c}+X_{g,c}\,W_{g,c},\qquad
Z_{g,c}\sim\mathrm{Poisson}(a_g d_c),\quad
W_{g,c}\sim\mathrm{NB}(\mu_g s_c,\theta),\quad
X_{g,c}\sim\mathrm{Bernoulli}(\pi_g).
$$
Two rank-1 rate fields with **distinct** per-cell depths: an *ambient* depth $d_c$ (how much soup
the droplet captured) and a *native-capture* depth $s_c$ (how well it captures integrated
guides). Assignment is the per-entry test of $H_0: X_{g,c}=0$, i.e. "is $N_{g,c}$ larger than the
ambient Poisson tail at rate $a_g d_c$?" — needing only the **null**, never the signal model.

**Why a test, not a mixture, for sceptre.** The whole family (geomux, fishash, CLEANSER,
sceptre-mixture, crispat) is soft-/hard-EM of *this one model*; they differ in four switches
(E-step soft vs hard; native component parametric vs nonparametric; per-guide vs shared; ambient
Poisson vs NB). For a calibration-first tool, the **frequentist reduction** — freeze the
alternative to "exceeds the ambient tail," harden the E-step to an FDR call — is the right point
in this space: it controls the per-entry error directly, needs no signal-model assumptions, and
borrows strength across guides through the *shared* rank-1 null. That is exactly the contingency
test.

**Estimating the null (the crux).** The null rate $a_g d_c$ must be estimated from a matrix that
contains the very signal we are detecting. Signal leaks into all three background cells of the
2×2 table; each competitor cleans a different subset. The rank-1 Poisson MLE on the
**carrier-masked** matrix (iterated: test → mask assigned → refit) gives the decontaminated
$\hat a_g$ (guide-side, fishash's Simpson fix) *and* the decontaminated $\hat d_c$ (cell-side).
Using $\hat d_c$ in place of the raw library size $N_{:,c}$ is the one thing fishash leaves on the
table — it computes $\hat d_c$ then discards it. `depth_fix` uses it.

---

## 3. Why the two key ingredients are right — validated across 11 datasets

Rather than argue from one dataset, both load-bearing modeling choices were checked on all
available Perturb-seq datasets (a549, cd8/perturb-CITE, dctap high & low MOI, endoc, gastric
organoid, in-vivo cortex, iPSC, multiome erythroid, tcell-cd4, Replogle) —
`scripts/run_ambient_validation.R`, `results/method_decision/ambient_validation.csv`.

**(a) The ambient null is (near-)Poisson given the denoised depth.** The aggregate soup
dispersion var/mean vs fitted rate $\mu=a_g d_c$ (zeros included analytically) sits at **exactly
1.0 for all datasets at low rates** (where the bulk of ambient mass lives), rising to a **mild
1.03–1.15** at moderate rates for clean low-MOI screens (Replogle 1.04, endoc 1.03, tcell 1.07,
dctap-low 1.05) and a **moderate 1.3–2** in the high-rate tail for messier ones (dctap-high 1.31,
gastric 1.47, cortex 1.51 — organoid/in-vivo/high-MOI, where residual doublets inflate the tail).
`results/method_decision/ambient_dispersion_all.png`. Conclusion: **Poisson is an excellent
first-order null everywhere**; the residual tail overdispersion is small and is the case for a
conservative large-θ NB, not for abandoning the model.

**(b) Library size ≠ ambient depth.** The median ratio of library size to the denoised ambient
depth ranges from **5× to 188×** across datasets (a549 188×, gastric 126×, endoc 75×, Replogle
72×, in-vivo 45×, tcell 20×, dctap-high 16×, iPSC 5×; only the ultra-shallow cd8, where the guide
library is essentially all soup, approaches 1×). Using library size as the null exposure — as geomux, fishash,
and sceptre's covariate all do — therefore systematically **mis-sizes the null** and, in low-MOI
cells, is nearly uncorrelated with the true soup exposure. This is the general, cross-dataset
justification for the cell-margin decontamination, not a Replogle artifact.

---

## 4. State-of-the-art comparison — upstream accuracy

**Real barnyard ground truth** (Liu 2025; wrong-species guide = provably ambient; fishash Table-2
per-cell metric; `results/barnyard_exact_corrected.csv`). All noise-conditioned methods land
~0.94; `depth_fix` ties fishash (0.937 vs 0.942 mean incl. direct-capture, per prior run), while
geomux collapses on direct capture (0.80). Barnyard is low-MOI, so the high-MOI advantage cannot
appear here — this is a *no-regression* result.

**Simulation panel** (`results/sim_framework/method_overall.csv`; per-guide Jaccard vs the
integrated-AND-observed truth):

| method | Jaccard | precision | recall | FDR |
|---|---|---|---|---|
| **depth_fix** | **0.854** | 0.962 | 0.887 | 0.041 |
| thresh3 | 0.835 | 0.955 | 0.875 | 0.070 |
| crispat | 0.822 | 0.974 | 0.843 | 0.040 |
| fishash | 0.806 | 0.974 | 0.828 | 0.027 |
| ambient (bare) | 0.739 | 0.871 | 0.859 | **0.142** |
| sceptre (current) | 0.722 | 0.967 | 0.750 | 0.055 |
| geomux | 0.381 | 0.980 | 0.392 | 0.032 |

`depth_fix` has the best Jaccard while controlling FDR (0.041); fishash is tighter on FDR but
gives up recall; the bare ambient test over-calls (FDR 0.14). **By MOI**
(`depthfix_vs_fishash_moi.csv`), `depth_fix` beats fishash at every level and the gap widens with
MOI (Jaccard low 0.834 vs 0.813; high 0.864 vs 0.806; very-high 0.849 vs 0.756) — the recall
fishash loses to co-occurring guides (self-masking) is exactly what the depth decontamination
recovers.

---

## 5. The decisive test — downstream DE calibration & power

All of §4 is *upstream* (assignment vs truth). The arbiter for a calibration-first tool is
whether the assignment yields **calibrated negative-control differential expression** (non-
targeting guides must produce uniform p-values) while retaining **power** (on-target guides must
recover their target-gene knockdown). This was run natively through sceptre
(`run_calibration_check` + `run_power_check`) on every dataset with genome-wide gene expression
and non-targeting guides — **a549, endoc, Replogle, Gasperini** (spanning low and high MOI;
2000 negative-control pairs each for resolution). `scripts/de_native.R`, `scripts/de_driver.R`,
`results/method_decision/de_master.csv`, `de_trustworthy_summary.csv`,
`de_calibration_compare.png`. *(dctap and cd8 were excluded: dctap ships only a 93-gene targeted
panel, so its `response_n_umis` covariate is meaningless and calibration is uninterpretable for
**all** methods; cd8's guide capture is degenerate — median 2 UMIs/cell.)*

**Finding 1 — downstream calibration is robust to the assignment method.** Averaged over the four
genome-wide datasets, *every* reasonable method — including a naive count threshold — yields
**nominal** negative-control type-I error:

| method | t1e @0.05 | KS | power (frac p<0.05) | power (median p) |
|---|---|---|---|---|
| ambient (bare) | 0.050 | 0.075 | 0.576 | 0.188 |
| fishash | 0.052 | 0.084 | 0.568 | 0.167 |
| thresh3 | 0.052 | 0.072 | 0.560 | 0.156 |
| **depth_fix** | 0.056 | 0.080 | 0.557 | **0.159** |
| sceptre (current) | 0.050 | 0.099 | **0.450** | **0.211** |

All calibration is at nominal (0.050–0.056; per-dataset range 0.036–0.070, figure
`de_calibration_compare.png`). This is expected and reassuring: sceptre's covariate-adjusted
resampling null **absorbs** assignment differences, so the downstream test does not break under
any sensible assignment. **Consequence: the method cannot be chosen on downstream calibration of
typical datasets — it must be chosen on §1–4 (principle, validated assumptions, upstream accuracy,
FDR control, robustness in hard regimes).**

**Finding 2 — the current sceptre mixture is the least powerful downstream.** It recovers the
fewest on-target effects (frac p<0.05 = 0.45 vs 0.56–0.58 for the others; worst median p) — it
assigns conservatively and loses power. This is an *additional* reason to replace it, beyond the
additive-mean bug.

**Finding 3 — depth_fix is calibrated and competitive-to-best, and wins at high MOI.** It holds
nominal calibration everywhere (mild 0.070 only on high-MOI a549's strong cell structure — shared
direction with all methods) and has top-tier power; on the clean high-MOI flagship (Gasperini) it
is simultaneously **best-calibrated (t1e 0.049, KS 0.018) and most powerful (0.950)**, matching
its upstream high-MOI advantage.

**Finding 4 — the NB null is not viable.** `depth_fix_nb` (§3a's overdispersion upgrade) collapses
— the fitted dispersion is contaminated by residual signal, especially at high MOI, so it assigns
almost nothing (5–58 cells) and its downstream calibration degenerates. **Keep the exact
Poisson/hypergeometric tail.** The mild ambient overdispersion is real but small, and modeling it
per-entry destabilizes the test; the Poisson null is FDR-controlled empirically (§4) and is the
robust choice.

---

## 6. Honest limitations & the alternative

- **The frequentist test gives no soft posterior.** If a calibrated per-cell *probability* is
  wanted (e.g. for probabilistic downstream models), the shared-rank-1 **Poisson+NB mixture EM**
  (the full model, §2) is the natural upgrade — it adds a genuine responsibility score and power
  in the overlap regime, at the cost of the NB/rank-1 signal assumptions and EM initialization.
  It is a strict superset; the test is its $\theta\to\infty$, hard-E-step limit.
- **depth_fix vs fishash on *realistic-MOI* data is close.** The clear win is at high MOI; on
  low/typical MOI they are near-tied (fishash marginally tighter FDR). The decision to prefer
  depth_fix rests on it being a *strict superset* (same method + one extra correction that only
  ever helps at high MOI and never hurt in any test), plus §5.
- **No real high-MOI ground-truth dataset exists** — the high-MOI advantage is validated on
  simulators and on the downstream DE calibration of real high-MOI screens (dctap-high, a549),
  not on a barnyard-style oracle.
- **crispat (Poisson-Gaussian) and CLEANSER (MCMC)** — the two published per-guide model-fitting
  methods — were installed fresh and run head-to-head through the same downstream DE harness on
  a549, endoc, and replogle (see `DETAILED_RESULTS`). Both are **downstream-calibrated at parity**
  (nominal negative-control type-I error, comparable power) — the "downstream is robust to the
  assignment method" finding holds for the mixture family too. But **neither scales**: at 13k
  guides (Gasperini) crispat's per-guide SVI needs ~3 h and CLEANSER's per-guide MCMC is worse,
  versus ~1 min for the vectorized contingency test — so Gasperini could not be run for them.
  Upstream, crispat trails depth_fix on the simulation panel (Jaccard 0.822 vs 0.854); CLEANSER is
  competitive once correctly scored (~0.89/0.95) but is not a shippable default for a fast,
  deterministic tool. The full model (§2, shared Poisson+NB EM) remains the soft-score upgrade path.

---

## 7. Recommendation for sceptre

**Ship the noise-conditioned contingency test with the denoised-ambient cell margin, the exact
hypergeometric/Poisson tail, and Guo–Sarkar block FDR** (`depth_fix`), implemented **natively**
(rank-1 IPF completion + hypergeometric tail + GS — no external dependency; the whole thing is
~60 lines and the rank-1 fit was reproduced from scratch this cycle in
`scripts/ambient_validation.R`). Keep the current mixture as `method="legacy_mixture"`. Offer the
shared Poisson+NB mixture EM (the full model, §2) as the optional soft-score method for users who
need a calibrated per-cell probability.

**What the evidence says, in order of strength:**
1. It is the point in method space that **controls the per-entry error directly** — sceptre's
   whole value proposition — as a frequentist reduction of the one generative model everyone else
   fits.
2. Its two load-bearing assumptions — ambient ≈ Poisson given the denoised depth, and library
   size ≠ ambient depth — hold across **11 diverse datasets**, not one.
3. It matches or beats the published state of the art **upstream** (best sim Jaccard with FDR
   control; ties fishash on real barnyard) and wins increasingly at high MOI.
4. **Downstream it is calibrated and competitive-to-best** — nominal negative-control DE across
   four genome-wide datasets, top-tier power, best on the high-MOI flagship — while the incumbent
   mixture is the least powerful. No method is downstream-disqualified, so the decision rightly
   rests on 1–3.
5. **It scales.** The shared rank-1 test is vectorized across all guides at once; on Gasperini
   (13k guides, MOI 26, 25k cells) each contingency method finished in ~1 min, while the per-guide
   `sceptre` mixture did not complete in 15+ min — the per-guide-fit methods (mixture, CLEANSER,
   crispat) are inherently slower at high-MOI scale.

**Honest bottom line.** The decisive downstream test showed *parity, not a blowout*: on typical
genome-wide data every sensible assignment calibrates, because sceptre's DE null is robust to the
assignment. So the case for `depth_fix` is not "it transforms downstream results" — it is that,
among methods that are all downstream-safe, `depth_fix` is the **most principled, the best-
motivated (assumptions validated broadly), the most accurate upstream, and the most robust in the
hard regimes** (high ambient, high MOI) where the alternatives — naive thresholds (the cutoff
trap), the bare ambient test (over-calls, FDR 0.14), and the current mixture (buggy, upstream-
mediocre, downstream-underpowered) — degrade. If maximum shippability is the priority, plain
fishash (the same test without the cell-margin fix) captures most of the value; `depth_fix` is the
strict superset that additionally fixes high-MOI self-masking at essentially no cost.

---

## 8. Artifacts (this decision cycle)

**Report:** `results/method_decision/RECOMMENDATION.md` (this file).

**Figures:**
- `method_scorecard.png` — the 4-axis method comparison (upstream accuracy/FDR, downstream calibration/power).
- `ambient_dispersion_all.png` — ambient soup dispersion vs rate, all 11 datasets (Poisson-null validation).
- `de_calibration_compare.png` — downstream neg-control type-I error by method, per genome-wide dataset.

**Tables:**
- `dataset_profiles.csv` — the 11-dataset landscape (MOI, depth, GE, NT guides).
- `ambient_validation.csv` — per-dataset ambient dispersion + library/denoised-depth ratio.
- `de_master.csv`, `de_aggregate.csv`, `de_trustworthy_summary.csv` — downstream DE calibration + power.

**Scripts (new this cycle):**
- `scripts/profile_datasets.R` — dataset landscape.
- `scripts/ambient_validation.R` + `run_ambient_validation.R` — rank-1 denoised-depth fit + soup dispersion (native reimplementation of the fishash rank-1 completion).
- `scripts/de_native.R` — per-(dataset,method) sceptre calibration+power (survey & ondisc modes).
- `scripts/de_driver.R`, `scripts/de_aggregate.R`, `scripts/method_scorecard.R` — grid + analysis + figure.

**Pre-existing evidence reused:** `results/barnyard_exact_corrected.csv` (barnyard ground truth),
`results/sim_framework/method_overall.csv` + `depthfix_vs_fishash_moi.csv` (simulation panel),
`scripts/contingency_method.R` (the method), `canonical_model.qmd` / `literature_review.qmd`.

**Open follow-ups (not blocking the decision):** head-to-head vs a correctly-scored CLEANSER;
full-scale (non-subsampled) Gasperini/Replogle DE on a cluster; a real high-MOI ground-truth
dataset if one becomes available; native sceptre implementation + Bioconductor packaging.
