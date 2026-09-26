> ## ⚠️ UNREVIEWED AI-GENERATED DRAFT — DO NOT CITE OR RELY ON
>
> This file was written by Claude as a summary of a working conversation. **Louis
> has not reviewed it.** Nothing here should be taken as fact, as a decision that
> has been made, or as an accurate description of what the code does. It is a work
> in progress and may contain errors, misattributed reasoning, and claims that were
> discussed but never agreed.
>
> **The authoritative account of the design is Louis's own:**
> [`manuscript/reports/meeting-2026-09-17/meeting-2026-09-17.Rmd`](../../../manuscript/reports/meeting-2026-09-17/meeting-2026-09-17.Rmd).
> Where this file and that one disagree, that one is right. This file exists to
> record what was built and measured *after* that report, and to note where the
> implementation has run into trouble.

# Computational-scaling datasets

Simulated gRNA count matrices for measuring how guide-assignment methods' runtime
and memory grow with dataset size. **Cost only** — see "What these are not".

Generator: `fishash::simulate_guidebender2` (GuideBender), used as-is. The code is
in the "GuideBender scaling datasets" section of `grna-sim-utils.R`;
`simulate-scaling-datasets.R` is a thin run script.

```
Rscript simulate-scaling-datasets.R gasperini                   # every rung
Rscript simulate-scaling-datasets.R replogle 1 2                # rungs 1 and 2
Rscript simulate-scaling-datasets.R gasperini --calibrate-only  # no datasets written
```

Generation is deterministic given its seeds, so running this on another machine
with the same real matrices reproduces the same datasets rather than copying them.

## The design, in brief

Argued in full in the meeting report; summarized here only so the code is legible.

Notation follows the report: **NNZ** is the total number of nonzero entries and
per-guide NNZ the number for a given guide; **NTP** ("number of true
perturbations") is the total number of perturbed (guide, cell) pairs, and per-guide
NTP the number of cells actually expressing a given guide.

Four realism requirements:

1. MOI stays fixed.
2. Per-guide NTP stays fixed.
3. The per-cell **signal** UMI total stays fixed.
4. The per-cell **error** UMI total stays fixed (endogenous ambient plus exogenous
   sources such as chimeras and barcode swaps).

Requirements 1 and 2 force **N/G fixed**, since counting perturbations two ways
gives `G x per-guide NTP = N x MOI`. No subset of a real matrix can hold both —
downsampling guides scales MOI by the fraction kept — which is why these are
simulated rather than downsampled.

Requirements 1 and 4 pin per-cell NNZ: MOI-many entries are nonzero from signal,
and the rest from error counts.

Everything being matched is a **first moment**.

## How the invariants map to GuideBender

Nothing is hardcoded. A regime is four inputs — the real dataset it mimics, the
count threshold at which an entry is called perturbed, the rung sizes, and the
number of cells used in calibration runs — and every target is measured from the
real matrix at run time, so changing the threshold changes all of them together.

| parameter | source |
|---|---|
| `n_cells` | the ladder (below) |
| `n_guides` | `round(n_cells / (N/G))`, N/G measured |
| `count_per_cell` | measured: median UMIs per cell |
| `hurdle_prob` | measured: fraction of cells with no perturbed guide |
| `guide_infection_alpha` | measured: from the spread of per-guide perturbed-cell counts |
| `moi` (lambda) | **calibrated per rung**, so realized MOI matches the real one |
| `snr` | **calibrated once**, so realized per-cell NNZ matches the real one |
| the remaining 8 | left at the values used in the fishash paper |

Two notes on the implementation:

`moi` is not MOI. It is the rate of the Poisson draw of *infections*, before zero
truncation, the hurdle, and repeated draws of the same guide, so it is solved for
numerically at every rung. `snr` is likewise solved for rather than derived,
because the closed forms assume a Dirichlet(1) guide abundance and single-UMI
ambient entries, and neither appears to hold: Gasperini's guide library measures a
Gini of 0.21 against Dirichlet(1)'s 0.50, and 17% of its sub-threshold nonzeros
carry more than one UMI.

`chunk_cells` is not only a memory setting. The exogenous-noise profile and the
depth rescaling are computed per chunk, so changing it changes the data. It is
pinned at 1000.

## The ladders

Sized by cleanser, the most expensive method. Its cost appears to be about
**11.4 ms per nonzero entry** — extrapolated from warm-up runs on small Replogle
guides and **not yet verified at these sizes**; cleanser's per-guide cost is
measurably superlinear, so the true rate may be higher. Under that assumption each
ladder spans **2 h to 8 h** of cleanser, i.e. 0.64M to 2.5M nonzeros.

| gasperini N | G | | replogle N | G |
|---|---|---|---|---|
| 11,000 | 694 | | 36,000 | 156 |
| 16,000 | 1,009 | | 49,000 | 212 |
| 22,000 | 1,388 | | 66,000 | 286 |
| 31,000 | 1,955 | | 89,000 | 385 |
| 43,000 | 2,712 | | 119,000 | 515 |

The 2-hour floor matters for Replogle. Error molecules collide onto the same guide
when guides are scarce, costing nonzeros; at a 30-minute floor its smallest rung
would be G = 31, where roughly 21 error molecules land among 31 guides and nonzeros
per cell fall by nearly half. At G = 156 the effect is mild. Gasperini is
insensitive either way — its error contribution is about 8 molecules per cell and
most of its nonzeros are signal.

## What has been measured

Everything fitted matches. The unfitted statistics do not. From
`sim_gasperini_g694_c11000`:

| | real | simulated |
|---|---|---|
| nonzeros per cell | 58.46 | 57.81 |
| MOI | 31.64 | 31.60 |
| perturbed cells per guide | 501.6 | 500.9 |
| median UMIs per cell | 564 | 563 |
| fraction of sub-threshold entries with 1 UMI | 0.828 | **0.518** |
| CV of perturbed cells per guide | 0.375 | **0.450** |
| true MOI (ground truth) | — | **49.7** |

Both mismatches appear to share a cause: targets are measured by thresholding real
counts, but parameters are derived as if the threshold were perfect.

- `Phi_cell = 1` makes the signal geometric, so many true perturbations fall below
  the threshold. Calibration then explains the real matrix's nonzeros as low-count
  *signal* rather than error, and the ground truth comes out denser than the
  thresholded MOI implies.
- `guide_infection_alpha` is derived from the spread of *thresholded* per-guide
  counts, which is inflated because lower-expression guides are detected less often.

Both could be calibrated numerically against the statistics they miss. Neither has
been, on the reasoning that neither moves the cost drivers and these datasets are
for cost. That reasoning is unreviewed.

Replogle additionally drifts about 19% in nonzeros per cell across its ladder, from
the collision effect above. Treat nonzeros per cell as a measured covariate rather
than a constant when fitting.

## What these are not

**Not valid for accuracy comparisons.** Density cannot match at small G — signal
alone puts MOI/G of entries nonzero — the per-entry error rate varies substantially
along each ladder, and the ground-truth matrix is denser than the real data's
thresholded MOI implies.

**Not needed for the fast methods.** fishash and fishash+ run on the full real
matrices in minutes, so they were timed there directly rather than on a ladder
(single runs on one laptop; about 9% run-to-run variation observed):

| | gasperini | replogle |
|---|---|---|
| fishash | 175 s, 5.6 GiB | 114 s, 6.6 GiB |
| fishash+ (Rcpp) | 29 s, 1.5 GiB | 16 s, 1.6 GiB |

These rungs are too small to time them usefully — fishash+ handles 2.5M nonzeros in
a few seconds, where container startup dominates.

## Outputs

Per dataset, under `<LOCAL_BENCHMARKING_DIR>/guide_assignment/input_data/`:

```
<dataset>/cleanser/grna_matrix.mtx    guides x cells; cleanser and fishash read this
<dataset>/crispat/grna_matrix.h5ad    cells x guides
<dataset>/pertpy/grna_matrix.h5ad     identical to crispat's
<dataset>/true_pert_matrix.rds        ground truth
<dataset>/sim_params.rds              every parameter used, plus the measured targets
```

plus `sim_scaling_calibration_<regime>.csv` (the solved lambda and snr) and
`sim_scaling_manifest_<regime>.csv` (the real matrix and every rung, measured with
the same function).

Methods and the input each reads are registered in `SCALING_METHODS`; adding one is
a single entry, and an input type with no writer raises rather than being skipped.

## Tests

```
Rscript tests/test-scaling-utils.R       # the scaling utilities
Rscript tests/test-simulate-scaling.R    # the generated data
```

52 checks, covering the format round-trip, the ray arithmetic, `simulate_regime`
being byte-identical to an explicit call with all 19 parameters named, calibration
and the real-data measurement defining MOI and nonzeros per cell identically, and
the method registry.

## Open questions, not decided

- Whether to calibrate `Phi_cell` and `guide_infection_alpha` numerically.
- How to handle Replogle's nonzeros-per-cell drift when fitting.
- Whether to add cell-subsets of the real matrices as a complementary,
  calibration-free way to vary N at fixed G.
- Verifying the 11.4 ms per nonzero cost constant before committing to the ladders.
- A test-driven rebuild was agreed but not started: six pure functions
  (`measure_counts`, `derive_params`, `build_rungs`, `solve_params`,
  `simulate_one`, `write_inputs`), with a round-trip test that recovers known
  parameters from simulated data.
