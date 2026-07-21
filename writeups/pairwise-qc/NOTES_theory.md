# Pairwise QC for sceptre — theoretical core

## Current method (from source)
For each (grna_target g, response gene r) pair, sceptre builds the 2×2 table
of {treatment vs control} × {gene nonzero vs zero} and keeps the pair iff

    n_nonzero_trt   >= 7   AND   n_nonzero_cntrl >= 7          (defaults, run_qc)

- `n_nonzero_trt`   = # treatment cells (carry g) with nonzero expression of r   → INTERIOR cell (a)
- `n_nonzero_cntrl` = # control cells (NT or complement) with nonzero expr of r  → ≈ column margin

Test statistic (src/low_level_full_test.cpp, precomputation_functions.R):
    z = ( Σ_{i∈trt} a_i ) / sqrt( Σ_{i∈trt} w_i − Σ_k (Σ_{i∈trt} D_{k,i})^2 )
    a_i = (y_i − μ_i)/(1 + μ_i/θ)   [NB working residual],  μ_i = fitted NB mean
Null distribution: resample the treatment index set (permutation low-MOI / CRT high-MOI);
gene expression y, μ, a, w, D are HELD FIXED. Only "which cells are treatment" varies.

## The flaw (why it's "circular")
n_nonzero_trt = Σ_{i∈trt} 1{y_i>0}.  This is the SAME functional form as the test numerator
Σ_{i∈trt} a_i — a sum over the treatment cells of a per-cell expression signal. Since 1{y_i>0}
is monotonically related to a_i (=residual), **n_nonzero_trt is a coarse proxy for the test
statistic itself.**
- A knockdown drives y_i→0 in treatment cells ⇒ n_nonzero_trt ↓ AND z ↓ (more negative) together.
  So the filter removes pairs *because* they show the effect. Strong knockdown ⇒ likely dropped.
- Directional: repression lowers n_nonzero_trt (drop risk); activation raises it (kept). The filter
  is biased AGAINST detecting repression — the dominant mode in CRISPRi.

## Independent-filtering criterion (Bourgon, Gentleman, Huber 2010 PNAS; DESeq2)
Two-stage testing preserves type-I/FDR control **iff the filter statistic is marginally
independent of the test statistic UNDER THE NULL** (then conditional null = unconditional null).
Under the alternative they SHOULD be correlated — that is where the power gain comes from.
DESeq2 filters on `baseMean` (mean normalized count over ALL samples) — a function of the data
that is independent of the per-gene test statistic under the null — and optimizes the threshold
to maximize rejections at target FDR.

Mapping to sceptre's resampling test: a statistic is ancillary (⇒ valid filter) iff it is a
function only of quantities held fixed under the resampling null = (gene expression vector y,
covariates Z, treatment-group size n_trt). Margins qualify; interior cells do not.
- `n_nonzero_trt` (interior) VARIES across the null ⇒ NOT independent ⇒ invalid filter. ✗
- gene expression margin (n_nonzero over control/all, baseMean) — function of y ⇒ ancillary. ✓
- treatment size n_trt — fixed ⇒ ancillary. ✓

## Proposed method — margin-based independent filter
Replace the OBSERVED treatment count by its NULL EXPECTATION (the DESeq2-baseMean analog for
this design). Let p̂ = baseline expression rate estimated from control cells (unaffected by g):
        p̂ = n_nonzero_cntrl / n_cntrl
Define the expected number of informative treatment cells
        m = n_trt · p̂
Keep the pair iff
        m >= 7            AND     n_nonzero_cntrl >= 7
Only the TREATMENT side changes (observed → expected). n_nonzero_cntrl is already a margin,
so the control side is unchanged. Threshold 7 keeps identical semantics: m is an unbiased
estimate of E[n_nonzero_trt | null], so "expected ≥ 7" is the independent-filtering version of
the current "observed ≥ 7".

Why m is the right statistic:
1. Ancillary (function of y over control cells × fixed n_trt) ⇒ exactly independent of z under
   the null ⇒ type-I/FDR preserved (Bourgon). Verified: cor(m, z_null)≈0 vs cor(n_nonzero_trt,z)>0.
2. Predicts POWER (few expected informative cells ⇒ genuinely underpowered) — so it still filters
   junk (e.g. rare gene × few cells) — but NOT effect (a strong knockdown no longer self-censors).
3. Governs the resampling null's discreteness (the null draws have ≈ m nonzero trt cells on
   average), so it is the correct CALIBRATION guard too — better than the single observed draw.

Caveat (transparent): a gene silent in control (p̂≈0 ⇒ m<7) that is strongly ACTIVATED by the
perturbation would be filtered. Rare in CRISPRi (repression-dominated); such pairs are near-
degenerate for the resampling null anyway. If activation of silent genes is of interest, base p̂
on the union or add a small observed-count OR-guard.

Optional (DESeq2-style): treat 7 as a default but allow choosing the m-threshold that maximizes
BH rejections at target FDR (independent_filtering=TRUE), since m ⟂ p under null makes this valid.

## Evidence
- SIM (real sceptre, ground truth): current-filter retention of true effects collapses with
  knockdown strength (→0 by κ=0.95); margin filter flat at 100% except genuinely underpowered
  cell (m=6). 70% of detectable effects dropped by current filter; margin rescues ~80%.
  BH<0.1 true discoveries 159→502 (two-sided) / 165→524 (left); realized FDR stays <10%.
  Null calibration identical (0.045/0.051). cor(n_nonzero_trt, z_null)=+0.26 vs cor(m,z)=−0.01.
- REAL (a549 CRISPRi): null pairs dropped by current filter (n_nonzero_trt<7) have mean z=−0.305
  (repression-leaning; −1.08 in the 0–2 bin) vs +0.079 for kept; margin filter shows no gradient
  (cor 0.106 vs 0.011). Calibration preserved under margin filter.
