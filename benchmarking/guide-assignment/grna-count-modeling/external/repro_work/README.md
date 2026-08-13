# Fishash Table 2 reproduction

This directory reproduces the four-sample Liu et al. barnyard benchmark used in
Table 2 of the Fishash preprint. The original reproduction covered Fishash and
Geomux. `run_extended_table2.sh` adds SCEPTRE mixture and both CLEANSER modes.

## Authorities and frozen versions

- Fishash analysis commit `ea8af91587d07f1b895e515fe24d5db10a0d539a`:
  `21_process_barnyard_data.R`, `bin/run_sceptre_mixture.R`, the CLEANSER
  process in `main.nf`, and `62_plot_barnyard_results.R`.
- Pinello CRISPR Pipeline commit `43c89c0806208a08a47c3e5bea5071443f9daf74`:
  `bin/assign_grnas_sceptre.R` and the native SCEPTRE/CLEANSER Nextflow modules.
- SCEPTRE 0.10.3, the version reported by the Fishash preprint.
- CLEANSER 1.2.1, the patch release used by the current Pinello container and
  locally available environment (the preprint reports 1.2).

The SCEPTRE runner uses the matched full gene-expression count matrix, `moi =
"high"`, default analysis parameters, and mixture assignment, matching the
Fishash authors' runner. The CLEANSER runner uses the authors' native CLI call:
both `--cs` and `--dc`, package-default sampling settings, and eight parallel
runs. The authors' command omitted `--seed`; CLEANSER 1.2.1 chooses a random
default seed at process import, so their exact MCMC draw is not recoverable from
the repository. This reproduction fixes seed `20260810` (override with
`CLEANSER_SEED`) and records any numerical discrepancy with the published row.

## CLEANSER threshold discrepancy

The preprint says that Table 2 uses posterior cutoffs 0.80 for CROP-seq
CLEANSER and 0.50 for direct-capture CLEANSER. The current public analysis runs
both modes on every sample but `subset_methods()` selects cutoff 0.95 for both
rows. A fresh first-sample run resolves the conflict empirically: cutoff 0.80
reproduces the published CROP-seq value closely (0.9412 versus 0.9414), whereas
0.95 does not (0.9332). The primary reproduction therefore uses the stated
0.80/0.50 Table 2 thresholds and reports 0.95 only as a public-code sensitivity
analysis.

## Run

First run `build_inputs.py`. Install SCEPTRE 0.10.3 into an isolated R library,
then execute:

```bash
RAW_DIR=/path/to/GSE272457 \
SCEPTRE_LIB=/path/to/r-library-with-sceptre-0.10.3 \
CPUS=8 \
./run_extended_table2.sh
```

The long-running posterior matrices and detailed logs are ignored by git. The
tracked `extended_table2_reproduction.csv` contains the exact-cohort accuracy,
published-value discrepancy, and counts of non-fatal, convergence, and
divergence warning blocks from each CLEANSER run.

## Reproduction result

Using the fixed inputs and settings above:

- SCEPTRE reproduces all four published accuracies at four decimals.
- CLEANSER reproduces six of eight published accuracies at four decimals; the
  two differences are 0.00018 (CROP-seq mode, 0-hour CROP-seq sample) and
  0.00049 (direct-capture mode, 0-hour Direct Capture sample).
- All 12 primary accuracies are within 0.001 of Table 2, and all 12 reported
  standard errors reproduce at four decimals.
- The public analysis helper's 0.95 CLEANSER cutoff does not regenerate the
  Table 2 rows, confirming that the table used the stated 0.80/0.50 cutoffs.

CLEANSER completed all guide fits, but CmdStan emitted non-fatal proposal-domain
warnings for 199 or 200 of 200 guides per run and convergence/divergence warning
blocks for a smaller subset. These counts are preserved in the summary and
should accompany any claim of numerical reproduction.
