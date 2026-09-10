# gRNA simulator

Generates synthetic gRNA UMI count matrices with **known ground-truth
perturbation labels**, used to benchmark guide-assignment methods (crispat,
pertpy, cleanser, sceptre) against truth. Each simulated dataset is written to
`guide_assignment/input_data/{dataset_name}/` in every method's input format,
alongside the true perturbation matrix that the methods are scored against
(see `../sceptre-nb/bugfix-sims-metrics.R`).

This is a standalone R data-generation stage — it *produces* the input data
that the Nextflow guide-assignment pipeline later consumes. It is not part of
the pipeline itself.

## The engine

### `sims-sum-of-three.R`
The simulator. Defines two functions plus a validation plot; sourced by every
driver script.

**Generative model** (`grna_simulator_iid_sum_process`) — each guide's UMI
count in each cell is a **sum of three independent pieces**:

1. **Exogenous background** — `Pois(lib_size · mu_exo)`. A low ambient floor
   present in every cell.
2. **Endogenous background** — in a random `prob_endo` fraction of cells, add
   `NB(size = theta_endo, mu = lib_size · mu_endo)`. Models sporadic high
   background (the contaminant that makes assignment hard).
3. **True perturbation** — in a random `prob_pert` fraction of cells, add
   `NB(size = theta_pert, mu = lib_size · mu_pert)`. These cells are the
   ground-truth positives, recorded in `true_perts`.

Everything scales by `lib_sizes_scaled` (real per-cell library sizes,
normalized to mean 1), so simulated depth tracks real cell-to-cell variation.
Counts are assembled as sparse-matrix triplets for efficiency.

This model is inspired by CellBender and fishash. We are moving to directly using the fishash simulation code, so this will likely be deprecated soon.

**Orchestrator** (`make_grna_simulation_sum_process`) —
- Samples `num_cells` real cells (with a guide detected) from a real sceptre
  object, reusing their library sizes and covariates.
- Loops over a **named `params_list`** (names become guide-name prefixes),
  generating `num_guides_per_param` guides per parameter set.
- Forces guide 1 to `non-targeting`; assignment formula optionally includes
  response covariates (`use_response_covariates`).
- Writes the matrix in four method formats — `.h5ad` (crispat, pertpy),
  `.mtx` (cleanser), and a full sceptre object — plus `true_pert_matrix.rds`,
  `cell_scaling_and_locs.rds`, and `sim_params.rds`.
- Prints an `scp` line to ship the dataset to HPC3.

**Validation** (`plot_umi_histogram_real_vs_sim`) — puts a real guide's UMI
histogram side-by-side with a simulated one (log1p x-axis, log10 y-axis,
simulated bars stacked perturbed/non-perturbed) as a visual check that the
simulated marginal count distribution resembles reality.

## The drivers

Each driver sources the engine, loads the real Replogle sceptre object
(`replogle-rd7/sceptre-pipeline`) for cell covariates/library sizes, defines a
`params_list`, calls `make_grna_simulation_sum_process`, and ends with a
real-vs-sim histogram sanity check. All three use 50,000 cells and 100 guides
per parameter set. These are exactly the three datasets scored by
`bugfix-sims-metrics.R`.

### `grna-sims-sum_2np_3p.R` → `sims_sum_2np_3p`
**2 non-perturbed regimes × 3 perturbation strengths** (6 guide groups).
- NP: low-variance (`theta_endo = 10`) and high-variance (`theta_endo = 1`).
- P: `mu_pert` = 500 / 750 / 1000 with `theta_pert` = 1 / 2.5 / 5.
- Sweeps how perturbation strength interacts with background variance.

### `grna-sims_sum_1np_3p_disp.R` → `sims_sum_1np_3p_disp`
**1 non-perturbed regime × 3 perturbation dispersions** (3 guide groups).
- NP: high-variance only.
- P: fixed `mu_pert = 1000`, varying `theta_pert` = 10 / 1 / 0.1 (low / med /
  high dispersion).
- Isolates the effect of perturbation-signal dispersion at fixed mean.

### `grna-sims_sum_repeat_old.R` → `sims_sum_repeat_old`
**Reproduces an older simulation setup** (3 guide groups), seed `123`.
- NP: Poisson-like (`prob_endo = 0`, so no endogenous NB spikes; only exogenous
  background with `mu_exo = 0.07`).
- P: fixed `theta_pert = 0.126` and `prob_pert = 0.004`, with `mu_pert` = 970 /
  242.5 / 121.25 (large / med / small).
- These replicate the only simulations presented in my defense. 

## Other files

- `grna-sims-old-regr-ver-writeup.Rmd` — writeup for an older regression-based
  version of the simulator (superseded by the sum-of-three model here).

## Usage

Run a driver directly in R/RStudio (they `rm(list=ls())` at the top and source
`~/.Rprofile`); each writes to
`guide_assignment/input_data/{dataset_name}/` and prints an `scp` command to
copy the dataset to HPC3 for the benchmarking pipeline.
