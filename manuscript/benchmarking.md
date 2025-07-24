# Data import from Cell Ranger

## Competitors
- **CLEANSER**: Direct Cell Ranger import via `cr2cleanser` utility
- **crispat**: Direct import via `create_anndata_from_cellranger()`
- **pySpade**: Direct import from Cell Ranger `.h5` files

## Datasets
- Gasperini (~200K cells); concern is that this dataset is getting old. Don't we have any more up-to-date examples?
- Replogle RPE1 essential-wide (~600K cells)?
- Replogle K562 essential-wide (~600K cells)?
- Replogle K562 genome-wide (~1.2M cells)?

## Evaluation criteria
- Not clear to me

# Guide assignment

## Competitors
- **pertpy**: Poisson-Gaussian mixture model with GPU acceleration
- **CLEANSER**: Mixture model with multi-CPU processing via `concurrent.futures`
- **crispat**: A representative mixture model method among the several considered (with multi-CPU processing via Dask)

## Datasets
- Simulated (scDesign3?)
- Gasperini (~200K cells); concern is that this dataset is getting old. Don't we have any more up-to-date examples?
- Replogle RPE1 essential-wide (~600K cells)?
- Replogle K562 essential-wide (~600K cells)?
- Replogle K562 genome-wide (~1.2M cells)?

## Evaluation criteria
- Accuracy (on simulated data)
- crispat's metrics (on real data)
- RAM and runtime (how exactly to measure across parallelization frameworks?)

Comment: The more we can borrow of crispat's evaluation framework, the better. Or, do we even have to benchmark guide assignment accuracy, given crispat's extensive efforts to this end? Could we get away with just benchmarking computational performance?

# DE testing

## Competitors
- **pertpy**: PyDESeq2 (multi-CPU processing with `joblib`) or EdgeR (no parallelization)
- **FR-Perturb**: elastic net regression plus permutation p-values (randomized partial SVD for large datasets and manual parallelization via batch processing functionality)
- **Mixscale**: Weighted negative binomial GLM, perturbation scoring (no parallelization)
- **pySpade**: Hypergeometric test, t-test, gamma distribution significance (multi-CPU processing with `multiprocessing.Pool`)

## Datasets
- Simulated?
- Gasperini (~200K cells); concern is that this dataset is getting old. Don't we have any more up-to-date examples?
- Replogle RPE1 essential-wide (~600K cells)?
- Replogle K562 essential-wide (~600K cells)?
- Replogle K562 genome-wide (~1.2M cells)?

## Evaluation criteria
- Type-I error control and power
- RAM and runtime (how exactly to measure across parallelization frameworks?)
