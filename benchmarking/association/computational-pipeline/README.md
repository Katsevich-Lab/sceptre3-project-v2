# Computational Benchmarking Pipeline

This pipeline benchmarks the **computational performance** (runtime and memory usage) of three CRISPR perturbation analysis methods:
- SCEPTRE
- Mixscale
- FR-Perturb

## Key Differences from neg-control-pipeline

Unlike the `neg-control-pipeline` which tests **statistical power** (false discovery rates), this pipeline focuses on **computational scalability**:

- **Fair comparison**: All methods test identical discovery pairs with equal CPU/memory allocation
- **Controlled parallelism**: Thread-level parallelism disabled for SCEPTRE/Mixscale (process-level only)
- **Multiple dataset sizes**: Tests scalability at 100, 500, 1000, 2000, 4000 genes
- **Performance metrics**: Nextflow trace reports capture runtime, memory, CPU utilization

## Fair Benchmarking Configuration

**Parallelism control:**
- **SCEPTRE**: `OMP_NUM_THREADS=1` (process-level parallelism via `NCPUS`)
- **Mixscale**: `OMP_NUM_THREADS=1` (process-level parallelism via `NCPUS`)
- **FR-Perturb**: Thread-level parallelism set to allocated CPUs (`OMP_NUM_THREADS=${task.cpus}`)

**Resource allocation:**
- All methods receive identical CPU and memory allocations per dataset size
- See `configs/*_config.csv` for specific allocations

## Directory Structure

```
computational-pipeline/
├── main.nf                    # Main Nextflow workflow
├── nextflow.config            # Configuration (paths, executor, containers)
├── modules/
│   ├── sceptre/main.nf       # SCEPTRE process definition
│   ├── mixscale/main.nf      # Mixscale process definition
│   └── frperturb/main.nf     # FR-Perturb process definition
├── bin/
│   ├── run_sceptre.R         # SCEPTRE wrapper script
│   ├── run_mixscale.R        # Mixscale wrapper script
│   └── run_frperturb.py      # FR-Perturb wrapper script
└── configs/
    ├── test_100genes_config.csv       # 100 genes, 10 CPUs, 8GB
    ├── small_500genes_config.csv      # 500 genes, 20 CPUs, 16GB
    ├── medium_1000genes_config.csv    # 1000 genes, 20 CPUs, 32GB
    ├── large_2000genes_config.csv     # 2000 genes, 40 CPUs, 64GB
    └── xlarge_4000genes_config.csv    # 4000 genes, 40 CPUs, 128GB
```

## Setup

### 1. Create Computational Datasets

The data preparation script creates datasets at multiple scales:

```bash
cd association/data-preparation
Rscript make_computational_replogle-rd7.R
```

This creates 5 datasets:
- `replogle-rd7_computational_100genes`
- `replogle-rd7_computational_500genes`
- `replogle-rd7_computational_1000genes`
- `replogle-rd7_computational_2000genes`
- `replogle-rd7_computational_4000genes`

### 2. Move Datasets to Computational Directory

The script outputs to `neg-control/input_data/` by default. Move them:

```bash
mkdir -p ~/data/projects/sceptre3/benchmarking/association/computational/input_data/
mv ~/data/projects/sceptre3/benchmarking/association/neg-control/input_data/replogle-rd7_computational_* \
   ~/data/projects/sceptre3/benchmarking/association/computational/input_data/
```

## Running the Pipeline

### Test Run (100 genes)

```bash
cd association/computational-pipeline
nextflow run main.nf --run_id test_100genes
```

### Full Runs

```bash
# Small scale (500 genes)
nextflow run main.nf --run_id small_500genes

# Medium scale (1000 genes)
nextflow run main.nf --run_id medium_1000genes

# Large scale (2000 genes)
nextflow run main.nf --run_id large_2000genes

# Extra large scale (4000 genes)
nextflow run main.nf --run_id xlarge_4000genes
```

### Local Testing (Laptop)

For testing on a local machine:

```bash
nextflow run main.nf --run_id test_100genes -profile local
```

## Output Files

Results are written to:
```
~/data/projects/sceptre3/benchmarking/association/computational/outputs/<run_id>/
```

**Method outputs:**
- `association_computational_sceptre_<dataset_id>.csv` - SCEPTRE results
- `association_computational_mixscale_<dataset_id>.csv` - Mixscale results
- `frperturb/frperturb_results_<dataset_id>.*` - FR-Perturb results (4 files)

**Performance metrics:**
- `trace.txt` - Detailed per-process runtime, memory, CPU metrics
- `timeline.html` - Visual timeline of task execution
- `report.html` - Nextflow execution report

## Analyzing Performance Metrics

The `trace.txt` file contains the computational benchmarking data:

```r
library(tidyverse)

# Read trace file
trace <- read.table("computational/outputs/test_100genes/trace.txt",
                    header=TRUE, sep="\t")

# Compare methods
trace %>%
  mutate(method = gsub("_COMPUTATIONAL.*", "", process)) %>%
  group_by(method) %>%
  summarise(
    median_runtime_sec = median(realtime) / 1000,
    max_memory_GB = max(peak_rss) / (1024^3),
    mean_cpu_util = mean(`%cpu`)
  )
```

**Key trace fields:**
- `realtime` - Wall-clock time (milliseconds)
- `peak_rss` - Peak memory usage (bytes)
- `%cpu` - CPU utilization percentage
- `cpus` - Number of CPUs allocated

## Fair Comparison Checklist

Before running benchmarks, verify:

1. ✅ All three methods test identical discovery pairs (full Cartesian product)
2. ✅ Equal CPU allocation across methods (check config CSV)
3. ✅ Equal memory allocation across methods (check config CSV)
4. ✅ Thread-level parallelism controlled (SCEPTRE/Mixscale: `OMP_NUM_THREADS=1`, FR-Perturb: set to `task.cpus`)
5. ✅ Same input data (assignments, expression matrices, covariates)

## Configuration Files

Each config CSV has format:
```csv
method,dataset,cpus,memory
sceptre,replogle-rd7_computational_100genes,10,8.GB
mixscale,replogle-rd7_computational_100genes,10,8.GB
frperturb,replogle-rd7_computational_100genes,10,8.GB
```

**Resource allocation strategy:**
- Test runs (100 genes): 10 CPUs, 8GB RAM
- Small (500 genes): 20 CPUs, 16GB RAM
- Medium (1000 genes): 20 CPUs, 32GB RAM
- Large (2000 genes): 40 CPUs, 64GB RAM
- Extra large (4000 genes): 40 CPUs, 128GB RAM

Adjust based on available cluster resources and observed memory usage.

## Troubleshooting

**Memory errors:**
- Increase memory allocation in config CSV
- Check `trace.txt` for `peak_rss` to see actual memory usage

**Slow execution:**
- Verify parallelism settings are correct
- Check CPU utilization in `trace.txt` (`%cpu` field)
- Ensure SGE job has access to requested CPUs

**Method failures:**
- Check Nextflow logs in `work/` directory
- Review SGE logs in `nf-logs/`
- Verify input data completeness

## Notes

- The pipeline reuses `make_neg_control_replogle()` function since it already creates full Cartesian product discovery pairs
- FR-Perturb outputs 4 files (results, log, model, scores) - all are published to `frperturb/` subdirectory
- Nextflow trace metrics are more accurate than self-reported method timings
- Use `work/` directory contents for detailed debugging only (automatically cleaned between runs)
