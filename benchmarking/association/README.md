# Association Testing Benchmarking Pipeline

Nextflow pipeline for benchmarking CRISPR perturbation association testing methods with controlled gRNA assignments.

## Overview

This pipeline tests different methods for detecting gene-perturbation associations GIVEN fixed/known gRNA assignments. This allows fair comparison of association testing methods independently of assignment quality.

**Supported methods:**
- **pertpy** (Python, conda)
- **sceptre** (R, Singularity)
- **mixscale** (R, Singularity)

## Directory Structure

```
association/
├── association-pipeline/          # Main Nextflow pipeline
│   ├── main.nf                   # Workflow with method routing
│   ├── nextflow.config           # Executor and resource config
│   ├── bin/                      # Method wrapper scripts
│   │   ├── run_pertpy.py
│   │   ├── run_sceptre.R
│   │   └── run_mixscale.R
│   ├── modules/                  # Method-specific Nextflow processes
│   │   ├── pertpy/
│   │   ├── sceptre/
│   │   └── mixscale/
│   └── configs/                  # Run-specific method/dataset configs
│       └── example_config.csv
├── input_data/                   # Association testing inputs
│   └── {dataset}/
│       └── {method}/
│           ├── response_matrix.*  # Gene expression
│           ├── grna_assignments.* # Fixed assignments
│           └── cell_covariates.csv
├── outputs/                      # Results by run_id
│   └── {run_id}/
│       ├── results_{method}_{dataset}.*
│       └── report.html
├── run-association-benchmarking.sh          # HPC launcher
└── local-run-association-benchmarking.sh    # Local launcher
```

## Data Input Organization

**Pattern:** `input_data/{dataset_id}/{method}/`

Each method subdirectory contains:
- **Gene expression data:** response_matrix.h5ad (pertpy) or .rds (R methods)
- **gRNA assignments:** Fixed/ground-truth assignments to test
- **Cell covariates:** Metadata for batch effects, QC, etc.
- **Target mapping:** (sceptre only) gRNA to target gene mapping

## Usage

### 1. Prepare Input Data

Create dataset directories following the structure:
```bash
mkdir -p input_data/my_dataset/{pertpy,sceptre,mixscale}
```

Place method-specific input files in each subdirectory.

### 2. Create Config File

Create `configs/my_run_config.csv`:
```csv
method,dataset,cpus,memory
pertpy,my_dataset,4,16GB
sceptre,my_dataset,2,8GB
mixscale,my_dataset,4,16GB
```

### 3. Run Pipeline

**On HPC cluster:**
```bash
./run-association-benchmarking.sh
# Edit RUN_ID variable in script first
```

**Locally:**
```bash
./local-run-association-benchmarking.sh
# Edit RUN_ID variable in script first
```

**Or directly:**
```bash
nextflow run association-pipeline/main.nf \
  --run_id "my_run" \
  --out_base_dir "$HOME/data/projects/sceptre3/benchmarking/association/pos-control/outputs"
```

## Adding New Methods

1. **Create module directory:** `modules/new_method/`
2. **Define process:** `modules/new_method/main.nf`
3. **Setup environment:** `environment.yml` (conda) or `.def` (Singularity)
4. **Create wrapper script:** `bin/run_new_method.{py,R}`
5. **Update main.nf:** Add import and branching logic
6. **Update nextflow.config:** Add label for environment

## Configuration

### Nextflow Config

- **Executor:** SGE (default) or local (with `-profile local`)
- **Paths:** `dataset_base_dir`, `out_base_dir` in nextflow.config
- **Containers:** Singularity/Apptainer enabled
- **Resources:** CPU/memory specified per method-dataset in CSV

### Method-specific Settings

Each method can have custom:
- CPU/memory requirements (in config CSV)
- Container/conda environment
- Input data format
- Output format

## Output

Results are published to `outputs/{run_id}/`:
- `results_{method}_{dataset}.{csv,rds}` - Association test results
- `report.html` - Nextflow execution report
- `trace.tsv` - Resource usage trace
- `timeline.html` - Execution timeline

## Development Status

**Current state:** Framework complete with placeholder implementations

**TODO:**
- Implement actual pertpy association testing in `bin/run_pertpy.py`
- Implement actual sceptre association testing in `bin/run_sceptre.R`
- Implement actual mixscale association testing in `bin/run_mixscale.R`
- Build mixscale Singularity container from `modules/mixscale/mixscale.def`
- Create or copy sceptre.sif container to `modules/sceptre/`
- Prepare test datasets in `input_data/`

## Related Pipelines

- **Guide assignment:** `benchmarking/guide-assignment/` - Tests gRNA assignment methods
- **Association testing:** `benchmarking/association/` (this pipeline) - Tests association methods with fixed assignments
