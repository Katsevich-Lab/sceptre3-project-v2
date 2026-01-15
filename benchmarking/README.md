# CRISPR Benchmarking Pipelines

This repository contains Nextflow pipelines for benchmarking CRISPR perturbation screening methods. The pipelines are designed to run on both local machines and HPC clusters (SGE) with support for conda and Singularity containers.

## Overview

Two complementary benchmarking pipelines:

1. **Guide Assignment Benchmarking** (`guide-assignment/`) - Tests methods that assign gRNAs to cells
2. **Association Testing Benchmarking** (`association/`) - Tests methods that detect gene-perturbation associations GIVEN fixed assignments

This separation allows independent evaluation of:
- How well methods assign gRNAs to cells (guide assignment)
- How well methods detect real associations (association testing)

---

# Guide Assignment Benchmarking

## Overview

Nextflow pipeline to benchmark guide assignment methods on CRISPR screening datasets. The pipeline provides a standardized framework for comparing different guide assignment algorithms with consistent input/output formats and resource management.

**Currently supported methods**: crispat
**Currently supported datasets**: example

**Framework ready for**: Adding new guide assignment methods and datasets

## Directory Structure

```
benchmarking/guide-assignment/
├── create-example-gRNA-data.py    # Script to create example data
├── example-data-requirements.txt  # Python requirements for data creation script
├── run-benchmarking.sh           # Launch script for pipeline
└── guide-assignment-pipeline/    # Benchmarking pipeline
    ├── main.nf                   # Main pipeline workflow
    ├── nextflow.config          # Base configuration parameters
    ├── bin/
    │   └── run_crispat.py       # Method-specific wrapper scripts
    ├── modules/
    │   └── crispat/
    │       ├── main.nf          # Nextflow process definition
    │       └── environment.yml  # Conda environment specification
    ├── configs/
    │   └── {run_id}_config.csv  # Run-specific method/dataset combinations
    └── conda-cache/             # Cached conda environments (gitignored)
```

## Pipeline Architecture

### 1. Configuration System
- **Run-based configs**: Each run uses `configs/{run_id}_config.csv`
- **Method-dataset matrix**: CSV defines which methods run on which datasets
- **Resource specification**: CPU/memory requirements per combination

### 2. Standardized Data Flow
1. **Input**: Method-specific input files, could be h5ad format or something else (`grna_counts_{method}.h5ad`)
2. **Processing**: Each method runs in isolated conda environment
3. **Output**: Standardized CSV format (`assignments_{method}_{dataset}.csv`)
4. **Publishing**: All results in single output directory for easy comparison

### 3. Method Integration Pattern
Each method follows a consistent pattern:
- **Wrapper script** in `bin/run_{method}.py`: Handles method-specific logic
- **Nextflow module** in `modules/{method}/main.nf`: Defines process and resources
- **Conda environment** in `modules/{method}/environment.yml`: Dependencies

## Usage

### Running the Pipeline

```bash
# Basic usage with run identifier
nextflow run main.nf --run_id "example"

# This uses configs/example_config.csv and outputs to:
# ~/data/projects/sceptre3/benchmarking/guide_assignment/outputs/example/
```

### Configuration File Format

Create `configs/{run_id}_config.csv`:
```csv
method,dataset,cpus,memory
crispat,example,1,4.GB
crispat,dataset2,2,8.GB
```

### Input Data Requirements

For each dataset and method combination, provide:
```
input_data/
└── {dataset_id}/
    └── grna_counts_{method}.h5ad  # Method-specific input format
                                   # (could be h5ad or something else)
```

### Output Structure

All results are published to a single directory:
```
outputs/{run_id}/
├── assignments_crispat_example.csv
├── assignments_crispat_dataset2.csv
└── assignments_newmethod_example.csv
```

**Standardized output format**:
- Filename: `assignments_{method}_{dataset}.csv`
- Columns: `cell_id`, `grna_id`

## Adding New Methods

### 1. Create Conda or Singularity Environment

Create `modules/{method}/environment.yml` with all dependencies your method needs. See `modules/crispat/environment.yml` as an example - it specifies Python version, pandas for data handling, and method-specific dependencies.

### 2. Create Method Wrapper Script

Create `bin/run_{method}.py` that runs your method and transforms output to the standardized format. See `bin/run_crispat.py` as an example - it loads input data, runs the assignment algorithm, and outputs a CSV with `cell_id` and `grna_id` columns. Develop and debug method wrapper script directly in Python on example data, within a conda environment created based on step 1.

### 3. Create Nextflow Module

Create `modules/{method}/main.nf` that defines the Nextflow process. See `modules/crispat/main.nf` as an example - it specifies resource requirements, conda environment, input/output structure, and calls your wrapper script.

### 4. Update Main Pipeline

In `main.nf`, add your method to the import statements and branching logic. Follow the pattern used for crispat to route your method through the appropriate processing channel.

## Adding New Datasets

1. **Prepare input data**: Create `input_data/{dataset_id}/` directory
2. **Method-specific formats**: Add input file in correct format for each method
3. **Update config**: Add dataset to your run's config CSV

---

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
  --out_base_dir "$HOME/data/projects/sceptre3/benchmarking/association/outputs"
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

---

## Pipeline Relationship

The two pipelines work together in the benchmarking workflow:

1. **Guide Assignment** produces: `assignments_{method}_{dataset}.csv` with cell-to-gRNA assignments
2. **Association Testing** uses fixed assignments to evaluate: How well methods detect real gene-perturbation associations

This separation enables:
- Testing assignment methods independently
- Testing association methods with controlled assignment quality
- Fair comparison across different methodological approaches
