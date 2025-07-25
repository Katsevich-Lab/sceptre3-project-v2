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

### 1. Create Conda Environment

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