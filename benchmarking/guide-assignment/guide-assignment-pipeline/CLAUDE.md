# CRISPAT Benchmarking Pipeline

## Overview

Nextflow pipeline to benchmark guide assignment methods on CRISPR screening datasets. Currently supports crispat with framework for adding more methods.

## Directory Structure

```
guide-assignment-pipeline/
├── main.nf                    # Main pipeline with method routing
├── nextflow.config           # Configuration and parameters
├── bin/
│   └── crispat.py           # Standalone Python scripts
├── modules/
│   └── crispat/
│       ├── crispat.nf       # Nextflow process definition
│       ├── environment.yml  # Conda environment (pinned versions)
│       └── venv/           # Local Python environment (gitignored)
├── configs/
│   └── config.csv          # Defines method/dataset combinations
└── conda-cache/            # Cached conda packages (gitignored)
```

## Data Flow

1. **Configuration**: `config.csv` defines method/dataset/resource combinations
2. **Input**: H5AD files from `dataset_base_dir/{dataset_id}/`
3. **Method Routing**: Pipeline branches to appropriate method based on CSV
4. **Processing**: Each method runs via `bin/{method}.py` script
5. **Output**: Results organized by run ID in `out_base_dir/{run_id}/{dataset_id}/{method}/`

## Usage

```bash
# Run pipeline with specific run identifier
nextflow run main.nf --run_id "experiment_name"

# Development: Test Python script directly
cd modules/crispat
source venv/bin/activate
python ../../bin/crispat.py input.h5ad output/ --dataset-id test
```

## Configuration

**Run-specific outputs** via `--run_id` parameter:
```
outputs/
├── baseline_experiment/
├── parameter_test_1/
└── final_analysis/
```

**Method/dataset combinations** in `configs/config.csv`:
```csv
method,dataset,cpus,memory
crispat,example,1,4.GB
```

## Environment Management

- **Module-specific environments**: Each module in `modules/{method}/` has its own venv for isolated dependencies
- **Testing modules**: Always use the module's venv when testing individual methods
- **Project-wide environment**: Use the overall project venv (e.g., sceptre3-env) for scripts outside the pipeline

## Development Pattern

1. **Prototype in Python**: Use local venv in `modules/{method}/venv/`
2. **Create standalone script**: Add to `bin/{method}.py` with CLI
3. **Nextflow integration**: Reference script from `modules/{method}/{method}.nf`
4. **Test pipeline**: Run with `--run_id` for isolated outputs