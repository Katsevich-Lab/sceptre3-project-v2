#!/bin/bash

# Nextflow invocation script for guide assignment benchmarking pipeline

# works for conda and mamba
export CONDA_ALWAYS_YES=1

# if you’re actually using micromamba under the hood, this also works:
export MAMBA_NO_CONFIRM=1

# Set run identifier (required)
RUN_ID="gasperini"

# Run the pipeline
nextflow run guide-assignment-pipeline/main.nf \
    --run_id "$RUN_ID"
