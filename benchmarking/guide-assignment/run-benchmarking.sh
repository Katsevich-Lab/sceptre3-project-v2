#!/bin/bash

# Nextflow invocation script for guide assignment benchmarking pipeline

# Set run identifier (required)
RUN_ID="example"

# Run the pipeline
nextflow run guide-assignment-pipeline/main.nf \
    --run_id "$RUN_ID"