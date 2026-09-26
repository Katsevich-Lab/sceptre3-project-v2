#!/bin/bash
# Build script for the sceptre-manuscript Singularity/Apptainer container.
#
# Usage:
#   ./build-sceptre-manuscript-container.sh
#
# Notes:
# - Requires sudo/root for building with singularity/apptainer.
# - Creates sceptre-manuscript.sif in this directory.
# - Building may take 15-30 minutes (installs tidyverse, bigstatsr, VGAM, sn, ...).

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/sceptre-manuscript.def"
SIF_FILE="$SCRIPT_DIR/sceptre-manuscript.sif"

echo "Building sceptre-manuscript Singularity container..."
echo "Definition file: $DEF_FILE"
echo "Output file: $SIF_FILE"
echo ""

if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "ERROR: Neither apptainer nor singularity found in PATH"
    exit 1
fi

echo "Using container command: $CONTAINER_CMD"
echo ""

sudo $CONTAINER_CMD build "$SIF_FILE" "$DEF_FILE"

echo ""
echo "Build complete!"
echo "Container saved to: $SIF_FILE"
echo ""
echo "Testing container..."
$CONTAINER_CMD run "$SIF_FILE" R --quiet -e 'suppressPackageStartupMessages(library(sceptre)); stopifnot(as.character(packageVersion("sceptre"))=="1.0.0"); stopifnot(exists("run_sceptre_using_precomp", mode="function")); cat("manuscript sceptre", as.character(packageVersion("sceptre")), "OK\n")'

echo ""
echo "Done! Point params.sceptre_manuscript_sif at this .sif in nextflow.config."
