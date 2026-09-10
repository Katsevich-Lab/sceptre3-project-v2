#!/bin/bash
# Build script for the fishash Singularity container
#
# Usage:
#   ./build-fishash-container.sh
#
# Notes:
# - Requires sudo/root for building with singularity/apptainer
# - This will create fishash.sif in the current directory
# - Building may take 15-30 minutes depending on network speed

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/fishash.def"
SIF_FILE="$SCRIPT_DIR/fishash.sif"

echo "Building fishash Singularity container..."
echo "Definition file: $DEF_FILE"
echo "Output file: $SIF_FILE"
echo ""

# Check if apptainer or singularity is available
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

# Build the container
sudo $CONTAINER_CMD build "$SIF_FILE" "$DEF_FILE"

echo ""
echo "Build complete!"
echo "Container saved to: $SIF_FILE"
echo ""
echo "Testing container..."
$CONTAINER_CMD run "$SIF_FILE" R --quiet -e 'suppressPackageStartupMessages({library(fishash); library(SingleCellExperiment)}); v <- as.character(packageVersion("fishash")); if (v != "0.99.3") stop("expected fishash 0.99.3, got ", v); data(tapseq_diffex); res <- fishash(counts(altExp(tapseq_diffex))); cat(R.version.string, "| fishash", v, "| example assigned", sum(assay(res, "assigned")), "pairs OK\n")'

echo ""
echo "Done! You can now use this container in your Nextflow pipeline."
