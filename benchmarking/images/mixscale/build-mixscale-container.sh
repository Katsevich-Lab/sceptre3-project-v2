#!/bin/bash
# Build script for Mixscale Singularity container
#
# Usage:
#   ./build-mixscale-container.sh
#
# Notes:
# - Requires sudo/root for building with singularity/apptainer
# - This will create mixscale.sif in the current directory
# - Building may take 15-30 minutes depending on network speed

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/mixscale.def"
SIF_FILE="$SCRIPT_DIR/mixscale.sif"

echo "Building Mixscale Singularity container..."
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
$CONTAINER_CMD run "$SIF_FILE" R --quiet -e 'suppressPackageStartupMessages({library(Mixscale); library(Seurat)}); sv <- as.character(packageVersion("Seurat")); so <- as.character(packageVersion("SeuratObject")); if (sv != "5.3.1" || so != "5.2.0") stop("expected Seurat 5.3.1 + SeuratObject 5.2.0, got ", sv, " + ", so); cat(R.version.string, "| Mixscale", as.character(packageVersion("Mixscale")), "| Seurat", sv, "| SeuratObject", so, "| presto", as.character(packageVersion("presto")), "OK\n")'

echo ""
echo "Done! You can now use this container in your Nextflow pipeline."
