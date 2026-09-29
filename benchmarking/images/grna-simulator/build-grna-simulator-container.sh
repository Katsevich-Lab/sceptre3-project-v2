#!/bin/bash
# Build script for the grna-simulator Apptainer container
#
# Usage:
#   ./build-grna-simulator-container.sh [--fakeroot]
#
#   --fakeroot  build without sudo (apptainer build --fakeroot), e.g. on Betty
#               via ../build-on-betty.sbatch. Without it, the build uses sudo.
#
# Notes:
# - This will create grna-simulator.sif in this directory
# - The build downloads the base image, R packages and basilisk's conda env, so
#   it needs internet access
# - Building may take 30-60 minutes depending on network speed

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/grna-simulator.def"
SIF_FILE="$SCRIPT_DIR/grna-simulator.sif"

FAKEROOT=0
if [ "${1:-}" = "--fakeroot" ]; then
    FAKEROOT=1
    shift
fi

# Check if apptainer or singularity is available
if command -v apptainer &> /dev/null; then
    CONTAINER_CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CONTAINER_CMD="singularity"
else
    echo "ERROR: Neither apptainer nor singularity found in PATH"
    exit 1
fi

echo "Building grna-simulator Apptainer container..."
echo "Definition file: $DEF_FILE"
echo "Output file: $SIF_FILE"
echo "Using container command: $CONTAINER_CMD"
echo ""

# --notest: %test runs once, below, on the finished image, rather than inside the
# build, where a failure discards the whole image. A failing test still fails this
# script (set -e), but leaves the .sif here so a fixed test can be re-run with
# `apptainer test --no-home grna-simulator.sif` instead of rebuilding. Only copy
# the .sif to pipeline-images after the test passes.
if [ "$FAKEROOT" -eq 1 ]; then
    $CONTAINER_CMD build --notest --fakeroot "$SIF_FILE" "$DEF_FILE"
else
    sudo $CONTAINER_CMD build --notest "$SIF_FILE" "$DEF_FILE"
fi

echo ""
echo "Build complete!"
echo "Container saved to: $SIF_FILE"
echo ""
echo "Testing container..."
# --no-home: the test must not be able to write to the real home directory.
$CONTAINER_CMD test --no-home "$SIF_FILE"

echo ""
echo "Done! Copy grna-simulator.sif to the pipeline-images directory and add its"
echo "sha256sum to benchmarking/images/SHA256SUMS."
