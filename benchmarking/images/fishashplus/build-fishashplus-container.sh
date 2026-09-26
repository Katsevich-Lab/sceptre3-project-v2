#!/bin/bash
# Build script for the fishashplus Apptainer container
#
# Usage:
#   ./build-fishashplus-container.sh [--fakeroot] [path/to/fishashplus clone]
#
#   --fakeroot  build without sudo (apptainer build --fakeroot), e.g. on Betty
#               via ../build-on-betty.sbatch. Without it, the build uses sudo.
#
# Notes:
# - fishashplus is a PRIVATE repo, so the image does not fetch it from GitHub.
#   This script packs the commit pinned in fishashplus.def out of a local clone
#   (default: a `fishashplus` checkout next to this repo) and builds from a
#   temporary directory OUTSIDE this repo, so the private source never lands in
#   it; only fishashplus.sif is written here. No credentials are involved, and
#   the commit is checked in %test.
# - Requires sudo/root for building with apptainer/singularity
# - This will create fishashplus.sif in this directory
# - Building may take 15-30 minutes depending on network speed

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEF_FILE="$SCRIPT_DIR/fishashplus.def"
SIF_FILE="$SCRIPT_DIR/fishashplus.sif"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"

FAKEROOT=0
if [ "${1:-}" = "--fakeroot" ]; then
    FAKEROOT=1
    shift
fi
SRC_REPO="${1:-$REPO_ROOT/../fishashplus}"

# The pin lives in one place: the commit the .def's %test expects.
COMMIT=$(grep -oE '\[ "\$c" = "[0-9a-f]{40}" \]' "$DEF_FILE" | grep -oE '[0-9a-f]{40}')
[ -n "$COMMIT" ] || { echo "ERROR: could not read the pinned commit from $DEF_FILE"; exit 1; }

if ! git -C "$SRC_REPO" cat-file -e "${COMMIT}^{commit}" 2>/dev/null; then
    echo "ERROR: commit $COMMIT not found in $SRC_REPO"
    echo "Pass the path to your fishashplus clone, and make sure it has that commit (git fetch)."
    exit 1
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

# Stage the recipe and the packed source in a throwaway directory; the .def's
# %files paths are relative to it. Removed on exit, even if the build fails.
BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

echo "Packing fishashplus @ $COMMIT from $SRC_REPO"
# git archive exports exactly the committed tree: uncommitted edits in the clone
# cannot leak into the image.
git -C "$SRC_REPO" archive --format=tar.gz --prefix=fishashplus/ \
    -o "$BUILD_DIR/fishashplus-src.tar.gz" "$COMMIT"
echo "$COMMIT" > "$BUILD_DIR/fishashplus-src.commit"
cp "$DEF_FILE" "$BUILD_DIR/"

echo "Building fishashplus Apptainer container..."
echo "Definition file: $DEF_FILE (staged in $BUILD_DIR)"
echo "Output file: $SIF_FILE"
echo "Using container command: $CONTAINER_CMD"
echo ""

cd "$BUILD_DIR"
if [ "$FAKEROOT" -eq 1 ]; then
    $CONTAINER_CMD build --fakeroot "$SIF_FILE" fishashplus.def
else
    sudo $CONTAINER_CMD build "$SIF_FILE" fishashplus.def
fi

echo ""
echo "Build complete!"
echo "Container saved to: $SIF_FILE"
echo ""
echo "Testing container..."
$CONTAINER_CMD test "$SIF_FILE"

echo ""
echo "Done! Copy fishashplus.sif to the pipeline-images directory and add its"
echo "sha256sum to benchmarking/images/SHA256SUMS."
