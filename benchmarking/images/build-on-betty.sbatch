#!/bin/bash
# Build a pipeline image on Betty, without sudo, as a Slurm job.
#
# Usage (from anywhere on the login node):
#   sbatch benchmarking/images/build-on-betty.sbatch <name> [build-script args...]
#   e.g. sbatch benchmarking/images/build-on-betty.sbatch fishashplus
#
# Runs images/<name>/build-<name>-container.sh --fakeroot, which must support that
# flag (so far: fishashplus, grna-simulator). The .sif is written to images/<name>/ as on a laptop
# build; then copy it to the pipeline-images directory and add its sha256sum to
# images/SHA256SUMS.
#
# Betty-specific: the apptainer module names below. Unprivileged --fakeroot builds
# were verified to work on Betty compute nodes (2026-09-26).
#SBATCH -J build-image
#SBATCH -c 4
#SBATCH -t 2:00:00
#SBATCH -o build-image-%j.out

# No `set -u` here: Betty's Lmod scripts read unset variables.
module load arch/zen4 && module load apptainer/1.4.4

set -eo pipefail

if [ $# -lt 1 ]; then
  echo "Usage: sbatch build-on-betty.sbatch <name> [build-script args...]" >&2
  exit 2
fi
NAME=$1; shift

# sbatch runs a copy of this script, so locate images/ via the submit directory's
# repo rather than this file's path.
IMAGES_DIR="$(git -C "$SLURM_SUBMIT_DIR" rev-parse --show-toplevel)/benchmarking/images"
BUILD_SCRIPT="$IMAGES_DIR/$NAME/build-$NAME-container.sh"
[ -x "$BUILD_SCRIPT" ] || { echo "No build script at $BUILD_SCRIPT" >&2; exit 2; }

# Keep the base-image download and build scratch on node-local disk, not $HOME.
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
export APPTAINER_CACHEDIR="$WORK/cache" APPTAINER_TMPDIR="$WORK/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "Building $NAME on $(hostname) with $(apptainer --version), ${SLURM_CPUS_PER_TASK} CPUs"
"$BUILD_SCRIPT" --fakeroot "$@"
