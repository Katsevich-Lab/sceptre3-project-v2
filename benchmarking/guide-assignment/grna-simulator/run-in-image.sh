#!/bin/bash
# Run an R script, or an interactive R session, inside grna-simulator.sif on
# Betty, with the home directory READ-ONLY.
#
#   sbatch run-in-image.sh <script.R> [args...]             # batch job
#   sbatch -c 4 -t 6:00:00 run-in-image.sh <script.R> ...   # more CPUs/memory/time
#   srun -p genoa-std-mem -c 2 -t 2:00:00 --pty run-in-image.sh   # interactive R
#
# Run from inside the repo; the script path is resolved from the current
# directory. The log (batch) is grna-sim-<jobid>.out in the directory submitted
# from, and starts with the git commit and any uncommitted changes.
#
# Inside the container:
#   - home is mounted read-only, so the repo, ~/.Rprofile and ~/.research_config
#     can be read but nothing can be written there: a stray cache or config file
#     fails with "Read-only file system" rather than landing in home;
#   - LOCAL_BENCHMARKING_DIR (from ~/.research_config) is mounted read-write, at
#     its resolved path -- a path through a symlink such as ~/data cannot be
#     mounted into a container;
#   - the lab's shared config is mounted read-only: every directory that
#     ~/.Rprofile, ~/.Renviron, ~/.research_config -- and any file
#     ~/.research_config sources -- passes through on its way to the real file.
#     Every hop of a symlink chain must exist inside the container, or the
#     chain breaks (~/.Rprofile -> config/.Rprofile -> lab-resources/...);
#   - --cleanenv: no host environment variables.
# Output must therefore go under LOCAL_BENCHMARKING_DIR (or /tmp, which is
# node-local and discarded).
#
# The image defaults to grna-simulator.sif in the pipeline-images directory;
# set SIF=/path/to/other.sif to use another.
#
#SBATCH -J grna-sim
#SBATCH -p genoa-std-mem
#SBATCH -c 2
#SBATCH -t 2:00:00
#SBATCH -o grna-sim-%j.out

# No `set -u` here: Betty's Lmod scripts read unset variables.
module load arch/zen4 && module load apptainer/1.4.4

set -eo pipefail

SIF=${SIF:-$(readlink -f "$HOME/data/projects/sceptre3/benchmarking/pipeline-images")/grna-simulator.sif}
[ -f "$SIF" ] || { echo "No image at $SIF (set SIF=...)" >&2; exit 2; }

. "$HOME/.research_config"
[ -n "$LOCAL_BENCHMARKING_DIR" ] || { echo "LOCAL_BENCHMARKING_DIR is not set by ~/.research_config" >&2; exit 2; }
DATA=$(readlink -f "$LOCAL_BENCHMARKING_DIR")

# The directory of every hop in a file's symlink chain, the file's own included.
chain_dirs() {
  local p=$1 t
  while [ -e "$p" ] || [ -L "$p" ]; do
    dirname "$p"
    [ -L "$p" ] || break
    t=$(readlink "$p")
    case $t in /*) p=$t ;; *) p=$(dirname "$p")/$t ;; esac
  done
}
STARTUP=("$HOME/.Rprofile" "$HOME/.Renviron" "$HOME/.research_config")
# Files ~/.research_config sources by absolute path (lines like `. /path/file`).
STARTUP+=($(sed -n 's#^[[:space:]]*\(\.\|source\)[[:space:]]\+\(/[^[:space:];]*\).*#\2#p' "$HOME/.research_config"))
CONFIG_DIRS=$(for f in "${STARTUP[@]}"; do chain_dirs "$f"; done |
              xargs -r -n1 readlink -f | sort -u | grep -v "^$HOME\(/\|$\)" || true)

BINDS=(--bind "$HOME:$HOME:ro" --bind "$DATA:$DATA")
for d in $CONFIG_DIRS; do BINDS+=(--bind "$d:$d:ro"); done
RUN=(apptainer exec --cleanenv --no-mount home,cwd "${BINDS[@]}" --pwd "$PWD" "$SIF")

echo "host $(hostname) | job ${SLURM_JOB_ID:-none} | $(date '+%F %T')"
echo "image $SIF"
echo "mounts: home (ro), $DATA (rw)$(for d in $CONFIG_DIRS; do printf ', %s (ro)' "$d"; done)"
echo "git $(git rev-parse HEAD 2>/dev/null || echo 'not in a git repo')"
git status --short 2>/dev/null | sed 's/^/  uncommitted: /' || true
echo

if [ $# -eq 0 ]; then
  exec "${RUN[@]}" R --no-save
fi
SCRIPT=$(readlink -f "$1"); shift
[ -f "$SCRIPT" ] || { echo "No script at $SCRIPT" >&2; exit 2; }
echo "Rscript $SCRIPT $*"
echo
exec "${RUN[@]}" Rscript "$SCRIPT" "$@"
