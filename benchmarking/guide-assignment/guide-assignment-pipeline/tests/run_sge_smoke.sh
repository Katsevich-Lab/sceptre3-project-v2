#!/bin/bash
#$ -N ga-sge-smoke
#$ -cwd
#$ -j y
#$ -l m_mem_free=2G
#$ -l h_rt=02:00:00      # driver must OUTLIVE the test; without this it lands on
                         # short.q (s_rt=4h) -- fine here, but we set it explicitly
                         # because that is the very behaviour under test
#$ -q hpc3.q             # long-job queue; driver is light, just orchestrates
#$ -pe openmp 1
#
# Graceful-failure smoke test for the guide-assignment pipeline, to be run ON
# WHARTON HPC3:
#
#     cd benchmarking/guide-assignment/guide-assignment-pipeline
#     qsub tests/run_sge_smoke.sh
#
# It answers the questions a laptop cannot:
#   1. Does a task killed by SGE (walltime) leave its siblings running?
#   2. Does `time` in the config actually become `-l h_rt`, and does the queue
#      routing in ~/.nextflow/config then send it to the right queue?
#   3. Does the trace record status / exit / native_id for tasks that died?
#   4. Does the workflow really survive to completion?
#
# Takes roughly 10-15 min: the walltime-killed task needs ~3 min to be killed,
# then Nextflow waits out its exitReadTimeout (~4.5 min default) before recording
# the failure, plus queue wait.

set -uo pipefail   # NOT -e: we want to inspect failures, not abort on them

# --------------------------------------------------------------------------
# Locate ourselves.
#
# DO NOT use "$0" here. SGE COPIES the job script into its spool directory
# before executing it, so $0 is /var/sge/.../job_scripts/<jobid>, NOT this file
# in the repo. Deriving paths from it puts every mkdir into a root-owned spool
# dir and the whole run dies with "Permission denied". "#$ -cwd" runs the job in
# the SUBMISSION directory, so resolve from there instead.
# --------------------------------------------------------------------------
SUBMIT_DIR="${SGE_O_WORKDIR:-$PWD}"

if   [ -f "${SUBMIT_DIR}/sge_smoke.nf" ];       then HERE="${SUBMIT_DIR}"
elif [ -f "${SUBMIT_DIR}/tests/sge_smoke.nf" ]; then HERE="${SUBMIT_DIR}/tests"
else
  echo "SETUP ERROR: cannot find sge_smoke.nf from submission dir '${SUBMIT_DIR}'." >&2
  echo "Submit from the pipeline dir or from tests/:" >&2
  echo "    cd <repo>/benchmarking/guide-assignment/guide-assignment-pipeline" >&2
  echo "    qsub tests/run_sge_smoke.sh" >&2
  exit 2
fi
PIPE_DIR="$(dirname "$HERE")"
TESTOUT="${HERE}/smoke-out"

# Fail LOUDLY and distinctly on setup problems. A setup failure that falls
# through to the verdict block prints "all checks FAILED", which reads as
# "graceful failure is broken" when the run never even started.
[ -f "${PIPE_DIR}/nextflow.config" ] || {
  echo "SETUP ERROR: no nextflow.config at ${PIPE_DIR}" >&2; exit 2; }
command -v nextflow >/dev/null 2>&1 || {
  echo "SETUP ERROR: nextflow not on PATH (module load / conda activate first)" >&2; exit 2; }
mkdir -p "$TESTOUT" 2>/dev/null || {
  echo "SETUP ERROR: cannot write to ${TESTOUT}" >&2; exit 2; }

echo "submission dir : ${SUBMIT_DIR}"
echo "test dir       : ${HERE}"
echo "pipeline dir   : ${PIPE_DIR}"
echo

# Guard the rm -rf. HERE is already validated (it must contain sge_smoke.nf),
# but never let an unexpected expansion reach a recursive delete on a shared box.
[ -n "${HERE}" ] && [ "${TESTOUT}" = "${HERE}/smoke-out" ] || {
  echo "SETUP ERROR: refusing to rm -rf unexpected path '${TESTOUT}'" >&2; exit 2; }
rm -rf "$TESTOUT" && mkdir -p "$TESTOUT" "${HERE}/nf-logs"

export NXF_OPTS="-Xms512m -Xmx2g"
export NXF_HOME="${PIPE_DIR}/.nextflow"

# Load the SAME configs the real pipeline uses, so this tests the real settings
# (errorStrategy, clusterOptions, queue routing) rather than a copy of them.
nextflow \
  -C ~/.nextflow/config \
  -C "${PIPE_DIR}/nextflow.config" \
  run "${HERE}/sge_smoke.nf" \
  --testout "$TESTOUT" \
  -with-trace "${TESTOUT}/trace.tsv"
NF_EXIT=$?

echo
echo "################################################################"
echo "# VERDICT"
echo "################################################################"

fail=0
note () { printf '  %-6s %s\n' "$1" "$2"; }

# --- 1. workflow itself survived ------------------------------------------
if [ "$NF_EXIT" -eq 0 ]; then
  note "PASS" "workflow exited 0 despite 3 tasks dying"
else
  note "FAIL" "workflow exited ${NF_EXIT} -- a dying task took the run down"; fail=1
fi

# --- 2. survivors actually produced output --------------------------------
for f in out_survivor1.txt out_survivor2.txt; do
  if [ -s "${TESTOUT}/${f}" ]; then note "PASS" "survivor published ${f}"
  else note "FAIL" "MISSING ${f} -- a healthy task was killed alongside its siblings"; fail=1; fi
done

# --- 3. the deaths were recorded, not silently swallowed ------------------
TR="${TESTOUT}/trace.tsv"
# Column lookup BY NAME -- trace.fields order is configurable, so never index
# columns positionally.
col () { head -1 "$TR" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }

if [ -s "$TR" ]; then
  C_STATUS=$(col status); C_EXIT=$(col exit); C_NATIVE=$(col native_id)
  nfailed=$(awk -F'\t' -v c="$C_STATUS" 'NR>1 && $c=="FAILED"' "$TR" | wc -l)
  if [ "$nfailed" -ge 3 ]; then note "PASS" "trace recorded ${nfailed} FAILED tasks"
  else note "FAIL" "trace shows only ${nfailed} FAILED (expected 3)"; fail=1; fi

  if [ -n "$C_STATUS" ] && [ -n "$C_EXIT" ]; then
    note "PASS" "trace carries status + exit columns"
    echo "         exit codes recorded for the dead tasks:"
    awk -F'\t' -v cs="$C_STATUS" -v ce="$C_EXIT" -v ct=4 \
        'NR>1 && $cs=="FAILED" {printf "           %-12s exit=%s\n", $ct, $ce}' "$TR"
    echo "         (a task SGE kills may show a sentinel exit, not 137/143, if it"
    echo "          died before the wrapper could write its .exitcode file)"
    # 124 == coreutils `timeout` fired. This is the mechanism the real modules
    # rely on to bound runtime, since SGE does not enforce h_rt on this cluster.
    if awk -F'\t' -v ce="$C_EXIT" 'NR>1 && $ce==124 {found=1} END{exit !found}' "$TR"; then
      note "PASS" "in-band timeout fired and surfaced as exit 124"
    else
      note "FAIL" "no exit 124 -- the in-band timeout did NOT fire; runtime is unbounded"; fail=1
    fi
  else
    note "FAIL" "trace is missing status/exit -- you cannot tell what failed"; fail=1
  fi
else
  note "FAIL" "no trace file at ${TR}"; fail=1
fi

# --- 4. did `time` become -l h_rt, and where did each task get routed? ----
echo
echo "--- SGE directives actually emitted (proves time -> h_rt + queue routing) ---"
find "${HERE}/work" -name .command.run 2>/dev/null | while read -r f; do
  nm=$(grep -m1 '^#\$ -N' "$f" | sed 's/^#\$ -N //')
  hrt=$(grep -m1 -o 'h_rt=[^ ,]*' "$f")
  q=$(grep -m1 -o '\-q [^ ]*' "$f")
  printf '  %-42s %-18s %s\n' "${nm:-?}" "${hrt:-NO_H_RT!}" "${q:-<no -q: routed by resource match>}"
done

if ! find "${HERE}/work" -name .command.run 2>/dev/null | xargs grep -l 'h_rt=' >/dev/null 2>&1; then
  note "FAIL" "no task emitted -l h_rt -- the time directive is not reaching SGE"; fail=1
else
  note "PASS" "tasks emitted -l h_rt"
fi

# --- 5. what SGE itself says about the killed jobs ------------------------
echo
echo "--- qacct for the tasks that died (ground truth from SGE) ---"
if [ -s "$TR" ]; then
  awk -F'\t' -v cs="$C_STATUS" -v cn="$C_NATIVE" \
      'NR>1 && $cs=="FAILED" {print $cn}' "$TR" | while read -r jid; do
    [ -z "$jid" ] || [ "$jid" = "-" ] && continue
    echo "  job ${jid}:"
    qacct -j "$jid" 2>/dev/null \
      | grep -iE 'exit_status|failed|maxvmem|ru_wallclock|deleted_by' \
      | sed 's/^/    /' || echo "    (qacct not available yet -- retry in a minute)"
  done
fi

echo
if [ "$fail" -eq 0 ]; then
  echo "  ==> ALL CHECKS PASSED: failures are isolated, recorded, and survivable."
else
  echo "  ==> SOME CHECKS FAILED (see above). Do NOT launch the full run yet."
fi
echo "  artifacts: ${TESTOUT}"
exit "$fail"
