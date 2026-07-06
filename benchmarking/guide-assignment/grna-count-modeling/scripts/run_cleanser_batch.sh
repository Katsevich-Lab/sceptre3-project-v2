#!/bin/zsh
# Run CLEANSER (direct-capture) on each exported small dataset, smallest-first, writing calls to
# results/ambient_ceiling/cleanser_calls/<ds>.csv . Best-effort: each dataset is independent, so
# partial coverage is fine for the writeup. Reduced sampling (n=200,w=100) for speed.
set -e
cd /Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling
BATCH=results/ambient_ceiling/cleanser_batch
OUT=results/ambient_ceiling/cleanser_calls
mkdir -p $OUT
PY=.venv_cleanser/bin/python
while read ds; do
  [ -f "$OUT/$ds.csv" ] && { echo "[$ds] already done, skip"; continue; }
  echo "[$(date +%H:%M:%S)] CLEANSER on $ds ..."
  CLNS_MODEL=--dc CLNS_SAMPLES=200 CLNS_WARMUP=100 CLNS_CHAINS=2 CLNS_PARALLEL=6 \
    $PY scripts/run_external_assign.py cleanser \
    $BATCH/$ds/m.mtx $BATCH/$ds/guides.txt $BATCH/$ds/cells.txt $OUT/$ds.csv 300 \
    > $BATCH/$ds.log 2>&1 && echo "[$(date +%H:%M:%S)] $ds DONE ($(wc -l < $OUT/$ds.csv) calls)" \
    || echo "[$(date +%H:%M:%S)] $ds FAILED (see $BATCH/$ds.log)"
done < $BATCH/order.txt
echo "CLEANSER batch complete."
