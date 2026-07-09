#!/bin/zsh
# Re-run CLEANSER with the CORRECT chemistry (--cs, CROP-seq/Poisson) on the CROP-seq datasets that
# are in / near the figure. Writes to a SEPARATE dir so the --dc calls are preserved for comparison.
set -e
cd /Users/ekatsevi/code/research/sceptre3-project-v2/benchmarking/guide-assignment/grna-count-modeling
BATCH=results/ambient_ceiling/cleanser_batch
OUT=results/ambient_ceiling/cleanser_calls_cs
mkdir -p $OUT
PY=.venv_cleanser/bin/python

for ds in gastric_organoid a549; do
  # export matrix if not present
  if [ ! -f "$BATCH/$ds/m.mtx" ]; then
    Rscript -e "suppressMessages(library(Matrix)); source('scripts/datasets.R');
      m<-load_grna_matrix('$ds'); d<-file.path('$BATCH','$ds'); dir.create(d,recursive=TRUE,showWarnings=FALSE);
      writeMM(m,file.path(d,'m.mtx')); writeLines(rownames(m),file.path(d,'guides.txt')); writeLines(colnames(m),file.path(d,'cells.txt'))"
  fi
  echo "[$(date +%H:%M:%S)] CLEANSER --cs on $ds ..."
  CLNS_MODEL=--cs CLNS_SAMPLES=200 CLNS_WARMUP=100 CLNS_CHAINS=2 CLNS_PARALLEL=6 \
    $PY scripts/run_external_assign.py cleanser \
    $BATCH/$ds/m.mtx $BATCH/$ds/guides.txt $BATCH/$ds/cells.txt $OUT/$ds.csv 300 \
    && echo "[$(date +%H:%M:%S)] $ds --cs DONE ($(wc -l < $OUT/$ds.csv) calls)"
done
echo "cs re-run complete."
