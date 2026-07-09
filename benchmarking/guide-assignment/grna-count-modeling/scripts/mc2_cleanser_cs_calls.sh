#!/bin/zsh
# Re-run CLEANSER with the CORRECT chemistry (--cs, CROP-seq/Poisson) on the CROP-seq datasets that
# are in / near the figure. Writes to a SEPARATE dir so the --dc calls are preserved for comparison.
#
# ⚠️ CAVEAT (see CLAUDE.md §Open threads): this covers ONLY a549 + gastric_organoid — the only CROP-seq
# datasets in the current mc2 concordance/guide-hist figures. run_cleanser_batch.sh still hardcodes --dc,
# and the mc2_ scripts prefer cleanser_calls_cs/ over cleanser_calls/. So any from-scratch regen, or
# adding another CROP-seq dataset (gasperini, barnyard_lrb100 x2) to those figures, must extend the loop
# below first — otherwise that dataset falls back to --dc and its CROP-seq CLEANSER ceilings inflate.
# Proper fix: generalize run_cleanser_batch.sh to choose --cs/--dc per dataset from the registry modality.
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
