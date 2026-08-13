#!/usr/bin/env bash
set -euo pipefail

# Run the missing Fishash Table 2 methods on the already-built four-sample
# inputs. The CLEANSER calls follow the Fishash analysis code (both --cs and
# --dc on all samples, package defaults, eight parallel chains). The SCEPTRE
# calls use the matched GEX matrix and version 0.10.3.

HERE=$(cd "$(dirname "$0")" && pwd)
RAW_DIR=${RAW_DIR:?Set RAW_DIR to the Liu GSE272457 directory}
SCEPTRE_LIB=${SCEPTRE_LIB:?Set SCEPTRE_LIB to an R library containing sceptre 0.10.3}
CLEANSER=${CLEANSER:-"$HERE/../../.venv_cleanser/bin/cleanser"}
CPUS=${CPUS:-8}
CLEANSER_SEED=${CLEANSER_SEED:-20260810}

samples=(mix0hr_Cropseq mix72hr_Cropseq mix0hr_DirectCapture mix72hr_DirectCapture)
prefixes=(
  GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_0hr_mix
  GSE272457_293T_LRB100_NTlib1-NIH3T3_LRB100_NTlib2_72hr_mix
  GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_0hr_mix
  GSE272457_293T_MCH2_NTlib1-NIH3T3_MCH2_NTlib2_72hr_mix
)

for idx in 0 1 2 3; do
  sample=${samples[$idx]}
  prefix=${prefixes[$idx]}
  sceptre_out="$HERE/${sample}_sceptre_mixture.mtx"
  if [[ ! -s "$sceptre_out" ]]; then
    R_LIBS_USER="$SCEPTRE_LIB" Rscript "$HERE/run_sceptre_table2.R" \
      --raw-prefix "$RAW_DIR/$prefix" --out "$sceptre_out" --cpus "$CPUS" \
      2>&1 | tee "$HERE/${sample}_sceptre_mixture.log"
  fi

  for mode in cs dc; do
    posterior="$HERE/${sample}_cleanser_${mode}_posterior.mtx"
    log="$HERE/${sample}_cleanser_${mode}.log"
    if [[ ! -s "$posterior" ]] || ! grep -q "Random seed:" "$log" 2>/dev/null; then
      echo "Running CLEANSER ${mode} for ${sample}; detailed sampler output goes to the log."
      "$CLEANSER" -i "$HERE/${sample}_grna_counts.mtx" \
        -o "$posterior" --"$mode" -p "$CPUS" -s "$CLEANSER_SEED" \
        >"$log" 2>&1
    fi
  done
done

Rscript "$HERE/score_extended_table2.R"
