#!/usr/bin/env bash
# Usage: pin_cores.sh N command [args...]
#
# Runs `command` restricted (via taskset) to exactly N of the CPU cores this task
# was allocated, so a method cannot use more cores than intended even when the
# task holds extra CPUs (e.g. requested only to get enough memory). Fails if the
# task has fewer than N cores, rather than silently oversubscribing them.
set -euo pipefail

n=$1; shift

# Cores the scheduler gave this task, e.g. "12-15,140-143"
allowed=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
cores=()
IFS=, read -ra ranges <<< "$allowed"
for r in "${ranges[@]}"; do
  if [[ $r == *-* ]]; then
    for ((c = ${r%-*}; c <= ${r#*-}; c++)); do cores+=("$c"); done
  else
    cores+=("$r")
  fi
done

if (( ${#cores[@]} < n )); then
  echo "pin_cores: need $n cores but this task was allocated ${#cores[@]} ($allowed)" >&2
  exit 1
fi

pick=$(IFS=,; echo "${cores[*]:0:n}")
echo "pin_cores: running on core(s) $pick of allocated $allowed;" \
     "OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset} XLA_FLAGS=${XLA_FLAGS:-unset}" >&2
exec taskset -c "$pick" "$@"
