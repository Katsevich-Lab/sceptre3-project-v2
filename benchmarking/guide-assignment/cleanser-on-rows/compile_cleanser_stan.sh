#!/bin/bash
#$ -N compile-cleanser-stan
#$ -cwd
#$ -j y
#$ -l m_mem_free=4G

source ~/.bashrc
conda activate cleanser-env

python - <<'PY'
from importlib.resources import files
from cmdstanpy import CmdStanModel
stan_file = files("cleanser").joinpath("dc-guide-mixture.stan")
m = CmdStanModel(stan_file=str(stan_file))  # compiles if needed
print("compiled exe:", m.exe_file)
PY

