#!/usr/bin/env python3
import sys
import pandas as pd
import subprocess
import os

if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} <grna_matrix.mtx> <dataset_id>")

input_mtx = sys.argv[1]
dataset_id = sys.argv[2]

output_dir = "cleanser_output/"
os.makedirs(output_dir, exist_ok=True)

# --dc = direct capture, --cs = CROP-seq: the guide capture chemistry. The models
# differ only in the ambient component (native is NB(nbMean * L, nbDisp) in both):
#   --cs  ambient ~ Poisson(lambda)             -- no library-size term
#   --dc  ambient ~ NB(n_nbMean * L, n_nbDisp)  -- scales with the cell
# Since --cs ambient ignores L, its decision boundary is flat in library size,
# which makes it close to a fixed marginal UMI threshold (~6 at CLEANSER's priors).
if "replogle" in dataset_id.lower():
    flag = "--dc"
elif "gasperini" in dataset_id.lower():
    flag = "--cs"
else:
    raise ValueError(f"Unknown dataset '{dataset_id}'. Expected 'replogle' or 'gasperini' in dataset name.")

# CLEANSER IS PARALLEL: it uses 4 cores, unlike crispat and pertpy, which are serial.
# Each guide runs 4 Stan chains at once -- cmdstanpy sets parallel_chains =
# min(cpu_count, chains) and CLEANSER's chain default is 4. That is intrinsic to HMC
# and is left alone. `cpus=4` in the config CSVs DESCRIBES this so the SGE request
# matches; it is not a knob, and nothing here reads it.
#
# -p 1 is the knob. -p is the outer loop over guides -- data parallelism the other
# methods don't use -- and its default is a bare mp.cpu_count() (constants.py:10),
# i.e. the whole node. Pinning it to 1 keeps the footprint at the 4 chain cores.
subprocess.run([
    "cleanser", "-i", input_mtx, "-o", f"{output_dir}/posteriors.csv", flag,
    "-p", "1"
], check=True)

# Process CLEANSER output to standardized format
df = pd.read_csv(f"{output_dir}/posteriors.csv", skiprows=1, sep='\t', 
                 names=['grna_id', 'cell_id', 'posterior'])

df.to_csv("assignments_cleanser.csv", index=False)

# # Assign each cell to gRNA with highest posterior probability
# assignments = df.loc[df.groupby('cell_id')['posterior'].idxmax()]

# # Write standardized output
# pd.DataFrame({
#     'cell_id': assignments['cell_id'],
#     'grna_id': assignments['grna_id']
# }).to_csv("assignments_cleanser.csv", index=False)
