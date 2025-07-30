#!/usr/bin/env python3
import sys
import pandas as pd
import subprocess
from pathlib import Path
import os

input_mtx = sys.argv[1]
output_dir = "cleanser_output/"
os.makedirs(output_dir, exist_ok=True)

# Run CLEANSER guide assignment
# Using --dc flag for direct capture (adjust based on your experiment type)
cmd = [
    "cleanser", 
    "-i", input_mtx,
    "-o", f"{output_dir}/posteriors.csv",
    "--dc"  # Use direct capture model - change to --cs for crop-seq if needed
]

try:
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print("CLEANSER stdout:", result.stdout)
except subprocess.CalledProcessError as e:
    print("CLEANSER failed:", e.stderr)
    raise

# Post-process output to standardized format
posteriors_output = Path(output_dir) / "posteriors.csv"

# CLEANSER output format: grna_id cell_id posterior_probability
# First line is header: n_grnas n_cells n_entries
# Skip the first line and read the data
df = pd.read_csv(posteriors_output, skiprows=1, sep='\t', 
                 names=['grna_id', 'cell_id', 'posterior'])

# For each cell, assign to gRNA with highest posterior probability
assignments = df.loc[df.groupby('cell_id')['posterior'].idxmax()]

# Convert to standardized format (cell_id, grna_id)
standardized_df = pd.DataFrame({
    'cell_id': assignments['cell_id'],
    'grna_id': assignments['grna_id']
}).sort_values('cell_id').reset_index(drop=True)

# Write standardized output
standardized_df.to_csv("assignments_cleanser.csv", index=False)