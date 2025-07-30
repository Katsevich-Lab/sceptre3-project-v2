#!/usr/bin/env python3
import sys
import pandas as pd
import subprocess
import os

input_mtx = sys.argv[1]
output_dir = "cleanser_output/"
os.makedirs(output_dir, exist_ok=True)

# Run CLEANSER guide assignment
subprocess.run([
    "cleanser", "-i", input_mtx, "-o", f"{output_dir}/posteriors.csv", "--dc"
], check=True)

# Process CLEANSER output to standardized format
df = pd.read_csv(f"{output_dir}/posteriors.csv", skiprows=1, sep='\t', 
                 names=['grna_id', 'cell_id', 'posterior'])

# Assign each cell to gRNA with highest posterior probability
assignments = df.loc[df.groupby('cell_id')['posterior'].idxmax()]

# Write standardized output
pd.DataFrame({
    'cell_id': assignments['cell_id'],
    'grna_id': assignments['grna_id']
}).to_csv("assignments_cleanser.csv", index=False)