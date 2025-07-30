#!/usr/bin/env python3
import sys
import pandas as pd
from crispat import ga_poisson_gauss

input_h5ad = sys.argv[1]
output_dir = "crispat_output/"

# Run CRISPAT guide assignment
ga_poisson_gauss(input_h5ad, output_dir)

# Post-process output to standardized format
crispat_output = f"{output_dir}/assignments.csv"
df = pd.read_csv(crispat_output)

# Convert to standardized format (cell_id, grna_id)
standardized_df = pd.DataFrame({
    'cell_id': df.iloc[:, 0],  # First column should be cell
    'grna_id': df.iloc[:, 1]   # Second column should be gRNA
})

# Write standardized output
standardized_df.to_csv("assignments_crispat.csv", index=False)