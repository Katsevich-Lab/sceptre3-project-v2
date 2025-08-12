#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    if len(sys.argv) != 2:
        print("Usage: python create_fake_sceptre_data.py <output_directory>")
        sys.exit(1)
    
    output_dir = Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Creating fake sceptre data...")
    
    # Parameters for simulation
    num_cells = 500
    num_grnas = 100
    num_NTs = 10
    num_targets = 10
    num_responses = 100
    
    # Create gRNA names
    grna_names = [f"gRNA_{i+1}" for i in range(num_grnas)] + [f"NT_{i+1}" for i in range(num_NTs)]
    cell_names = [f"CELL_{i+1:04d}" for i in range(num_cells)]
    
    # 1. Simulate gRNA matrix (gRNAs x cells) following your R code
    np.random.seed(42)  # For reproducibility
    grna_matrix = np.zeros((num_grnas + num_NTs, num_cells))
    
    # Each cell gets one primary gRNA expressed
    grnas_expressed = np.random.choice(num_grnas + num_NTs, num_cells, replace=True)
    
    for i in range(num_cells):
        # Primary gRNA gets high expression
        grna_matrix[grnas_expressed[i], i] = np.random.uniform(10, 15, 1)
        
        # Other gRNAs get background expression (sparse)
        other_grnas = [j for j in range(num_grnas + num_NTs) if j != grnas_expressed[i]]
        background_counts = np.random.uniform(1,8, len(other_grnas))
        background_mask = np.random.binomial(1, 0.5, len(other_grnas))
        for j, grna_idx in enumerate(other_grnas):
            grna_matrix[grna_idx, i] = background_counts[j] * background_mask[j]
    
    grna_df = pd.DataFrame(
        grna_matrix,
        index=grna_names,  # gRNA names as rows
        columns=cell_names  # cell names as columns
    )
    
    grna_output = output_dir / "grna_matrix.csv"
    grna_df.to_csv(grna_output, index_label="grna_id")
    print(f"✓ Saved gRNA matrix: {grna_output} ({grna_df.shape[0]} gRNAs x {grna_df.shape[1]} cells)")
    
    # 2. Create response matrix (genes x cells)
    response_names = [f"response_{i+1}" for i in range(num_responses)]
    response_matrix = np.random.poisson(3, size=(num_responses, num_cells))
    
    response_df = pd.DataFrame(
        response_matrix,
        index=response_names,
        columns=cell_names
    )
    response_output = output_dir / "response_matrix.csv"
    response_df.to_csv(response_output, index_label="response_id")
    print(f"✓ Saved response matrix: {response_output} ({num_responses} genes x {num_cells} cells)")
    
    # 3. Create gRNA target mapping
    grnas_per_target = num_grnas // num_targets
    targets = []
    
    # Assign targeting gRNAs to targets
    for target_idx in range(num_targets):
        for grna_in_target in range(grnas_per_target):
            targets.append(f"target_{target_idx + 1}")
    
    # Add non-targeting gRNAs
    for nt_idx in range(num_NTs):
        targets.append("non-targeting")
    
    grna_target_df = pd.DataFrame({
        'grna_id': grna_names,
        'grna_target': targets
    })
    target_output = output_dir / "grna_target_data_frame.csv"
    grna_target_df.to_csv(target_output, index=False)
    print(f"✓ Saved gRNA targets: {target_output} ({num_grnas} targeting + {num_NTs} non-targeting)")
    
    print(f"\nCreated fake sceptre data with:")
    print(f"- {num_cells} cells")
    print(f"- {num_grnas} targeting gRNAs ({grnas_per_target} per target)")
    print(f"- {num_NTs} non-targeting gRNAs")
    print(f"- {num_targets} targets")
    print(f"- {num_responses} response genes")

if __name__ == "__main__":
    main()