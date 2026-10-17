#!/usr/bin/env python3
"""
Generate a purely synthetic, simulated vertex overlap dataset for Brain_Spatial_Overlap_Plot.
Contains NO clinical or participant-level data.
"""

import gzip
import numpy as np
import pandas as pd

def generate_simulated_dataset(output_path="example_overlap_data.tsv.gz", n_cortical=10000, n_subcortical=4000, seed=42):
    np.random.seed(seed)
    records = []
    
    # 1. Cortical vertices (Left & Right; Thickness & Area)
    for hemi in ["lh", "rh"]:
        for moda in ["thickness", "area"]:
            for i in range(1, n_cortical + 1):
                probe_id = f"v_{hemi}_{moda[:2]}_{i:05d}"
                records.append({
                    "Probe": probe_id,
                    "hemi": hemi,
                    "moda": moda,
                    "ROIlabel": "Cortex",
                    "structure_type": "cortical"
                })
                
    # 2. Subcortical vertices (Left & Right; Local thickness & LogJacs)
    sub_rois = ["Thalamus", "Caudate", "Putamen", "Pallidum", "Hippocampus", "Amygdala", "Accumbens"]
    for hemi in ["lh", "rh"]:
        for moda in ["thick", "LogJacs"]:
            for i in range(1, n_subcortical + 1):
                probe_id = f"sub_{hemi}_{moda[:2]}_{i:05d}"
                roi = sub_rois[(i - 1) % len(sub_rois)]
                records.append({
                    "Probe": probe_id,
                    "hemi": hemi,
                    "moda": moda,
                    "ROIlabel": roi,
                    "structure_type": "subcortical"
                })
                
    df = pd.DataFrame(records)
    total_n = len(df)
    
    # Simulate focal cluster hotspots of overlap
    # 92% vertices have 0 overlap (background)
    # 8% vertices have 1 to 4 overlapping traits
    overlap_counts = np.zeros(total_n, dtype=int)
    
    # Create focal simulated clusters
    n_hotspots = int(total_n * 0.08)
    hotspot_indices = np.random.choice(total_n, size=n_hotspots, replace=False)
    
    # Assign simulated trait counts: 1 (55%), 2 (25%), 3 (15%), 4 (5%)
    weights = [0.55, 0.25, 0.15, 0.05]
    overlap_counts[hotspot_indices] = np.random.choice([1, 2, 3, 4], size=n_hotspots, p=weights)
    df["overlap_count"] = overlap_counts
    
    # Simulate continuous test statistic (e.g. z-score / local rGM)
    z_scores = np.random.normal(0, 0.6, size=total_n)
    z_scores[hotspot_indices] = np.random.normal(2.2, 0.7, size=n_hotspots)
    df["z_score"] = np.round(z_scores, 4)
    
    # Binary significance flag
    df["significant"] = df["overlap_count"] > 0
    
    df.to_csv(output_path, sep="\t", index=False, compression="gzip")
    print(f"Generated synthetic simulated dataset: {output_path} ({total_n} vertices)")
    return df

if __name__ == "__main__":
    generate_simulated_dataset()
