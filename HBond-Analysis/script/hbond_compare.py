#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hbond_compare.py

Author: Jocelyn Solorza
Assisted by: ChatGPT (OpenAI)

Description:
------------
Processes cpptraj hydrogen bond summary files (.avg.dat) to generate a unified CSV file
describing hydrogen bond occupancies across multiple simulation conditions
(APO, APO–Ca²⁺, Icilin, Icilin–Ca²⁺).

It excludes backbone–backbone interactions (N···O/OXT) and reindexes residues
to a global numbering scheme using a fixed offset. Pairwise comparisons between
conditions are computed to classify each residue pair as:

- "Gained": absent in the reference condition (occupancy = 0) and ≥ 0.20 in the treated condition.
- "Disrupted": present with occupancy ≥ 0.20 in the reference condition, but absent (occupancy = 0) in the treated condition.
- "Rearranged": present in both conditions (occupancy > 0) with an absolute difference ≥ 0.20.
- "Maintained": present in both conditions (occupancy ≥ 0.50) with an absolute difference < 0.20.

Output:
-------
Creates a single CSV file named:
    Hbond_comparaciones_offset557_filteredBB.csv

This CSV is used as input for `build_hbond_networks.py` to generate network visualizations.

Requirements:
-------------
    pip install pandas
"""

import os
import re
import pandas as pd
from itertools import combinations


def parse_cpptraj_avg_filtered_offset(file_path, offset=557):
    """
    Parse a cpptraj .avg.dat hydrogen bond summary file and exclude
    backbone–backbone interactions (donor=N and acceptor in {O, OXT}).

    Expected line format (space-separated, comment lines start with '#'):
        0: Acceptor heavy atom (e.g., PHE_317@O)
        1: DonorH atom       (e.g., MET_321@H)
        2: Donor heavy atom  (e.g., MET_321@N)
        3: Frames
        4: Frac (occupancy)
        5: AvgDist
        6: AvgAng

    Parameters
    ----------
    file_path : str
        Path to the .avg.dat file.
    offset : int, optional
        Integer offset applied to residue indices to map to the global numbering.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by residue pair ("pair"), containing the maximum
        occupancy observed for that pair.
    """
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 5:
                continue  # malformed line

            acceptor_str = parts[0]
            donor_str = parts[2]
            try:
                occupancy = float(parts[4])
            except ValueError:
                continue

            # Extract residue indices and atom names (Amber format: RES_####@ATOM)
            acc_m = re.search(r"_(\d+)@([A-Za-z0-9]+)$", acceptor_str)
            don_m = re.search(r"_(\d+)@([A-Za-z0-9]+)$", donor_str)
            if not (acc_m and don_m):
                continue

            acc_res, acc_atom = acc_m.groups()
            don_res, don_atom = don_m.groups()

            # === Exclude backbone–backbone N···O/OXT ===
            is_backbone_acceptor = acc_atom in {"O", "OXT"}
            is_backbone_donor = don_atom == "N"
            if is_backbone_acceptor and is_backbone_donor:
                continue  # skip backbone-only H-bonds

            # Keep all other combinations (side-chain, mixed, etc.)
            r1 = int(don_res) + offset
            r2 = int(acc_res) + offset
            pair = f"{min(r1, r2)}-{max(r1, r2)}"
            data.append((pair, occupancy))

    df = pd.DataFrame(data, columns=["pair", "occupancy"])
    if df.empty:
        return pd.DataFrame(columns=["occupancy"])
    # Keep maximum occupancy per pair across the file
    return df.groupby("pair").agg({"occupancy": "max"})


# ================================================================
# Main execution block
# ================================================================
if __name__ == "__main__":
    # Input files (edit paths if necessary)
    input_files = {
        "APO": "hbond_APO.avg.dat",
        "APO_Ca": "hbond_APO_Ca.avg.dat",
        "Icilin": "hbond_Icilin.avg.dat",
        "Icilin_Ca": "hbond_Icilin_Ca.avg.dat",
    }

    # If files are in another directory, uncomment and modify:
    # BASE_DIR = "/path/to/your/files"
    # input_files = {k: os.path.join(BASE_DIR, v) for k, v in input_files.items()}

    df_all_offset = []
    missing = []

    # Parse all input datasets
    for label, file in input_files.items():
        if not os.path.isfile(file):
            print(f"File not found: {file}")
            missing.append(file)
            continue
        df = parse_cpptraj_avg_filtered_offset(file, offset=557)
        df.columns = [label]
        df_all_offset.append(df)

    if not df_all_offset:
        raise SystemExit(f"No input files found. Missing: {missing}")

    # Merge on 'pair'
    df_merged_offset = pd.concat(df_all_offset, axis=1).fillna(0.0)

    # Pairwise comparisons between conditions
    comparisons = list(combinations([c for c in input_files.keys() if c in df_merged_offset.columns], 2))
    results_offset = []

    for s1, s2 in comparisons:
        subset = df_merged_offset[[s1, s2]].copy()
        subset["Status"] = ""

        # Classification according to thresholds defined in the Methods
        subset.loc[(subset[s1] == 0) & (subset[s2] >= 0.20), "Status"] = "Gained"
        subset.loc[(subset[s1] >= 0.20) & (subset[s2] == 0), "Status"] = "Disrupted"
        subset.loc[(subset[s1].gt(0) & subset[s2].gt(0) &
                    (subset[s1].sub(subset[s2]).abs() >= 0.20)), "Status"] = "Rearranged"
        subset.loc[(subset[s1].sub(subset[s2]).abs() < 0.20) &
                   (subset[s1] >= 0.50) & (subset[s2] >= 0.50), "Status"] = "Maintained"

        subset = subset[subset["Status"] != ""].copy()
        subset["Comparison"] = f"{s1} vs {s2}"
        results_offset.append(subset.reset_index())

    # Concatenate all comparisons
    final_offset_df = (pd.concat(results_offset, ignore_index=True)
                       if results_offset else
                       pd.DataFrame(columns=["pair"] + list(input_files.keys()) + ["Status", "Comparison"]))

    out_csv = "Hbond_comparaciones_offset557_filteredBB.csv"
    final_offset_df.to_csv(out_csv, index=False)

    # Minimal console report
    print("Hydrogen bond comparison completed successfully.")
    print(f"   → Non-zero pairs per condition: " +
          ", ".join(f"{c}={int((df_merged_offset[c] > 0).sum())}" for c in df_merged_offset.columns))
    print(f"   → Total comparisons generated: {len(comparisons)}")
    print(f"   → Output file saved as: {out_csv}")

