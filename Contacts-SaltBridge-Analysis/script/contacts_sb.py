#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRPM8 – Contact and Salt-Bridge Analysis
Author: Jocelyn Solorza
Assisted by: ChatGPT (OpenAI)

-------------------------------------------------------------------------------
DESCRIPTION
-------------------------------------------------------------------------------
This script analyzes residue–residue contacts and salt bridges in TRPM8
molecular dynamics trajectories. It reproduces the exact behavior of the
original analysis pipeline used in Jocelyn Solorza’s TRPM8 project.

Specifically, it:
  • Calculates residue–residue contact frequencies using all-atom distances
    ≤ 4.5 Å between residues within a predefined set.
  • Calculates salt bridge frequencies between charged residues
    (ARG/LYS ↔ ASP/GLU) using all-atom distances ≤ 4.0 Å, excluding
    adjacent residues (|Δresid| < 2).
  • Generates frequency comparison CSV files for multiple systems.
  • Creates symmetric 21×21 heatmaps (SVG format) showing contact
    and salt-bridge occurrence per system.

The script is intentionally “legacy-exact”:
  – No periodic boundary corrections (no PBC imaging)
  – Uses any-atom distance matrices (not Cα-only)
  – Fixed matrix dimensions and coloring (vmin=0, vmax=1, cmap=rocket_r)

-------------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------------
Run from the command line in a folder containing your topology (.prmtop)
and trajectory (.nc) files for each system:

    python contacts_sb.py

The script expects the following files:
    apo.prmtop          apo.nc
    apo.prmtop          apo_Ca.nc
    icilin.prmtop       icilin.nc
    icilin_ca.prmtop    icilin_Ca.nc

Output files:
    Contacts_Comparison_TRPM8.csv
    SaltBridges_Comparison_TRPM8.csv
    heatmap_Contacts_<system>.svg
    heatmap_SaltBridges_<system>.svg

Dependencies:
    MDAnalysis ≥ 2.0, NumPy, Pandas, Matplotlib, Seaborn
Install with:
    pip install MDAnalysis numpy pandas matplotlib seaborn
-------------------------------------------------------------------------------
"""

import MDAnalysis as mda
import numpy as np
import pandas as pd
import os
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

# === SYSTEMS (same file names) ===
systems = {
    "APO": ("apo.prmtop", "apo.nc"),
    "APO_Ca": ("apo.prmtop", "apo_Ca.nc"),
    "Icilin": ("icilin.prmtop", "icilin.nc"),
    "Icilin_Ca": ("icilin_ca.prmtop", "icilin_Ca.nc"),
}

# === RESIDUES (original numbering) ===
resids_original = [
    738, 741, 742, 745, 778, 781, 782, 784, 785,
    796, 799, 802, 805, 839, 842, 845, 846,
    1001, 1004, 1005, 1008
]
OFFSET = 557
resids_internal = [r - OFFSET for r in resids_original]
resids_str = "resid " + " ".join(str(r) for r in resids_internal) + " and protein"

# === CUTOFFS (identical to original logic) ===
CONTACT_CUTOFF = 4.5  # Å, any atom
SALT_CUTOFF = 4.0     # Å, any atom between charged residues

# --- Helper functions ---
def pair_key(i, j):
    return tuple(sorted((i, j)))

def build_df(data_dict, data_type):
    """Build and export comparative CSV identical to the original output."""
    all_pairs_internal = set()
    for d in data_dict.values():
        all_pairs_internal |= set(d.keys())

    all_pairs_orig = sorted((i + OFFSET, j + OFFSET) for (i, j) in all_pairs_internal)
    index = [f"{i}-{j}" for (i, j) in all_pairs_orig]
    df = pd.DataFrame(index=index, columns=list(systems.keys()), dtype=float)
    df[:] = 0.0

    for label in systems:
        for (ri, rj), freq in data_dict[label].items():
            key = f"{min(ri, rj)+OFFSET}-{max(ri, rj)+OFFSET}"
            df.loc[key, label] = freq

    out = f"{data_type}_Comparison_TRPM8.csv"
    df.to_csv(out)
    return df

def plot_heatmap(df, data_type):
    """Generate 21×21 heatmaps identical to the legacy format."""
    labels = sorted(resids_original)
    idx_map = {res: i for i, res in enumerate(labels)}
    for col in df.columns:
        M = np.zeros((21, 21), dtype=float)
        for pair_str, val in df[col].items():
            try:
                r1, r2 = map(int, pair_str.split('-'))
            except Exception:
                continue
            if r1 in idx_map and r2 in idx_map:
                i, j = idx_map[r1], idx_map[r2]
                M[i, j] = val
                M[j, i] = val
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            M, xticklabels=labels, yticklabels=labels,
            cmap='rocket_r', square=True, vmin=0, vmax=1,
            cbar_kws={'label': 'Contact frequency'}
        )
        plt.title(f"{data_type} - {col}", fontsize=14)
        plt.xlabel("Residue", fontsize=12)
        plt.ylabel("Residue", fontsize=12)
        plt.tight_layout()
        plt.savefig(f"heatmap_{data_type}_{col}.svg", format="svg")
        plt.close()

# --- Contact frequencies (legacy) ---
def get_contact_freq_legacy(u, selection, cutoff=4.5):
    """Compute any-atom contact frequency (legacy behavior)."""
    sel = u.select_atoms(selection)
    pairs = list(combinations(sel.residues, 2))
    n_frames = 0
    counts = {}
    for ts in u.trajectory:
        n_frames += 1
        for r1, r2 in pairs:
            if len(r1.atoms) == 0 or len(r2.atoms) == 0:
                continue
            A = r1.atoms.positions
            B = r2.atoms.positions
            d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
            if np.any(d <= cutoff):
                key = pair_key(r1.resid, r2.resid)
                counts[key] = counts.get(key, 0) + 1
    return {k: v / max(1, n_frames) for k, v in counts.items()}

# --- Salt-bridge frequencies (legacy) ---
def get_salt_bridge_freq_legacy(u, resids_int, cutoff=4.0):
    """Compute any-atom salt bridge frequency between charged residues."""
    resid_str = " ".join(str(r) for r in resids_int)
    pos_res = u.select_atoms("resname ARG LYS and resid " + resid_str).residues
    neg_res = u.select_atoms("resname ASP GLU and resid " + resid_str).residues
    n_frames = 0
    counts = {}
    for ts in u.trajectory:
        n_frames += 1
        for r1 in pos_res:
            for r2 in neg_res:
                if abs(r1.resid - r2.resid) < 2:
                    continue
                if len(r1.atoms) == 0 or len(r2.atoms) == 0:
                    continue
                A = r1.atoms.positions
                B = r2.atoms.positions
                d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
                if np.any(d <= cutoff):
                    key = pair_key(r1.resid, r2.resid)
                    counts[key] = counts.get(key, 0) + 1
    return {k: v / max(1, n_frames) for k, v in counts.items()}

# --- Main execution ---
def main():
    contact_results = {}
    salt_results = {}

    for label, (top, traj) in systems.items():
        print(f"\nProcessing {label}...")
        if not (os.path.exists(top) and os.path.exists(traj)):
            print(f"Missing files: {top} / {traj}. Skipping.")
            continue
        try:
            u = mda.Universe(top, traj)
        except ValueError as e:
            if "n_atoms" in str(e) and "!=" in str(e):
                print(f"Atom-count mismatch in {label}, loading without bond guessing...")
                u = mda.Universe(top, traj, guess_bonds=False, topology_guess=False)
            else:
                raise

        contact_results[label] = get_contact_freq_legacy(u, resids_str, cutoff=CONTACT_CUTOFF)
        salt_results[label] = get_salt_bridge_freq_legacy(u, resids_internal, cutoff=SALT_CUTOFF)

    df_contacts = build_df(contact_results, "Contacts")
    df_salts = build_df(salt_results, "SaltBridges")

    plot_heatmap(df_contacts, "Contacts")
    plot_heatmap(df_salts, "SaltBridges")

    print("\nAll done. Files generated:")
    print("- Contacts_Comparison_TRPM8.csv")
    print("- SaltBridges_Comparison_TRPM8.csv")
    print("- heatmap_*.svg (1 per system and type)")

if __name__ == "__main__":
    main()

