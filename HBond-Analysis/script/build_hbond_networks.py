#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_hbond_networks.py

Author: Jocelyn Solorza
Assisted by: ChatGPT (OpenAI)

Description:
------------
Unified pipeline to generate hydrogen-bond interaction networks from a master
CSV produced by cpptraj post-processing. The script filters to functional
residue ranges, encodes interaction type by color/line style, and scales edge
widths by occupancy in the treated condition.

Input (master CSV):
-------------------
"Hbond_comparaciones_offset557_filteredBB.csv" containing at least:
    ['pair', 'Comparison', 'Status', 'APO', 'APO_Ca', 'Icilin', 'Icilin_Ca']

Output:
-------
- One SVG network per comparison (color/style/width encodings)
- One edges_*.csv and nodes_*.csv per comparison (for reproducibility)

Install:
--------
    pip install pandas matplotlib networkx

Example:
--------
    python build_hbond_networks.py \
        -i Hbond_comparaciones_offset557_filteredBB.csv \
        -o out/svg_networks

Notes:
------
- Nodes are colored by structural domain: S1 (cyan), S2 (green),
  S2–S3 (yellow), S3 (dark green), S4 (blue), S4–S5 (orange), TRP (red).
- Edge color & style by Status:
    Gained     → green,  dashed
    Disrupted  → red,    dashed
    Rearranged → orange, solid
    Maintained → gray,   solid
- Edge width scales with treated occupancy:
    width = 1.0 + 5.0 * occupancy
"""

import os
import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# ==========================================================
# Defaults
# ==========================================================

DEFAULT_INPUT = "Hbond_comparaciones_offset557_filteredBB.csv"

# Comparisons (reference vs treated)
DEFAULT_COMPARISONS = [
    "APO vs APO_Ca",
    "APO vs Icilin",
    "APO vs Icilin_Ca",
    "APO_Ca vs Icilin",
    "APO_Ca vs Icilin_Ca",
    "Icilin vs Icilin_Ca",
]

# Map comparison → treated column (for width scaling)
DEST_COL = {
    "APO vs APO_Ca": "APO_Ca",
    "APO vs Icilin": "Icilin",
    "APO vs Icilin_Ca": "Icilin_Ca",
    "APO_Ca vs Icilin": "Icilin",
    "APO_Ca vs Icilin_Ca": "Icilin_Ca",
    "Icilin vs Icilin_Ca": "Icilin_Ca",
}

# Functional residue ranges (both residues must fall in any of these)
DEFAULT_FUNC_RANGES = [(740, 757), (770, 812), (835, 865), (990, 1008)]

# Domain → color (node colors)
DOMAIN_RANGES = [
    ("S1",      (740, 757),  "#00B5FF"),  # cyan
    ("S2",      (770, 785),  "#16A34A"),  # green
    ("S2-S3",   (786, 798),  "#FACC15"),  # yellow
    ("S3",      (799, 820),  "#065F46"),  # dark green
    ("S4",      (821, 845),  "#1D4ED8"),  # blue
    ("S4-S5",   (846, 865),  "#FB923C"),  # orange
    ("TRP",     (990, 1008), "#EF4444"),  # red
]

# Edge color by Status
EDGE_COLOR = {
    "Gained": "green",
    "Disrupted": "red",
    "Rearranged": "orange",
    "Maintained": "gray",
}

# Edge line style by Status
EDGE_STYLE = {
    "Gained": "dashed",
    "Disrupted": "dashed",
    "Rearranged": "solid",
    "Maintained": "solid",
}

# Width scaling: width = base + scale * treated_occupancy
WIDTH_BASE  = 1.0
WIDTH_SCALE = 5.0  # 0.20→2.0 ; 0.50→3.5 ; 0.70→4.5 ; 1.0→6.0


# ==========================================================
# Helpers
# ==========================================================

def parse_args():
    ap = argparse.ArgumentParser(description="Build H-bond networks from a master CSV.")
    ap.add_argument("-i", "--input", default=DEFAULT_INPUT,
                    help="Path to master CSV (default: %(default)s)")
    ap.add_argument("-o", "--outdir", default="svg_networks",
                    help="Output directory for SVGs and CSVs (default: %(default)s)")
    ap.add_argument("--comparisons", nargs="*", default=DEFAULT_COMPARISONS,
                    help="List of comparisons 'Ref vs Treated'")
    ap.add_argument("--ranges", nargs="*", type=str, default=[],
                    help="Functional ranges 'a-b' separated by spaces (default: predefined)")
    ap.add_argument("--figsize", type=str, default="12x10",
                    help="Figure size in inches, e.g. '12x10'")
    return ap.parse_args()


def build_ranges(range_strs):
    if not range_strs:
        return DEFAULT_FUNC_RANGES
    out = []
    for s in range_strs:
        a, b = s.split("-")
        out.append((int(a), int(b)))
    return out


def in_any_range(resi, ranges):
    return any(lo <= resi <= hi for (lo, hi) in ranges)


def pair_both_in_ranges(pair_str, ranges):
    try:
        r1, r2 = map(int, pair_str.split("-"))
        return in_any_range(r1, ranges) and in_any_range(r2, ranges)
    except Exception:
        return False


def domain_of(resi):
    for name, (lo, hi), _color in DOMAIN_RANGES:
        if lo <= resi <= hi:
            return name
    return "Other"


def domain_color(resi):
    for _name, (lo, hi), color in DOMAIN_RANGES:
        if lo <= resi <= hi:
            return color
    return "#93C5FD"  # default light blue for undefined residues


def compute_edge_width(occ):
    try:
        return round(WIDTH_BASE + WIDTH_SCALE * float(occ), 2)
    except Exception:
        return WIDTH_BASE


def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)


def save_nodes_edges_csv(outdir, comp, edges_df, nodes_df):
    tag = comp.replace(" ", "_")
    edges_path = os.path.join(outdir, f"edges_{tag}.csv")
    nodes_path = os.path.join(outdir, f"nodes_{tag}.csv")
    edges_df.to_csv(edges_path, index=False)
    nodes_df.to_csv(nodes_path, index=False)
    print(f"   ↳ CSV exported:\n      - {edges_path}\n      - {nodes_path}")


def draw_and_save_graph(G, pos, edge_groups, node_df, title, out_svg):
    plt.figure(figsize=FIGSIZE)
    node_colors = [node_df.loc[n, "node_color"] if n in node_df.index else "#93C5FD" for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=600, edgecolors='black', linewidths=1)
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold')

    # Draw edges grouped by style (so we can mix dashed and solid)
    for _, group in edge_groups.items():
        if not group["edges"]:
            continue
        nx.draw_networkx_edges(
            G, pos,
            edgelist=group["edges"],
            edge_color=group["colors"],
            width=group["widths"],
            style=group["style"]
        )

    # Manual legend by Status/color (line style implied)
    for status, color in EDGE_COLOR.items():
        plt.plot([], [], color=color, label=status, linewidth=3,
                 linestyle='--' if EDGE_STYLE[status] == 'dashed' else 'solid')

    plt.legend(loc="upper left")
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_svg, format="svg")
    plt.close()
    print(f"   ↳ SVG exported: {out_svg}")


# ==========================================================
# Main
# ==========================================================

if __name__ == "__main__":
    args = parse_args()
    FUNC_RANGES = build_ranges(args.ranges)
    ensure_outdir(args.outdir)

    # Figure size parsing
    try:
        w, h = args.figsize.lower().split("x")
        FIGSIZE = (float(w), float(h))
    except Exception:
        FIGSIZE = (12, 10)

    # Load master CSV
    if not os.path.isfile(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    df = pd.read_csv(args.input)

    # Required columns (now in English)
    required_cols = {"pair", "Comparison", "Status"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns {required_cols}. Found: {df.columns.tolist()}")

    # Optional: occupancy columns (for width scaling)
    occ_cols = [c for c in ["APO", "APO_Ca", "Icilin", "Icilin_Ca"] if c in df.columns]
    if len(occ_cols) < 2:
        print("Warning: insufficient occupancy columns detected. Proceeding with limited width scaling.")

    # Filter to functional residue ranges (both residues must be inside)
    df_func = df[df["pair"].apply(lambda p: pair_both_in_ranges(p, FUNC_RANGES))].copy()
    if df_func.empty:
        print("No residue pairs fall within the functional ranges provided. Check your ranges.")
    else:
        print(f"Functional pairs selected: {len(df_func)}")

    # Process each comparison
    for comp in args.comparisons:
        if comp not in df_func["Comparison"].unique():
            print(f"• {comp}: no rows for this comparison after filtering.")
            continue

        treated_col = DEST_COL.get(comp, None)
        if treated_col is None or treated_col not in df_func.columns:
            print(f"• {comp}: treated occupancy column '{treated_col}' not found. Using base edge width.")

        sub = df_func[df_func["Comparison"] == comp].copy()
        if sub.empty:
            print(f"• {comp}: no edges after functional filtering.")
            continue

        # Build graph
        G = nx.Graph()
        edges_records = []
        nodes_records = set()

        for _, row in sub.iterrows():
            try:
                r1, r2 = map(int, row["pair"].split("-"))
            except Exception:
                continue

            status = str(row["Status"]) if pd.notna(row["Status"]) else "Rearranged"
            color = EDGE_COLOR.get(status, "black")
            style = EDGE_STYLE.get(status, "solid")
            occ_treated = row.get(treated_col, float("nan"))
            width = compute_edge_width(occ_treated) if pd.notna(occ_treated) else WIDTH_BASE

            G.add_node(r1)
            G.add_node(r2)
            G.add_edge(r1, r2, status=status, color=color, style=style,
                       width=width, treated_occupancy=occ_treated)

            edges_records.append({
                "source": r1,
                "target": r2,
                "Comparison": comp,
                "Status": status,
                "Color": color,
                "Style": style,
                "Treated_Occupancy": round(float(occ_treated), 4) if pd.notna(occ_treated) else None,
                "Line_Width": width
            })
            nodes_records.update([r1, r2])

        if G.number_of_edges() == 0:
            print(f"• {comp}: graph has no edges.")
            continue

        # Node table with domain & color
        nodes_df = pd.DataFrame({"node": sorted(nodes_records)})
        nodes_df["domain"] = nodes_df["node"].apply(domain_of)
        nodes_df["node_color"] = nodes_df["node"].apply(domain_color)
        nodes_df.set_index("node", inplace=True)

        edges_df = pd.DataFrame(edges_records)

        # Degree-based shell layout
        high_degree = [n for n, d in G.degree() if d > 3]
        medium_degree = [n for n, d in G.degree() if 2 <= d <= 3]
        low_degree = [n for n, d in G.degree() if d == 1]
        shells = [lst for lst in [high_degree, medium_degree, low_degree] if lst]
        pos = nx.shell_layout(G, nlist=shells) if shells else nx.spring_layout(G, seed=42)

        # Group edges by style
        solid_edges, dashed_edges = [], []
        solid_colors, dashed_colors = [], []
        solid_widths, dashed_widths = [], []

        for u, v, data in G.edges(data=True):
            if data.get("style", "solid") == "dashed":
                dashed_edges.append((u, v))
                dashed_colors.append(data.get("color", "black"))
                dashed_widths.append(data.get("width", WIDTH_BASE))
            else:
                solid_edges.append((u, v))
                solid_colors.append(data.get("color", "black"))
                solid_widths.append(data.get("width", WIDTH_BASE))

        edge_groups = {
            "solid":  {"edges": solid_edges,  "colors": solid_colors,  "widths": solid_widths,  "style": "solid"},
            "dashed": {"edges": dashed_edges, "colors": dashed_colors, "widths": dashed_widths, "style": "dashed"},
        }

        tag = comp.replace(" ", "_")
        out_svg = os.path.join(args.outdir, f"Hbond_Network_{tag}.svg")
        draw_and_save_graph(G, pos, edge_groups, nodes_df, f"H-bond Network ({comp})", out_svg)

        # Save edges/nodes tables
        save_nodes_edges_csv(args.outdir, comp, edges_df, nodes_df.reset_index())

    print("\nH-bond network generation completed.")

