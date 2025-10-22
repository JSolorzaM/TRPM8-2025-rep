# TRPM8 Contact and Salt-Bridge Analysis

### Author
**Jocelyn Solorza**  
_Assisted by ChatGPT (OpenAI)_

---

##  Description
This repository contains a Python script designed to reproduce the original TRPM8 contact and salt-bridge analysis using legacy parameters.

The script performs all-atom distance-based frequency calculations between selected residues of the TRPM8 ion channel under multiple simulation conditions. It produces both numerical data and visual heatmaps for comparison across systems.

This implementation exactly replicates the behavior of the original workflow:
- **Contacts:** Any-atom minimum distance ≤ 4.5 Å between residues of interest.  
- **Salt bridges:** Any-atom minimum distance ≤ 4.0 Å between positively (ARG/LYS) and negatively (ASP/GLU) charged residues, excluding adjacent ones (|Δresid| < 2).  
- **No PBC correction** (pure Cartesian coordinates).  
- **Fixed 21×21 matrix output** using consistent coloring (`rocket_r` colormap, vmin=0, vmax=1).  

---

## Input Files
The script expects the following molecular dynamics files in the working directory:

| System       | Topology (`.prmtop`) | Trajectory (`.nc`)         |
|---------------|----------------------|-----------------------------|
| APO           | `apo.prmtop`         | `apo.nc`              |
| APO + Ca²⁺    | `apo.prmtop`         | `apo_ca.nc`           |
| Icilin        | `icilin.prmtop`      | `icilin.nc`           |
| Icilin + Ca²⁺ | `icilin_ca.prmtop`   | `icilin_ca.nc`        |

---

##  Output Files
After execution, the following files are generated:

| Type | File Name | Description |
|------|------------|-------------|
| CSV  | `Contacts_Comparison_TRPM8.csv` | Matrix of contact frequencies per residue pair per system |
| CSV  | `SaltBridges_Comparison_TRPM8.csv` | Matrix of salt bridge frequencies |
| SVG  | `heatmap_Contacts_<system>.svg` | Heatmap visualization for contacts |
| SVG  | `heatmap_SaltBridges_<system>.svg` | Heatmap visualization for salt bridges |

All heatmaps are symmetric 21×21 matrices showing residue–residue frequency values.

---
## Acknowledgments
This repository was developed by **Jocelyn Solorza** with technical assistance 
from **ChatGPT (OpenAI)** for documentation and Python optimization.
All scientific design, analysis, and interpretation were performed by the author.

##  How to Run

1. Install required Python packages:
   ```bash
   pip install MDAnalysis numpy pandas matplotlib seaborn
   python contacts_sb.py
