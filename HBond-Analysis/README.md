
# TRPM8 Hydrogen-Bond Network Analysis

### Author
**Jocelyn Solorza**  
_Assisted by ChatGPT (OpenAI)_

---

## Overview
This repository provides a unified and transparent workflow to analyze **hydrogen-bond networks** from molecular dynamics (MD) simulations of the human **TRPM8 tetramer** under multiple experimental conditions:

- **APO**
- **APO–Ca²⁺**
- **Icilin**
- **Icilin–Ca²⁺**

The workflow includes:
1. Extraction and classification of hydrogen bonds from cpptraj output (`.avg.dat` files).  
2. Consolidation of results into a single comparison CSV file.  
3. Construction and visualization of condition-specific hydrogen-bond interaction networks.

The analysis quantifies structural rearrangements in H-bond connectivity across states and generates network visualizations suitable for figures and supplementary materials.

---

##  Methodological Basis
Residue pairs were classified according to occupancy changes between states:

| Category     | Definition |
|---------------|-------------|
| **Gained**     | Absent in the reference condition (occupancy = 0), but present (≥ 0.20) in the treated condition. |
| **Disrupted**  | Present with occupancy ≥ 0.20 in the reference condition, but absent (occupancy = 0) in the treated condition. |
| **Rearranged** | Present in both states (occupancy > 0) with an absolute difference ≥ 0.20. |
| **Maintained** | Present in both states with occupancy ≥ 0.50 and a difference < 0.20. |

The resulting networks encode:
- **Edge color & style:** interaction type  
- **Edge thickness:** treated occupancy (persistence)  
- **Node color:** structural domain (S1–S4, TRP, etc.)

## Acknowledgments
This repository was developed by **Jocelyn Solorza** with technical assistance 
from **ChatGPT (OpenAI)** for documentation and Python optimization.
All scientific design, analysis, and interpretation were performed by the author.

## Dependencies
Install via pip:
```bash
pip install pandas matplotlib networkx
