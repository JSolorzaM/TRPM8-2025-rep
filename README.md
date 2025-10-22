# TRPM8-2025-rep
**Author:** Jocelyn Solorza  
**Affiliation:** PhD(c) DoMoSQB – University of Talca, Chile  
**Repository version:** 1.0 (2025)
# TRPM8 Multiscale Framework: Icilin–Calcium Binding Dynamics and Electronic Interaction Analysis 
# Abstract:

The transient receptor potential melastatin 8 (TRPM8) channel has been established as the principal molecular sensor of cold in mammals and has emerged as a promising pharmacological target for pain therapy. Whereas menthol activates TRPM8 independently of intracellular Ca<sup>2+</sup>, the synthetic agonist icilin strictly requires this ion as a binding cofactor. The mechanistic basis of this ion-dependent agonism remains incompletely understood. In this work, a multiscale computational framework integrating comparative modeling, molecular dynamics simulations (MD), binding free energy decomposition (MMPBSA), and hybrid QM/MM calculations was applied to elucidate the structural and electronic determinants of icilin recognition in the human TRPM8 channel. The results indicate that the Ca<sup>2+</sup> ion does not substantially alter overall binding affinity but reorganizes the hydrogen-bonding and coordination network, incorporating icilin into its coordination sphere while redistributing electron density among critical residues. This reorganization stabilizes coupling between the voltage sensor-like domain (VSLD) and the TRP domain, supporting an allosteric model in which Ca<sup>2+</sup> functions as an enhancer of conformational communication rather than as a simple affinity booster. Non-covalent interaction (NCI) and natural bond orbital (NBO) analyses identified residues such as E782, D802, and E1004 as key sites coordinating Ca<sup>2+</sup> within the TRPM8 cavity. Overall, these findings advance the mechanistic understanding of TRPM8 activation and emphasize the role of ion-mediated polarization. Importantly, although Ca<sup>2+</sup> modulation is a conserved feature across several thermoTRPs, the underlying structural mechanisms are specific to each channel. Beyond fundamental insights, this framework provides a mechanistic basis for the rational design of selective TRPM8 modulators. 

## Overview
This repository contains all computational workflows, scripts, and data associated with the multiscale analysis of the **TRPM8 ion channel**, including structural modeling, molecular dynamics (MD) analyses, hydrogen-bond and contact network evaluation, and hybrid QM/MM calculations.

Each module is self-contained with its own data, outputs, and documentation.

---

##  Repository Structure

###  `Contacts-SaltBridge-Analysis/`
Analysis of **residue–residue contact frequencies** and **salt-bridge formation** across TRPM8 systems (APO, APO–Ca²⁺, Icilin, Icilin–Ca²⁺).

- `data/` – Input topology and trajectory files (.prmtop, .nc).  
- `outputs/` – Generated CSVs and heatmaps (`Contacts_Comparison_TRPM8.csv`, `SaltBridges_Comparison_TRPM8.csv`).  
- `script/` – Python script (`contacts_sb.py`) for contact and salt-bridge frequency computation.  
- `README.md` – Detailed description of methods and execution instructions.

---

### `HBond-Analysis/`
Quantitative analysis of **hydrogen bond networks** from TRPM8 simulations.

- `data/` – Input files and trajectory selections used for hbond analysis.  
- `outputs/` – CSV tables and occupancy summaries for hydrogen bonds per system.  
- `script/` – Python and cpptraj scripts used to compute hydrogen-bond persistence.  
- `README.md` – Workflow and usage documentation.

---

### `QM-MM/`
Hybrid **Quantum Mechanics / Molecular Mechanics (QM/MM)** calculations focused on the TRPM8 binding pocket.

- `Electronic-energy-estimation/` – Energy decomposition and convergence data from QM/MM runs.  
- `NBO/` – Natural Bond Orbital (NBO) charge-transfer and orbital interaction analysis.  
- `NCI/` – Non-Covalent Interaction (NCI) surfaces and interaction maps.

---

###  `Structures_TRPM8/`
All structural models of the TRPM8 systems.

- `apo/` – TRPM8 in apo state.  
- `apo_ca/` – Apo system with bound Ca²⁺.  
- `icilin/` – TRPM8 with Icilin agonist bound.  
- `icilin_ca/` – Icilin + Ca²⁺ bound complex.  

Each folder contains the starting `.pdb` files used for MD simulations.

---

## Purpose
The aim of this repository is to provide a **reproducible and transparent computational framework** for understanding the structural determinants of TRPM8 activation and allosteric modulation through multiscale modeling.

---

## Requirements
All workflows rely on the following main software:
- **AMBER** (MD simulations)
- **MDAnalysis** (trajectory analysis)
- **cpptraj**
- **Gaussian** (QM/MM calculations)
- **Python 3.10+**

---

## Citation
If you use this repository, please cite:

> Solorza, J. *TRPM8 Multiscale Framework: Icilin–Calcium Binding Dynamics and Electronic Interaction Analysis*  
> University of Talca, 2025.

---

## Contact
For questions or collaborations:
**Jocelyn Solorza**  
jsolorza11@alumnos.utalca.cl - jocesolorza@gmail.com  


