# FAZENA Tools: Modular Drug Discovery Platform

[![Python Version](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.x-orange.svg)](https://pytorch.org/)
[![Streamlit](https://img.shields.io/badge/Streamlit-1.30%2B-red.svg)](https://streamlit.io/)
[![AutoDock Vina](https://img.shields.io/badge/Docking-AutoDock%20Vina-brightgreen.svg)](https://vina.scripps.edu/)
[![RDKit](https://img.shields.io/badge/Cheminformatics-RDKit-green.svg)](https://www.rdkit.org/)
[![License](https://img.shields.io/badge/License-Proprietary%20%2F%20FAZENA-darkgrey.svg)](https://fazena.org)

**FAZENA Tools** is an enterprise-grade bio-computational platform designed for computational chemistry and modern small-molecule drug discovery. The software unifies contrastive deep learning co-embeddings, structure-agnostic RNA pharmacophore modeling, Lipinski and PAINS ADMET safety filtering, and physics-based AutoDock Vina molecular docking inside a responsive, modern web interface.

---

## Table of Contents

- [Overview and Capabilities](#overview-and-capabilities)
- [Core Engines](#core-engines)
  - [RNA Hunter](#1-rna-hunter-rna-target-screening)
  - [Prozena](#2-prozena-protein-target-screening)
  - [Refazenation](#3-refazenation-target-prediction-and-reverse-docking)
- [Quick Start](#quick-start)
- [Installation Guide](#installation-guide)
- [Running the Application](#running-the-application)
- [Repository Structure](#repository-structure)
- [Input Formats and Compatibility](#input-formats-and-compatibility)
- [Automated Test Suite](#automated-test-suite)
- [Authors and Citation](#authors-and-citation)

---

## Overview and Capabilities

FAZENA Tools brings together specialized modules for both oligonucleotide and protein drug discovery:

```
                          +--------------------------+
                          |       FAZENA TOOLS       |
                          +-------------+------------+
         +------------------------------+------------------------------+
         |                              |                              |
         v                              v                              v
+------------------+           +------------------+           +------------------+
|    RNA Hunter    |           |     Prozena      |           |   Refazenation   |
|   (RNA Targets)  |           | (Protein Targets)|           | (Target Predict) |
+------------------+           +------------------+           +------------------+
```

Key technical features:

- **End-to-End Pipeline:** From raw sequence or structure inputs to filtered hit tables, 2D interaction plots, and 3D docked complexes.
- **Deep Learning Co-Embeddings:** ConPLex contrastive neural networks map 2048-bit Morgan fingerprints and 1024-dimensional target representations into a shared latent space for sub-second library evaluation.
- **Physics-Grounding:** Integrated AutoDock Vina engine for automated grid-box calculation, receptor/ligand PDBQT preparation, energetic ranking, and ProLIF residue contact analysis.
- **Pre-Screening Safety Gates:** Built-in Lipinski Rule of 5, Veber parameters, PAINS filter, and toxicophore alert detection.
- **Hardware Agnostic:** Automatic detection of modern NVIDIA GPUs (including NVIDIA Blackwell RTX 50-series and Ada Lovelace architectures) with transparent fallback to multi-threaded CPU routines.

---

## Core Engines

### 1. RNA Hunter (RNA Target Screening)

Designed for RNA-targeted therapeutics, riboswitches, and viral RNA secondary motifs:

- **Virtual Target Screening:** Evaluates chemical libraries against target RNA sequences using a scoring model based on planar base-stacking, phosphate backbone electrostatics, and groove hydrogen-bonding networks.
- **Binding Site Discovery:** Resolves residue-level interaction propensities across the entire RNA sequence and highlights potential binding pockets.
- **Automated RNA Docking:** Automated rigid receptor conversion and AutoDock Vina docking for top candidate molecules.

### 2. Prozena (Protein Target Screening)

Engineered for high-throughput screening against macromolecular protein targets:

- **ConPLex Contrastive Co-Embedding:** Screen thousands of candidate compounds against protein sequences in seconds using the pre-trained ConPLex model (`weights/models/ConPLex_v1_BindingDB.pt`).
- **ADMET Pre-Screening Gate:** Enforces user-configurable thresholds for Molecular Weight, LogP, H-Bond Donors/Acceptors, TPSA, reactive functional groups, and PAINS alerts.
- **Automated Top-5 Docking and Ki Calculation:** Docks the top 5 ranked hits into the uploaded target PDB file, reporting binding affinity (kcal/mol) and estimated inhibition constants (Ki).
- **Interactive Inspection:** Integrated 3D structure viewer (py3Dmol) and 2D interaction fingerprint tags (ProLIF).

### 3. Refazenation (Target Prediction and Reverse Docking)

Engineered for polypharmacology, off-target identification, and drug repurposing:

- **Fast Offline Target Prediction:** Fast bioactivity profiling across major therapeutic families (Kinases, GPCRs, Nuclear Receptors, Ion Channels, Proteases, Epigenetic targets) based on curated bioactive chemotype fingerprints.
- **Multi-Target Reverse Docking:** Concurrently docks a query molecule against a panel of uploaded protein structures to detect binding preferences and repurposing avenues.
- **Binding Affinity Matrices:** Computes association constants (Ka), logKa, and estimated Ki across all target structures.

---

## Quick Start

If you already have Conda installed:

```bash
# 1. Clone the repository
git clone https://github.com/MrEgorin/fazena_tools.git
cd fazena_tools

# 2. Create and activate environment
conda create -n fazena_tools python=3.10 -y
conda activate fazena_tools

# 3. Install dependencies
pip install -r requirements.txt

# 4. Launch the web platform
streamlit run app.py
```

The web interface will be accessible at `http://localhost:8501`.

---

## Installation Guide

### Prerequisites

- Linux or Windows with WSL2 (Ubuntu 22.04 LTS or newer recommended)
- Python 3.10 or higher
- Optional: NVIDIA GPU with CUDA drivers (CUDA 12.0+)

### Step-by-Step Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/MrEgorin/fazena_tools.git
   cd fazena_tools
   ```

2. **Set up a dedicated environment:**
   ```bash
   conda create -n fazena_tools python=3.10 -y
   conda activate fazena_tools
   ```

3. **Install platform requirements:**
   ```bash
   pip install -r requirements.txt
   ```

4. **(Optional) PyTorch with GPU Acceleration:**
   For modern NVIDIA GPUs (including RTX 30, RTX 40, and RTX 50-series):
   ```bash
   pip install --upgrade torch --extra-index-url https://download.pytorch.org/whl/cu128
   ```

---

## Running the Application

### Starting the Web Interface

Run the application using Streamlit:

```bash
streamlit run app.py
```

Alternatively, you can launch via the entrypoint wrapper:

```bash
streamlit run app/app.py
```

### Navigating the Interface

Use the sidebar navigation to select your product workflow:
- **RNA Hunter:** Select *Virtual Target Screening* or *Binding Site Search*.
- **Prozena:** Select *HTS Protein Screening & Docking*.
- **Refazenation:** Select *Fast Target Prediction* or *Multi-Target Reverse Docking*.

Each screen provides interactive parameters, file uploaders, data tables, and export options (.CSV, .PDB).

---

## Repository Structure

```
fazena_tools/
|-- app.py                      # Primary Streamlit web application entrypoint
|-- app/
|   `-- app.py                  # Module redirect wrapper
|-- modules/
|   |-- config.py               # Hardware discovery, GPU detection, and paths
|   |-- input_parser.py         # Multi-format parser (SMILES, CSV, SDF, PDB, FASTA)
|   |-- admet.py                # Consolidated Lipinski, Veber, PAINS, Brenk, toxicophore engine
|   |-- docking_engine.py       # AutoDock Vina execution, grid box calc, ProLIF analysis
|   |-- conplex_screener.py     # ConPLex deep learning screening engine (GPU/CPU)
|   |-- rna_screener.py         # Structure-agnostic RNA virtual screening and binding sites
|   |-- target_prediction.py    # Fast offline bioactive target prediction engine
|   |-- ui_components.py        # FAZENA design system styling, dark theme, and cards
|   |-- prozena/                # Prozena HTS Protein Screening & Docking package
|   |   |-- __init__.py         # Package exports (ProteinScreening, render_page)
|   |   |-- protein_screening.py# ConPLex programmatic API wrapper
|   |   `-- ui.py               # Interactive Streamlit UI page
|   |-- refazenation/           # Refazenation Target Prediction & Reverse Docking package
|   |   |-- __init__.py         # Package exports (DrugRepurposing, render_page)
|   |   |-- drug_repurposing.py # Target prediction programmatic API wrapper
|   |   `-- ui.py               # Interactive Streamlit UI page
|   |-- virtual_screening.py    # RNA Hunter virtual screening UI page
|   |-- binding_prediction.py   # RNA Hunter binding site discovery UI page
|   `-- visualization/
|       |-- __init__.py         # Package exports for visualization utilities
|       `-- interaction_plots.py# 2D molecule rendering, RNA heatmaps, 3D py3Dmol viewer
|-- weights/
|   `-- models/
|       `-- ConPLex_v1_BindingDB.pt # Pre-trained ConPLex neural network weights
|-- tests/
|   `-- test_platform.py        # Automated test suite (unit and integration tests)
|-- requirements.txt            # Python dependencies
`-- README.md                   # Platform documentation
```

---

## Input Formats and Compatibility

FAZENA Tools supports standard computational chemistry formats:

| Data Type | Supported Formats | Example Content |
|:---|:---|:---|
| **Small Molecules** | SMILES text, .smi, .sdf, .mol, .csv, .tsv | One SMILES per line or multi-column table with compound names |
| **Protein Sequences** | FASTA, raw single-letter text | Standard 20 amino acid sequence (e.g., `MKTVRQERL...`) |
| **RNA Sequences** | FASTA, raw single-letter text | RNA nucleotide sequence (`AUCG...`) |
| **3D Structures** | .pdb, .pdbqt | Standard Brookhaven PDB format with ATOM/HETATM records |

---

## Automated Test Suite

The platform includes an automated unit and integration test suite to verify module integrity, model loading, and computational pipelines.

To execute the test suite:

```bash
python3 -m unittest discover tests
```

Expected output:
```
Ran 6 tests in 4.5s

OK
```

---

## Authors and Citation

Developed by **Yehor Franchuk, MrEgorin, and the FAZENA Team**.  
Official Website: [https://fazena.org](https://fazena.org)

If you use FAZENA Tools in your academic or industrial research, please cite:

```bibtex
@software{fazena_tools_2026,
  author = {Yehor Franchuk and MrEgorin},
  title = {FAZENA Tools: Modular Bio-Computational Drug Discovery Platform},
  year = {2026},
  url = {https://github.com/MrEgorin/fazena_tools}
}
```
