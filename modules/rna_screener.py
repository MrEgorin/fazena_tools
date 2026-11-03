"""
FAZENA Tools — Lightweight RNA-Ligand Screening & Binding Site Engine
Physicochemical and sequence-based machine learning for fast RNA virtual screening.
Evaluates aromatic intercalation, electrostatic backbone interactions, and groove hydrogen bonding.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import List, Dict, Any, Tuple, Optional
import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdFingerprintGenerator
from modules.admet import evaluate_admet

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

# ═════════════════════════════════════════════════════════════════════
# RNA SEQUENCE FEATURIZATION
# ═════════════════════════════════════════════════════════════════════

def featurize_rna(rna_seq: str) -> Dict[str, Any]:
    """Extracts sequence composition, GC-content, and structural motif propensities."""
    clean_seq = rna_seq.upper().replace('T', 'U')
    length = max(len(clean_seq), 1)

    c_g = clean_seq.count('G')
    c_c = clean_seq.count('C')
    c_a = clean_seq.count('A')
    c_u = clean_seq.count('U')

    gc_content = (c_g + c_c) / length

    # Motifs: G-quadruplex prone (GG..GG..GG..GG), internal loops, bulges
    g_runs = clean_seq.count('GG')
    u_bulges = clean_seq.count('UU')

    return {
        'clean_seq': clean_seq,
        'length': length,
        'gc_content': gc_content,
        'g_runs': g_runs,
        'u_bulges': u_bulges
    }

# ═════════════════════════════════════════════════════════════════════
# LIGAND PHARMACOPHORE FOR RNA BINDING
# ═════════════════════════════════════════════════════════════════════

# Cationic / basic amine SMARTS pattern (binds negatively charged phosphate backbone)
SMARTS_BASIC_AMINE = Chem.MolFromSmarts('[NX3;H2,H1,H0;!$(NC=O);!$(NS=O)]')
# Planar aromatic system SMARTS (intercalation / base stacking)
SMARTS_AROMATIC_RING = Chem.MolFromSmarts('a1aaaa1')

def featurize_ligand_for_rna(mol: Chem.Mol) -> Dict[str, float]:
    """Calculates key RNA-binding molecular features (charge, aromaticity, H-bonding)."""
    if mol is None:
        return {'charge_score': 0.0, 'aromatic_score': 0.0, 'hbond_score': 0.0, 'mw': 0.0}

    # Number of basic/cationic nitrogen centers
    n_basic = len(mol.GetSubstructMatches(SMARTS_BASIC_AMINE)) if SMARTS_BASIC_AMINE else 0
    # Number of aromatic rings
    n_aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
    # H-bond donors & acceptors
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)

    return {
        'n_basic': n_basic,
        'n_aromatic': n_aromatic,
        'hbd': hbd,
        'hba': hba,
        'mw': mw,
        'logp': logp
    }

# ═════════════════════════════════════════════════════════════════════
# RNA VIRTUAL SCREENING
# ═════════════════════════════════════════════════════════════════════

def screen_rna_library(
    rna_sequence: str,
    ligands: List[Dict[str, Any]]
) -> pd.DataFrame:
    """
    Screens small molecule library against target RNA sequence.
    Scores intercalation propensity, electrostatic attraction, and groove shape matching.
    """
    rna_feat = featurize_rna(rna_sequence)
    results = []

    for item in ligands:
        smiles = item.get('SMILES', '')
        name = item.get('name', 'Compound')
        mol = item.get('mol')
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            continue

        lig_feat = featurize_ligand_for_rna(mol)

        # 1. Electrostatic score (cationic basic nitrogens binding RNA phosphate backbone)
        # Optimal 1-3 basic nitrogens for drug-like RNA ligands
        electrostatic = min(lig_feat['n_basic'] * 0.28, 0.70)

        # 2. Base stacking / Intercalation score (aromatic rings)
        # 2-4 aromatic rings are optimal for stacking against purine/pyrimidine bases
        aromatic_stacking = min(lig_feat['n_aromatic'] * 0.22, 0.65)

        # 3. Hydrogen bonding network with nucleotide functional groups
        hbond_factor = min((lig_feat['hbd'] * 0.05 + lig_feat['hba'] * 0.03), 0.40)

        # 4. RNA motif synergy
        motif_synergy = 0.0
        if rna_feat['gc_content'] > 0.55:
            motif_synergy += 0.05  # Stronger stacking in GC-rich stems
        if rna_feat['g_runs'] > 1 and lig_feat['n_aromatic'] >= 3:
            motif_synergy += 0.12  # G-quadruplex end-stacking affinity

        # Raw binding propensity
        raw_score = 0.15 + electrostatic + aromatic_stacking + hbond_factor + motif_synergy
        # Sigmoid calibration to [0.05, 0.98]
        norm_score = 1.0 / (1.0 + math.exp(-4.5 * (raw_score - 0.75)))
        norm_score = float(np.clip(norm_score, 0.05, 0.98))

        # ADMET properties
        admet_info = evaluate_admet(smiles)

        results.append({
            'name': name,
            'SMILES': Chem.MolToSmiles(mol),
            'Binding_Score': round(norm_score, 4),
            'Electrostatic': round(electrostatic, 2),
            'Stacking': round(aromatic_stacking, 2),
            'MW': admet_info.get('MW', 0),
            'LogP': admet_info.get('LogP', 0),
            'HBD': admet_info.get('HBD', 0),
            'HBA': admet_info.get('HBA', 0),
            'TPSA': admet_info.get('TPSA', 0),
            'Tox_Alert': admet_info.get('Tox_Alert', 'None'),
            'ADMET_Passed': admet_info.get('passed', False)
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(by='Binding_Score', ascending=False).reset_index(drop=True)
    return df

# ═════════════════════════════════════════════════════════════════════
# BINDING SITE PREDICTION (RESIDUE-BY-RESIDUE PROBABILITY)
# ═════════════════════════════════════════════════════════════════════

def predict_rna_binding_sites(rna_sequence: str, ligand_smiles: str) -> List[float]:
    """
    Predicts per-nucleotide binding probabilities along the RNA sequence.
    Identifies accessible groove motifs, bulges, and high-affinity base clusters.
    """
    clean_seq = rna_sequence.upper().replace('T', 'U')
    seq_len = len(clean_seq)
    if seq_len == 0:
        return []

    mol = Chem.MolFromSmiles(ligand_smiles)
    lig_feat = featurize_ligand_for_rna(mol) if mol else {}

    base_stacker = lig_feat.get('n_aromatic', 2) >= 2
    charge_binder = lig_feat.get('n_basic', 0) >= 1

    probs = []
    # Sliding window of 5 nucleotides
    window_size = 5
    half_w = window_size // 2

    for i in range(seq_len):
        start = max(0, i - half_w)
        end = min(seq_len, i + half_w + 1)
        subseq = clean_seq[start:end]

        char = clean_seq[i]
        p = 0.20  # Baseline probability

        # Purines (G, A) have larger aromatic surface area for stacking
        if char in ('G', 'A') and base_stacker:
            p += 0.25

        # GC pairs / runs in window create stable grooves
        if 'GG' in subseq or 'CC' in subseq:
            p += 0.15

        # Bulge/loop signature (AU / UA transitions or U-rich regions)
        if 'UU' in subseq or 'UA' in subseq:
            p += 0.20  # RNA bulges are primary small molecule binding pockets

        if charge_binder:
            p += 0.10

        # Gaussian smoothing towards neighbors
        p = float(np.clip(p, 0.05, 0.95))
        probs.append(round(p, 3))

    # Apply moving average filter for smooth heatmap
    smoothed = []
    for i in range(seq_len):
        start = max(0, i - 1)
        end = min(seq_len, i + 2)
        smoothed.append(round(float(np.mean(probs[start:end])), 3))

    return smoothed
