"""
FAZENA Tools — Fast Target Prediction Engine (SwissTargetPrediction Alternative)
Lightweight target prediction using 2D Morgan fingerprints and curated bioactive target chemotypes.
Zero heavy download, fast offline inference on Linux/cross-platform.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import List, Dict, Any, Optional
import math
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _get_fp(mol):
    return _MORGAN_GEN.GetFingerprint(mol)

# ═════════════════════════════════════════════════════════════════════
# REFERENCE BIOACTIVE TARGET KNOWLEDGE BASE
# Curated target profiles with representative clinical/tool compounds
# ═════════════════════════════════════════════════════════════════════

TARGET_PROFILES = [
    # Kinases
    {
        'target_name': 'Epidermal Growth Factor Receptor (EGFR)',
        'class': 'Kinase',
        'uniprot': 'P00533',
        'pdb_id': '1M17',
        'compounds': [
            ('Gefitinib', 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1'),
            ('Erlotinib', 'COCCOC1=C(OCCOC)C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C'),
            ('Osimertinib', 'CN(C)CC=CC(=O)Nc1cc(Nc2ncccn2)c(OC)cc1N(C)CCN(C)C')
        ]
    },
    {
        'target_name': 'Cyclin-Dependent Kinase 2 (CDK2)',
        'class': 'Kinase',
        'uniprot': 'P24941',
        'pdb_id': '1HCK',
        'compounds': [
            ('Roscovitine', 'CCN(C1=NC=NC2=C1N=C(N2CC(C)C)NCC3=CC=CC=C3)C(C)CO'),
            ('Dinaciclib', 'CCN1CCN(CC1)C(=O)C2=CC=C(C=C2)NC3=NC=C4C(=N3)N(N=C4)C5CCCC5'),
            ('Flavopiridol', 'CN1CCC(C(C1)O)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=CC=C4Cl)O)O')
        ]
    },
    {
        'target_name': 'BRAF Kinase (V600E)',
        'class': 'Kinase',
        'uniprot': 'P15056',
        'pdb_id': '4MNE',
        'compounds': [
            ('Vemurafenib', 'CCCS(=O)(=O)NC1=C(C=C(C(=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F'),
            ('Dabrafenib', 'CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C=C(C(=C3)F)NS(=O)(=O)C4=C(C=CC=C4F)F)F')
        ]
    },
    {
        'target_name': 'Tyrosine-Protein Kinase ABL1',
        'class': 'Kinase',
        'uniprot': 'P00519',
        'pdb_id': '1IEP',
        'compounds': [
            ('Imatinib', 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'),
            ('Nilotinib', 'Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1'),
            ('Dasatinib', 'Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1')
        ]
    },
    # GPCRs
    {
        'target_name': 'Serotonin 2A Receptor (5-HT2A)',
        'class': 'GPCR (Class A)',
        'uniprot': 'P28223',
        'pdb_id': '6A93',
        'compounds': [
            ('Ketanserin', 'O=C1NC(=O)N(CCN2CCC(CC2)C(=O)C3=CC=C(F)C=C3)C1=O'),
            ('Risperidone', 'CC1=C(CCN2CCC(CC2)C3=NOC4=C3C=CC(=C4)F)C(=O)N5CCCC5=N1'),
            ('Clozapine', 'CN1CCN(CC1)C2=NC3=CC=CC=C3NC4=C2C=C(C=C4)Cl')
        ]
    },
    {
        'target_name': 'Dopamine D2 Receptor (DRD2)',
        'class': 'GPCR (Class A)',
        'uniprot': 'P14416',
        'pdb_id': '6CM4',
        'compounds': [
            ('Haloperidol', 'OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1'),
            ('Aripiprazole', 'O=C1CCc2ccc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)cc2N1'),
            ('Sulpiride', 'CCN1CCCC1CNC(=O)c2cc(S(=O)(=O)N)ccc2OC')
        ]
    },
    {
        'target_name': 'Cannabinoid Receptor 1 (CB1)',
        'class': 'GPCR (Class A)',
        'uniprot': 'P21554',
        'pdb_id': '5XRA',
        'compounds': [
            ('Rimonabant', 'Clc1ccc(cc1)-c1nn(c(Cl)c1C(=O)NN1CCCCC1)-c1ccc(Cl)cc1Cl'),
            ('THC', 'CCCCCC1=CC(O)=C2C(=C1)OC(C)(C)C1CC=C(C)CC21')
        ]
    },
    {
        'target_name': 'Beta-2 Adrenergic Receptor (ADRB2)',
        'class': 'GPCR (Class A)',
        'uniprot': 'P07550',
        'pdb_id': '2RH1',
        'compounds': [
            ('Salbutamol', 'CC(C)(C)NCC(O)c1ccc(O)c(CO)c1'),
            ('Propranolol', 'CC(C)NCC(O)COc1cccc2ccccc12'),
            ('Isoproterenol', 'CC(C)NCC(O)c1ccc(O)c(O)c1')
        ]
    },
    # Proteases & Enzymes
    {
        'target_name': 'Cyclooxygenase-2 (COX-2)',
        'class': 'Enzyme',
        'uniprot': 'P35354',
        'pdb_id': '1CX2',
        'compounds': [
            ('Celecoxib', 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1'),
            ('Rofecoxib', 'O=S(=O)(N)c1ccc(cc1)C1=C(c2ccccc2)COC1=O'),
            ('Valdecoxib', 'Cc1onc(c1-c1ccccc1)c1ccc(cc1)S(N)(=O)=O')
        ]
    },
    {
        'target_name': 'Angiotensin-Converting Enzyme (ACE)',
        'class': 'Protease',
        'uniprot': 'P12821',
        'pdb_id': '1O86',
        'compounds': [
            ('Captopril', 'CC(CS)C(=O)N1CCCC1C(=O)O'),
            ('Enalapril', 'CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(=O)O'),
            ('Lisinopril', 'NCCCCC(NC(CCc1ccccc1)C(=O)O)C(=O)N1CCCC1C(=O)O')
        ]
    },
    {
        'target_name': 'SARS-CoV-2 Main Protease (Mpro)',
        'class': 'Viral Protease',
        'uniprot': 'P0DTD1',
        'pdb_id': '6LU7',
        'compounds': [
            ('Nirmatrelvir', 'CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C'),
            ('Ensitrelvir', 'CN1C(=NN=C1C(=O)NC2=CC(=C(C=C2)F)F)C3=C(C=CC(=C3)F)F')
        ]
    },
    # Nuclear Receptors
    {
        'target_name': 'Estrogen Receptor Alpha (ESR1)',
        'class': 'Nuclear Receptor',
        'uniprot': 'P03372',
        'pdb_id': '1ERR',
        'compounds': [
            ('Tamoxifen', 'CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1'),
            ('Raloxifene', 'Oc1ccc(cc1)C(=O)c1c(sc2cc(O)ccc12)c1ccc(OCCN2CCCCC2)cc1'),
            ('Estradiol', 'CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O')
        ]
    },
    {
        'target_name': 'Peroxisome Proliferator-Activated Receptor Gamma (PPARG)',
        'class': 'Nuclear Receptor',
        'uniprot': 'P37231',
        'pdb_id': '2PRG',
        'compounds': [
            ('Rosiglitazone', 'CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1'),
            ('Pioglitazone', 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')
        ]
    },
    # Ion Channels & Toxicity Off-Targets
    {
        'target_name': 'hERG Potassium Channel (KCNH2)',
        'class': 'Ion Channel (Cardiac safety)',
        'uniprot': 'Q12809',
        'pdb_id': '5VA1',
        'compounds': [
            ('Astemizole', 'COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1'),
            ('Dofetilide', 'CNS(=O)(=O)c1ccc(OCCN(CCc2ccc(NS(C)(=O)=O)cc2)C)cc1'),
            ('Terfenadine', 'CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1')
        ]
    }
]

# Pre-compute reference fingerprints
_PRECOMPUTED = False

def _ensure_precomputed():
    global _PRECOMPUTED
    if _PRECOMPUTED:
        return
    for t in TARGET_PROFILES:
        fps = []
        for name, smi in t['compounds']:
            m = Chem.MolFromSmiles(smi)
            if m:
                fp = _get_fp(m)
                fps.append((name, smi, fp))
        t['cached_fps'] = fps
    _PRECOMPUTED = True

# ═════════════════════════════════════════════════════════════════════
# TARGET PREDICTION INFERENCE
# ═════════════════════════════════════════════════════════════════════

def predict_targets(smiles_or_mol: Any, top_k: int = 10) -> pd.DataFrame:
    """
    Predicts biological macromolecular targets for a small molecule query.
    Calculates Tanimoto similarity to known bioactive chemotypes and models target probability.
    """
    _ensure_precomputed()

    if isinstance(smiles_or_mol, str):
        query_mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        query_mol = smiles_or_mol

    if query_mol is None:
        return pd.DataFrame()

    query_fp = _get_fp(query_mol)
    records = []

    for target in TARGET_PROFILES:
        cached = target.get('cached_fps', [])
        if not cached:
            continue

        best_sim = 0.0
        best_comp_name = ""
        best_comp_smi = ""

        sims = []
        for name, smi, ref_fp in cached:
            sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            sims.append(sim)
            if sim > best_sim:
                best_sim = sim
                best_comp_name = name
                best_comp_smi = smi

        avg_sim = float(np.mean(sims)) if sims else 0.0

        # Calibrated logistic probability model based on chemical similarity:
        # P = 1 / (1 + exp(-k * (sim - x0)))
        # Tuned so high similarity (>0.5) has high confidence (>80%)
        # and background random similarity (~0.15-0.20) has low probability (<15%)
        k_factor = 14.0
        x0_threshold = 0.38
        effective_sim = 0.7 * best_sim + 0.3 * avg_sim
        prob = 1.0 / (1.0 + math.exp(-k_factor * (effective_sim - x0_threshold)))
        prob = min(max(prob, 0.05), 0.99)

        records.append({
            'Target Name': target['target_name'],
            'Class': target['class'],
            'UniProt ID': target['uniprot'],
            'PDB ID': target['pdb_id'],
            'Probability': round(prob * 100, 1),
            'Tanimoto Similarity': round(best_sim, 3),
            'Reference Compound': best_comp_name,
            'Reference SMILES': best_comp_smi
        })

    df = pd.DataFrame(records)
    df = df.sort_values(by=['Probability', 'Tanimoto Similarity'], ascending=[False, False]).reset_index(drop=True)
    df.index = df.index + 1
    return df.head(top_k)
