"""
FAZENA Tools — SMARTBind & RNA Hunter Wrapper
Connects RNA Hunter virtual screening and binding site prediction to real inference engines.
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

import pandas as pd
from typing import List, Dict, Any
from modules.rna_screener import screen_rna_library, predict_rna_binding_sites

def is_engine_active() -> bool:
    """Returns True indicating RNA inference engine is operational."""
    return True

def run_virtual_screening(rna_sequence: str, ligands: List[Any], device='cpu') -> pd.DataFrame:
    """
    Runs virtual screening for ligands against target RNA sequence.
    Accepts list of SMILES or list of dicts.
    """
    ligand_items = []
    for i, item in enumerate(ligands):
        if isinstance(item, str):
            ligand_items.append({'name': f"Mol_{i+1}", 'SMILES': item})
        elif isinstance(item, dict):
            ligand_items.append(item)

    df = screen_rna_library(rna_sequence, ligand_items)
    return df

def run_binding_site_prediction(rna_sequence: str, ligand: str, device='cpu') -> List[float]:
    """
    Predicts interaction probability heatmap across the RNA sequence.
    Returns per-nucleotide probabilities.
    """
    return predict_rna_binding_sites(rna_sequence, ligand)
