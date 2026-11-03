"""
FAZENA Tools — Prozena Protein Screening Wrapper
Connects Prozena to ConPLex neural network engine and BIND.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import List, Dict, Any, Tuple
import pandas as pd
from modules.conplex_screener import run_conplex_screening

def run_bind_real(protein_seq: str, ligands: List[Any], device: str = "cpu") -> Tuple[pd.DataFrame, bool]:
    """
    Executes high-throughput protein-ligand screening using ConPLex neural network.
    Accepts SMILES list or ligand dicts.
    Returns: (df, is_real)
    """
    ligand_items = []
    for i, item in enumerate(ligands):
        if isinstance(item, str):
            ligand_items.append({'name': f"Compound_{i+1}", 'SMILES': item})
        elif isinstance(item, dict):
            ligand_items.append(item)

    df, is_real = run_conplex_screening(protein_seq, ligand_items)
    return df, is_real
