"""
FAZENA Tools — ProteinScreening Compatibility Class
Wraps ConPLex protein screening engine in modules/conplex_screener.py
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

from typing import List, Dict, Optional
from modules.conplex_screener import run_conplex_screening

class ProteinScreening:
    """Programmatic API for protein screening using ConPLex neural network."""
    def __init__(self, model_path: Optional[str] = None):
        self.model_path = model_path

    def predict(self, protein_seq: str, smiles_list: Optional[List[str]] = None) -> Dict[str, float]:
        """Predicts pKi binding affinities for smiles list against target protein sequence."""
        if not smiles_list:
            smiles_list = ["CCO"]
        ligand_dicts = [{'name': f"Mol_{i+1}", 'SMILES': s} for i, s in enumerate(smiles_list)]
        df, _ = run_conplex_screening(protein_seq, ligand_dicts)
        return dict(zip(df['SMILES'], df['pKi']))
