"""
FAZENA Tools — RNAHunter Compatibility Class
Wraps RNA screening & binding site engine in modules/rna_screener.py
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

from typing import Dict, List, Optional
from modules.rna_screener import screen_rna_library, predict_rna_binding_sites

class RNAHunter:
    """Programmatic API for RNA virtual screening and binding site prediction."""
    def __init__(self, model_path: Optional[str] = None):
        self.model_path = model_path

    def predict(self, rna_sequence: str, smiles_list: Optional[List[str]] = None) -> Dict[str, float]:
        """Predicts binding scores for smiles list against target RNA sequence."""
        if not smiles_list:
            smiles_list = ["CCO"]
        ligands = [{'name': f"Mol_{i+1}", 'SMILES': s} for i, s in enumerate(smiles_list)]
        df = screen_rna_library(rna_sequence, ligands)
        return dict(zip(df['SMILES'], df['Binding_Score']))

    def predict_binding_sites(self, rna_sequence: str, ligand_smiles: str) -> List[float]:
        return predict_rna_binding_sites(rna_sequence, ligand_smiles)
