"""
FAZENA Tools — AutoDocking Class Wrapper
Wraps the real AutoDock Vina engine in modules/docking_engine.py
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import Dict, List, Optional
import os
from modules.docking_engine import dock_single_compound, calculate_grid_box

class AutoDocking:
    """Provides high-level programmatic docking API for protein and RNA targets."""
    def __init__(self, vina_path: str = "vina"):
        self.vina_path = vina_path

    def dock_protein(self, protein_file: str, ligand_file: str, **kwargs) -> Dict[str, float]:
        """Docks ligand to protein using AutoDock Vina."""
        with open(protein_file, 'r') as f:
            prot_str = f.read()
        with open(ligand_file, 'r') as f:
            lig_str = f.read()
        
        # If ligand is SMILES string
        smi = lig_str.strip().split()[0]
        res = dock_single_compound(prot_str, smi, is_rna=False)
        return {'affinity': res.get('affinity', 0.0), 'ki': res.get('ki_estimated', 'N/A')}

    def dock_rna(self, rna_file: str, ligand_file: str, **kwargs) -> Dict[str, float]:
        """Docks ligand to RNA using AutoDock Vina."""
        with open(rna_file, 'r') as f:
            rna_str = f.read()
        with open(ligand_file, 'r') as f:
            lig_str = f.read()

        smi = lig_str.strip().split()[0]
        res = dock_single_compound(rna_str, smi, is_rna=True)
        return {'affinity': res.get('affinity', 0.0), 'ki': res.get('ki_estimated', 'N/A')}

    def auto_configure_grid(self, protein_file: str) -> Dict[str, float]:
        with open(protein_file, 'r') as f:
            prot_str = f.read()
        return calculate_grid_box(prot_str)
