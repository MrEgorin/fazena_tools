"""
FAZENA Tools — AutoDocking Class Wrapper
Wraps the real AutoDock Vina engine in modules/docking_engine.py
Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

from typing import Dict, Any, Optional
import os
from modules.docking_engine import dock_single_compound, calculate_grid_box

class AutoDocking:
    """Provides high-level programmatic docking API for protein and RNA targets."""
    def __init__(self, vina_path: str = "vina"):
        self.vina_path = vina_path

    def auto_configure_grid(self, receptor_file: str) -> Dict[str, float]:
        """Calculates docking search box coordinates and dimensions from receptor structure."""
        with open(receptor_file, 'r') as f:
            rec_str = f.read()
        return calculate_grid_box(rec_str)

    def dock_protein(self, protein_file: str, ligand_file: str, **kwargs) -> Dict[str, Any]:
        """
        Docks ligand to protein target using AutoDock Vina.
        Accepts SMILES file, SDF library, or molecule file.
        """
        with open(protein_file, 'r') as f:
            prot_str = f.read()

        ligand_input = ligand_file
        if os.path.isfile(ligand_file):
            # If it's a short 1-line SMILES file
            ext = os.path.splitext(ligand_file)[1].lower()
            if ext in ('.smi', '.smiles', '.txt'):
                with open(ligand_file, 'r') as f:
                    content = f.read().strip()
                if content and '\n' not in content and ' ' not in content:
                    ligand_input = content

        res = dock_single_compound(prot_str, ligand_input, is_rna=False, **kwargs)
        return {
            'affinity': res.get('affinity', 0.0),
            'ki': res.get('ki_estimated', 'N/A'),
            'status': res.get('status', 'Failed'),
            'complex_pdb': res.get('complex_pdb', ''),
            'top_pose_pdbqt': res.get('top_pose_pdbqt', ''),
            'error': res.get('error', None)
        }

    def dock_rna(self, rna_file: str, ligand_file: str, **kwargs) -> Dict[str, Any]:
        """
        Docks ligand to RNA target using AutoDock Vina with nucleic-acid atom typing.
        Accepts SMILES, SDF library, or 3D molecule file.
        """
        with open(rna_file, 'r') as f:
            rna_str = f.read()

        ligand_input = ligand_file
        if os.path.isfile(ligand_file):
            ext = os.path.splitext(ligand_file)[1].lower()
            if ext in ('.smi', '.smiles', '.txt'):
                with open(ligand_file, 'r') as f:
                    content = f.read().strip()
                if content and '\n' not in content and ' ' not in content:
                    ligand_input = content

        res = dock_single_compound(rna_str, ligand_input, is_rna=True, **kwargs)
        return {
            'affinity': res.get('affinity', 0.0),
            'ki': res.get('ki_estimated', 'N/A'),
            'status': res.get('status', 'Failed'),
            'complex_pdb': res.get('complex_pdb', ''),
            'top_pose_pdbqt': res.get('top_pose_pdbqt', ''),
            'error': res.get('error', None)
        }

    def dock(self, target_file: str, ligand_file: str, target_type: str = 'protein', **kwargs) -> Dict[str, Any]:
        """Automated docking based on target type ('protein' or 'rna')."""
        t_type = str(target_type).lower()
        if t_type == 'rna':
            return self.dock_rna(target_file, ligand_file, **kwargs)
        elif t_type == 'protein':
            return self.dock_protein(target_file, ligand_file, **kwargs)
        else:
            raise ValueError(f"Invalid target type '{target_type}'. Supported types: 'protein', 'rna'")
