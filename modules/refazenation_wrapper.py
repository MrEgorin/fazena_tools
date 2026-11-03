"""
FAZENA Tools — Refazenation Wrapper
Wraps real AutoDock Vina docking and Target Prediction engines.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import List, Dict, Any, Optional
import os
import math
import pandas as pd
from rdkit import Chem
from modules.docking_engine import (
    dock_single_compound,
    prepare_ligand_pdbqt,
    prepare_receptor_pdbqt,
    calculate_grid_box,
    analyze_interactions
)
from modules.target_prediction import predict_targets

def run_docking_real(
    ligand_mol_or_smiles: Any,
    pdb_files: List[Any],
    grid_box: Optional[Dict[str, float]] = None
) -> List[Dict[str, Any]]:
    """
    Coordinates reverse docking of a ligand across multiple uploaded target protein structures.
    """
    results = []
    
    if isinstance(ligand_mol_or_smiles, str):
        smiles = ligand_mol_or_smiles
    else:
        smiles = Chem.MolToSmiles(ligand_mol_or_smiles)

    for pdb_file in pdb_files:
        prot_name = getattr(pdb_file, 'name', 'Target_Protein')
        
        # Read content
        if hasattr(pdb_file, 'getvalue'):
            pdb_str = pdb_file.getvalue().decode('utf-8', errors='ignore')
        elif hasattr(pdb_file, 'read'):
            pdb_str = pdb_file.read()
            if isinstance(pdb_str, bytes):
                pdb_str = pdb_str.decode('utf-8', errors='ignore')
        else:
            pdb_str = str(pdb_file)

        dock_res = dock_single_compound(
            receptor_pdb=pdb_str,
            ligand_smiles=smiles,
            compound_name=prot_name,
            grid_box=grid_box,
            exhaustiveness=6,
            n_poses=3
        )

        if dock_res.get('status') == 'Success':
            score = dock_res.get('affinity', 0.0)
            ka = math.exp(-score / (1.987e-3 * 298.15))
            logka = round(math.log10(max(ka, 1.0)), 2)
            results.append({
                "Protein": prot_name,
                "Score": score,
                "Ka": f"{ka:.2e}",
                "logKa": logka,
                "Ki": dock_res.get('ki_estimated', 'N/A'),
                "RMSD": 0.0,
                "Status": "Docked Successfully",
                "complex_pdb": dock_res.get('complex_pdb', '')
            })
        else:
            results.append({
                "Protein": prot_name,
                "Score": 0.0,
                "Ka": "0.0",
                "logKa": 0.0,
                "Ki": "N/A",
                "RMSD": 0.0,
                "Status": f"Failed: {dock_res.get('error', 'Unknown')}",
                "complex_pdb": ""
            })

    # Sort by best binding score (most negative)
    results.sort(key=lambda x: x.get('Score', 0.0))
    return results

def create_complex_pdb(protein_pdb_str: str, ligand_mol: Chem.Mol) -> str:
    """Merges protein PDB and ligand PDB into one standard PDB string."""
    ligand_pdb = Chem.MolToPDBBlock(ligand_mol)
    clean_prot = "\n".join([l for l in protein_pdb_str.splitlines() if not l.startswith(('END', 'CONECT'))])
    clean_lig = "\n".join([l for l in ligand_pdb.splitlines() if l.startswith(('ATOM', 'HETATM', 'CONECT'))])
    return clean_prot + "\n" + clean_lig + "\nEND\n"
