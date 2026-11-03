"""
FAZENA Tools — Flexible Multi-Format Input & Output Parser
Handles multi-line SMILES, SDF, CSV, TSV, FASTA, PDB, and PDBQT with auto-detection.
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

from typing import List, Tuple, Dict, Any, Optional
import io
import re
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# ═════════════════════════════════════════════════════════════════════
# LIGAND / CHEMICAL LIBRARY PARSER
# ═════════════════════════════════════════════════════════════════════

def sanitize_mol_safe(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Applies relaxed sanitization to avoid failing on non-standard hypervalent atoms."""
    if mol is None:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return mol
    except Exception:
        return mol

def parse_ligand_text(text: str) -> List[Dict[str, Any]]:
    """
    Parses pasted text containing one or many SMILES strings.
    Supports formats:
    - One SMILES per line
    - SMILES with molecule name: "CCO Ethanol" or "CC(=O)O\tAspirin"
    - Comma-separated: "CCO, c1ccccc1, CC(=O)O"
    """
    if not text or not text.strip():
        return []

    results = []
    # If text contains commas and no newlines, split by comma
    if ',' in text and '\n' not in text:
        tokens = [t.strip() for t in text.split(',') if t.strip()]
        for i, token in enumerate(tokens):
            mol = Chem.MolFromSmiles(token)
            if mol:
                results.append({
                    'id': f"Mol_{i+1}",
                    'name': f"Compound_{i+1}",
                    'SMILES': Chem.MolToSmiles(mol),
                    'mol': mol
                })
        return results

    # Standard line-by-line parsing
    lines = text.strip().splitlines()
    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        smiles_candidate = parts[0]
        name = parts[1] if len(parts) > 1 else f"Compound_{i+1}"
        
        mol = Chem.MolFromSmiles(smiles_candidate)
        if mol:
            results.append({
                'id': f"Mol_{i+1}",
                'name': name,
                'SMILES': Chem.MolToSmiles(mol),
                'mol': mol
            })

    return results

def parse_ligand_file(file_bytes: bytes, filename: str) -> List[Dict[str, Any]]:
    """
    Parses an uploaded chemical library file.
    Supports: .sdf, .mol, .smi, .smiles, .csv, .tsv, .txt
    """
    results = []
    fname = filename.lower()

    # 1. SDF format
    if fname.endswith(('.sdf', '.mol')):
        try:
            suppl = Chem.ForwardSDMolSupplier(io.BytesIO(file_bytes), sanitize=False)
            for i, mol in enumerate(suppl):
                if mol is None:
                    continue
                mol = sanitize_mol_safe(mol)
                name = mol.GetProp('_Name') if mol.HasProp('_Name') and mol.GetProp('_Name').strip() else f"Mol_{i+1}"
                try:
                    smi = Chem.MolToSmiles(mol)
                    results.append({
                        'id': f"Mol_{i+1}",
                        'name': name,
                        'SMILES': smi,
                        'mol': mol
                    })
                except Exception:
                    continue
            return results
        except Exception as e:
            pass

    # 2. CSV / TSV format
    if fname.endswith(('.csv', '.tsv', '.txt')) and (b',' in file_bytes[:1000] or b'\t' in file_bytes[:1000]):
        try:
            sep = '\t' if (fname.endswith('.tsv') or b'\t' in file_bytes[:1000]) else ','
            df = pd.read_csv(io.BytesIO(file_bytes), sep=sep)
            
            # Find SMILES column
            smiles_col = None
            potential_smi_cols = ['smiles', 'SMILES', 'Smiles', 'canonical_smiles', 'structure', 'Structure', 'mol', 'MOL']
            for col in df.columns:
                if str(col).strip() in potential_smi_cols:
                    smiles_col = col
                    break
            
            # If not found by name, check first 3 columns for valid SMILES
            if smiles_col is None:
                for col in df.columns[:3]:
                    sample = df[col].dropna().head(5).astype(str)
                    if any(Chem.MolFromSmiles(s) is not None for s in sample):
                        smiles_col = col
                        break

            # Find Name/ID column
            name_col = None
            potential_name_cols = ['name', 'Name', 'id', 'ID', 'compound_id', 'Compound_ID', 'Molecule', 'title']
            for col in df.columns:
                if str(col).strip() in potential_name_cols and col != smiles_col:
                    name_col = col
                    break

            if smiles_col is not None:
                for i, row in df.iterrows():
                    smi_str = str(row[smiles_col]).strip()
                    mol = Chem.MolFromSmiles(smi_str)
                    if mol:
                        name = str(row[name_col]).strip() if (name_col and pd.notna(row[name_col])) else f"Compound_{i+1}"
                        results.append({
                            'id': f"Mol_{i+1}",
                            'name': name,
                            'SMILES': Chem.MolToSmiles(mol),
                            'mol': mol
                        })
                return results
        except Exception:
            pass

    # 3. Plain SMI / TXT (one per line)
    try:
        text = file_bytes.decode('utf-8', errors='ignore')
        return parse_ligand_text(text)
    except Exception:
        return []

# ═════════════════════════════════════════════════════════════════════
# BIOLOGICAL SEQUENCE PARSER (RNA & PROTEIN)
# ═════════════════════════════════════════════════════════════════════

def parse_sequence_input(text: str, seq_type: str = 'rna') -> Tuple[str, str]:
    """
    Parses a sequence from FASTA or raw string.
    Returns: (clean_sequence, sequence_name)
    """
    if not text or not text.strip():
        return "", ""

    lines = text.strip().splitlines()
    name = "Target"
    seq_lines = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:].strip()
        else:
            seq_lines.append(line)

    raw_seq = "".join(seq_lines).replace(" ", "").upper()

    if seq_type.lower() == 'rna':
        # Auto-convert Thymine (T) to Uracil (U) for RNA workflows
        clean_seq = raw_seq.replace('T', 'U')
        # Filter to standard RNA bases
        clean_seq = re.sub(r'[^AUCGN]', '', clean_seq)
        return clean_seq, name
    else:
        # Standard 20 amino acids + ambiguous
        clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWYX]', '', raw_seq)
        return clean_seq, name

# ═════════════════════════════════════════════════════════════════════
# 3D STRUCTURE / PDB PARSER
# ═════════════════════════════════════════════

def clean_pdb_and_calc_box(pdb_content: str, padding: float = 12.0) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    """
    Cleans a PDB structure by removing water and solvent ions,
    and automatically calculates the binding pocket / grid box coordinates.
    """
    lines = pdb_content.splitlines()
    clean_lines = []
    coords = []
    
    # Common solvent and counter-ion residue names to strip
    SOLVENT_RESIDUES = {
        'HOH', 'WAT', 'SOL', 'TIP3', 'TIP4', 'CL', 'NA', 'MG', 'K', 'CA', 'ZN', 'PO4', 'SO4', 'ACT', 'EDO', 'DMS'
    }

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            res_name = line[17:20].strip()
            if res_name in SOLVENT_RESIDUES:
                continue
            clean_lines.append(line)
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except Exception:
                pass
        elif line.startswith(("TER", "ENDMDL")):
            clean_lines.append(line)

    if not coords:
        return None, None

    coords_arr = np.array(coords)
    center = np.mean(coords_arr, axis=0)
    min_coord = np.min(coords_arr, axis=0)
    max_coord = np.max(coords_arr, axis=0)
    size = (max_coord - min_coord) + padding

    box = {
        'center_x': round(float(center[0]), 3),
        'center_y': round(float(center[1]), 3),
        'center_z': round(float(center[2]), 3),
        'size_x': round(float(min(size[0], 60.0)), 3),
        'size_y': round(float(min(size[1], 60.0)), 3),
        'size_z': round(float(min(size[2], 60.0)), 3),
    }

    clean_pdb_str = "\n".join(clean_lines) + "\nEND\n"
    return clean_pdb_str, box
