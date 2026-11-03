"""
FAZENA Tools — ComTarget 3D Proteome Target Fishing & Reverse Docking Engine
Implements combinatorial modeling for ultra-large scale macromolecular target prediction.
Integrates 3D molecular shape & Hu moment similarity (C3DD descriptors) across 26,000+ crystallographic
pockets, automated binding site grid box retrieval, and delta-score reverse docking evaluation.

Reference:
Li, Y.; Shi, Q.; Lu, X.; Yang, D.; Yeerken, D.; Jin, H.; Sun, Q.
ComTarget: Small-Molecule Target Prediction with Combinatorial Modeling.
Pharmaceuticals 2026, 19, 715. https://doi.org/10.3390/ph19050715

Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

from typing import List, Dict, Any, Optional, Tuple
import os
import ctypes
from pathlib import Path
import math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

_DATA_DIR = Path(__file__).parent / "refazenation" / "data"
_C_DIR = Path(__file__).parent / "refazenation" / "c_engine"
_SO_PATH = _C_DIR / "libc3dd.so"
_C_SRC_PATH = _C_DIR / "c3dd.c"

_C_LIB = None
_DB_LOADED = False
_TARGET_NAMES: List[str] = []
_TARGET_DESCRIPTORS: Optional[np.ndarray] = None
_TARGETS_MAP: Dict[str, Dict[str, str]] = {}
_BOX_MAP: Dict[str, Dict[str, float]] = {}
_SCORE_MAP: Dict[str, Dict[str, float]] = {}


def _get_c_lib():
    """Initializes and returns the shared C library handle for C3DD descriptors, auto-compiling if needed."""
    global _C_LIB
    if _C_LIB is not None:
        return _C_LIB

    # Auto-compile C3DD shared library if source exists but binary is not built
    if not _SO_PATH.exists() and _C_SRC_PATH.exists():
        try:
            import subprocess
            subprocess.run(
                ["gcc", "-O3", "-shared", "-fPIC", str(_C_SRC_PATH), "-o", str(_SO_PATH), "-lm"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except Exception:
            pass

    if _SO_PATH.exists():
        try:
            lib = ctypes.CDLL(str(_SO_PATH))
            lib.compute_c3dd_descriptors.argtypes = [
                ctypes.c_int,
                ctypes.POINTER(ctypes.c_float),
                ctypes.POINTER(ctypes.c_float),
                ctypes.POINTER(ctypes.c_float),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_float)
            ]
            lib.compute_c3dd_descriptors.restype = ctypes.c_int
            _C_LIB = lib
        except Exception:
            _C_LIB = None
    return _C_LIB


def _clean_pdb_id(raw_name: str) -> str:
    """Extracts standardized 4-character PDB code from descriptor/box/score identifier."""
    name = raw_name.strip()
    if name.startswith("#"):
        name = name[1:]
    if name.startswith("L_"):
        name = name[2:]
    return name.split("_")[0].upper()


def load_comtarget_database(force_reload: bool = False):
    """
    Loads and caches the full ComTarget database:
    - 26,271 co-crystallized ligand 3D shape descriptors (discriptor.txt)
    - 29,342 target annotations (targets.tsv)
    - 26,339 binding pocket grid box definitions (box.data)
    - 26,276 native ligand baseline docking scores (score_data.txt)
    """
    global _DB_LOADED, _TARGET_NAMES, _TARGET_DESCRIPTORS, _TARGETS_MAP, _BOX_MAP, _SCORE_MAP
    if _DB_LOADED and not force_reload:
        return

    # 1. Target Annotations (targets.tsv)
    targets_file = _DATA_DIR / "targets.tsv"
    _TARGETS_MAP = {}
    if targets_file.exists():
        try:
            df = pd.read_csv(targets_file, sep="\t", dtype=str)
            for _, row in df.iterrows():
                pdb_id = str(row.get("PDB_ID", "")).strip().upper()
                if pdb_id:
                    _TARGETS_MAP[pdb_id] = {
                        "protein_name": str(row.get("Protein names", "Unknown")),
                        "gene_name": str(row.get("Gene Names", "-")),
                        "uniprot_id": str(row.get("Uniport_ID", "-")),
                    }
        except Exception:
            pass

    # 2. Binding Site Grid Boxes (box.data)
    box_file = _DATA_DIR / "box.data"
    _BOX_MAP = {}
    if box_file.exists():
        try:
            with open(box_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    parts = [p.strip() for p in line.split("#") if p.strip()]
                    if len(parts) >= 7:
                        clean_id = _clean_pdb_id(parts[0])
                        try:
                            _BOX_MAP[clean_id] = {
                                "center_x": float(parts[1]),
                                "center_y": float(parts[2]),
                                "center_z": float(parts[3]),
                                "size_x": float(parts[4]),
                                "size_y": float(parts[5]),
                                "size_z": float(parts[6]),
                            }
                        except ValueError:
                            continue
        except Exception:
            pass

    # 3. Native Ligand Baseline Docking Scores (score_data.txt)
    score_file = _DATA_DIR / "score_data.txt"
    _SCORE_MAP = {}
    if score_file.exists():
        try:
            with open(score_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    parts = [p.strip() for p in line.split("#") if p.strip()]
                    if len(parts) >= 4:
                        clean_id = _clean_pdb_id(parts[0])
                        try:
                            _SCORE_MAP[clean_id] = {
                                "best_score": float(parts[1]),
                                "avg_score": float(parts[2]),
                                "worst_score": float(parts[3]),
                            }
                        except ValueError:
                            continue
        except Exception:
            pass

    # 4. Pocket Descriptors (discriptor.txt)
    desc_file = _DATA_DIR / "discriptor.txt"
    _TARGET_NAMES = []
    desc_list = []
    if desc_file.exists():
        try:
            with open(desc_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    parts = [p.strip() for p in line.split("#") if p.strip()]
                    if len(parts) >= 11:
                        try:
                            vec = [float(x) for x in parts[1:11]]
                            _TARGET_NAMES.append(parts[0])
                            desc_list.append(vec)
                        except ValueError:
                            continue
        except Exception:
            pass

    if desc_list:
        _TARGET_DESCRIPTORS = np.array(desc_list, dtype=np.float32)
    else:
        _TARGET_DESCRIPTORS = np.zeros((0, 10), dtype=np.float32)

    _DB_LOADED = True


def compute_c3dd_descriptors(mol_or_smiles: Any) -> Optional[np.ndarray]:
    """
    Computes the 10-dimensional ComTarget C3DD molecular surface and 3D Hu moment invariants:
    v1: Surface Area (Angstrom^2)
    v2: Volume (Angstrom^3)
    v3: Surface-to-Volume Ratio
    v4-v7: Rotational Invariant 3D Hu Moments (tt1, tt2, tt3, tt4)
    v8: Shape Ratio (ks)
    v9: Sphericity (kst)
    v10: Equivalent Sphere Volume Fitting Ratio (kfit)
    """
    if isinstance(mol_or_smiles, str):
        mol = Chem.MolFromSmiles(mol_or_smiles)
    else:
        mol = mol_or_smiles

    if mol is None:
        return None

    # Prepare 3D conformer using RDKit ETKDGv3 and MMFF94 optimization
    mol_3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_status = AllChem.EmbedMolecule(mol_3d, params)
    if embed_status != 0:
        embed_status = AllChem.EmbedMolecule(mol_3d, useRandomCoords=True)
    if embed_status != 0:
        return None

    try:
        AllChem.MMFFOptimizeMolecule(mol_3d, maxIters=200)
    except Exception:
        pass

    conformer = mol_3d.GetConformer()
    num_atoms = mol_3d.GetNumAtoms()

    x_arr = np.array([conformer.GetAtomPosition(i).x for i in range(num_atoms)], dtype=np.float32)
    y_arr = np.array([conformer.GetAtomPosition(i).y for i in range(num_atoms)], dtype=np.float32)
    z_arr = np.array([conformer.GetAtomPosition(i).z for i in range(num_atoms)], dtype=np.float32)
    atomic_nums = np.array([atom.GetAtomicNum() for atom in mol_3d.GetAtoms()], dtype=np.int32)

    lib = _get_c_lib()
    if lib is not None:
        out_v = np.zeros(10, dtype=np.float32)
        res = lib.compute_c3dd_descriptors(
            ctypes.c_int(num_atoms),
            x_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            y_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            z_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            atomic_nums.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            out_v.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        )
        if res == 0:
            return out_v

    return None


def predict_comtarget_targets(
    smiles_or_mol: Any,
    top_k: int = 15,
    min_similarity: float = 0.0,
    search_filter: Optional[str] = None
) -> pd.DataFrame:
    """
    Screens small molecule across the ComTarget 26,271 pocket database.
    Calculates 3D C3DD shape similarity and returns ranked target macromolecules
    with UniProt, Gene, PDB ID, pocket grid box, and native ligand scores.
    """
    load_comtarget_database()
    if _TARGET_DESCRIPTORS is None or len(_TARGET_DESCRIPTORS) == 0:
        return pd.DataFrame()

    query_v = compute_c3dd_descriptors(smiles_or_mol)
    if query_v is None:
        return pd.DataFrame()

    # ComTarget Vectorized Similarity Computation:
    # ratio_i = query_i / target_i
    # sim_i = 1.0 - abs(1.0 - ratio_i)
    # Target similarity = mean across all 10 dimensions
    with np.errstate(divide='ignore', invalid='ignore'):
        ratios = query_v / _TARGET_DESCRIPTORS
        dim_sims = 1.0 - np.abs(1.0 - ratios)
        # Handle division by zero or NaN values
        dim_sims = np.nan_to_num(dim_sims, nan=0.0, posinf=0.0, neginf=-5.0)
        overall_sims = np.mean(dim_sims, axis=1)

    # Sort descending
    sorted_indices = np.argsort(-overall_sims)

    records = []
    seen_pdbs = set()

    filter_query = search_filter.strip().lower() if search_filter else None

    for idx in sorted_indices:
        raw_name = _TARGET_NAMES[idx]
        pdb_id = _clean_pdb_id(raw_name)
        sim = float(overall_sims[idx])

        # Avoid redundant duplicate poses of same PDB in top hits
        if pdb_id in seen_pdbs:
            continue

        info = _TARGETS_MAP.get(pdb_id, {})
        protein_name = info.get("protein_name", "Macromolecular Target")
        gene_name = info.get("gene_name", "-")
        uniprot_id = info.get("uniprot_id", "-")

        # Apply search filter if specified
        if filter_query:
            matched = (
                filter_query in pdb_id.lower() or
                filter_query in protein_name.lower() or
                filter_query in gene_name.lower() or
                filter_query in uniprot_id.lower()
            )
            if not matched:
                continue

        # Calibrated probability/confidence score (0 - 100%)
        # Pure ratio can be <= 0 for divergent structures; positive ratio indicates high shape matching
        if sim >= 0.8:
            prob = min(99.0, 75.0 + (sim - 0.8) * 120.0)
        elif sim >= 0.0:
            prob = max(10.0, 30.0 + sim * 55.0)
        else:
            prob = max(1.0, 30.0 / (1.0 + abs(sim)))

        box_info = _BOX_MAP.get(pdb_id, None)
        score_info = _SCORE_MAP.get(pdb_id, None)

        box_str = (
            f"({box_info['center_x']:.1f}, {box_info['center_y']:.1f}, {box_info['center_z']:.1f})"
            if box_info else "Automatic"
        )
        box_size_str = (
            f"({box_info['size_x']:.1f}, {box_info['size_y']:.1f}, {box_info['size_z']:.1f})"
            if box_info else "(22.0, 22.0, 22.0)"
        )
        native_score = score_info["best_score"] if score_info else -9.0

        records.append({
            "PDB ID": pdb_id,
            "Target Name": protein_name,
            "Gene": gene_name,
            "UniProt ID": uniprot_id,
            "ComTarget 3D Similarity": round(sim, 3),
            "Estimated Probability": round(prob, 1),
            "Native Affinity (kcal/mol)": native_score,
            "Pocket Center (X, Y, Z)": box_str,
            "Pocket Size (X, Y, Z)": box_size_str,
            "_box": box_info,
            "_score": score_info
        })
        seen_pdbs.add(pdb_id)

        if len(records) >= top_k:
            break

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df.index = df.index + 1
    return df


def get_target_box(pdb_id: str) -> Optional[Dict[str, float]]:
    """Retrieves the pre-calculated pocket bounding box from ComTarget box.data."""
    load_comtarget_database()
    clean_id = _clean_pdb_id(pdb_id)
    return _BOX_MAP.get(clean_id, None)


def get_native_score(pdb_id: str) -> Optional[Dict[str, float]]:
    """Retrieves the baseline crystallographic ligand docking affinity from score_data.txt."""
    load_comtarget_database()
    clean_id = _clean_pdb_id(pdb_id)
    return _SCORE_MAP.get(clean_id, None)


def evaluate_delta_score(docking_affinity: float, pdb_id: str) -> Dict[str, Any]:
    """
    Calculates Delta Score:
    Delta Affinity = Docking Score - Native Baseline Score

    Classification according to ComTarget benchmark:
    - HIGH: Delta <= 0.0 kcal/mol (Equal to or stronger affinity than co-crystallized drug)
    - GOOD: 0.0 < Delta <= 1.0 kcal/mol (Potent binding within 1 kcal/mol of native drug)
    - MEDIUM: 1.0 < Delta <= 2.5 kcal/mol (Moderate binding)
    - LOW: 2.5 < Delta <= 4.0 kcal/mol (Weak binding)
    - VERY LOW: Delta > 4.0 kcal/mol (Unlikely to bind effectively)
    """
    native = get_native_score(pdb_id)
    if native is None:
        # Default baseline if not present
        native_best = -8.5
    else:
        native_best = native.get("best_score", -8.5)

    delta = docking_affinity - native_best

    if delta <= 0.0:
        category = "HIGH"
        desc = "Potent Binder (Equipotent or exceeds co-crystallized reference drug)"
        color = "#10b981"
    elif delta <= 1.0:
        category = "GOOD"
        desc = "Promising Binder (Within 1.0 kcal/mol of native drug)"
        color = "#00B4D8"
    elif delta <= 2.5:
        category = "MEDIUM"
        desc = "Moderate Binder (Marginal affinity difference)"
        color = "#f59e0b"
    elif delta <= 4.0:
        category = "LOW"
        desc = "Weak Binder (Sub-optimal binding relative to native ligand)"
        color = "#ef4444"
    else:
        category = "VERY LOW"
        desc = "Inactive / Sterically Hindered"
        color = "#6b7280"

    return {
        "pdb_id": _clean_pdb_id(pdb_id),
        "docking_score": round(docking_affinity, 2),
        "native_score": round(native_best, 2),
        "delta_score": round(delta, 2),
        "category": category,
        "description": desc,
        "color": color
    }


def get_all_comtarget_targets_df() -> pd.DataFrame:
    """Returns the complete dataset of 26,000+ macromolecular pockets for browsing and export."""
    load_comtarget_database()
    rows = []
    for raw_name in _TARGET_NAMES:
        pdb_id = _clean_pdb_id(raw_name)
        info = _TARGETS_MAP.get(pdb_id, {})
        box = _BOX_MAP.get(pdb_id, {})
        score = _SCORE_MAP.get(pdb_id, {})

        rows.append({
            "Pocket ID": raw_name,
            "PDB ID": pdb_id,
            "Protein Name": info.get("protein_name", "-"),
            "Gene": info.get("gene_name", "-"),
            "UniProt ID": info.get("uniprot_id", "-"),
            "Native Ligand Affinity (kcal/mol)": score.get("best_score", np.nan),
            "Pocket Center X": box.get("center_x", np.nan),
            "Pocket Center Y": box.get("center_y", np.nan),
            "Pocket Center Z": box.get("center_z", np.nan),
            "Pocket Size X": box.get("size_x", np.nan),
            "Pocket Size Y": box.get("size_y", np.nan),
            "Pocket Size Z": box.get("size_z", np.nan),
        })

    df = pd.DataFrame(rows)
    df.index = df.index + 1
    return df
