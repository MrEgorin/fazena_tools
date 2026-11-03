"""
FAZENA Tools — Consolidated ADMET and Safety Filtering Engine
Provides Lipinski Rule of 5, Veber parameters, RDKit PAINS and Brenk catalogs,
and toxicophore/reactive functional group screening.
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

from typing import Dict, List, Tuple, Any, Optional
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# ═════════════════════════════════════════════════════════════════════
# TOXICOPHORES & REACTIVE GROUPS (SMARTS)
# ═════════════════════════════════════════════════════════════════════
REACTIVE_SMARTS = {
    'Epoxide': '[OX2]1[CX4][CX4]1',
    'Quinone': 'O=C1[#6]~[#6]C(=O)[#6]~[#6]1',
    'Nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',
    'Michael_Acceptor': '[#6]=[#6]-[#6](=[OX1])',
    'Aliphatic_Halide': '[CX4][F,Cl,Br,I]',
    'Acyl_Halide': '[#6](=[OX1])[Cl,Br,I,F]',
    'Acid_Anhydride': '[#6](=[OX1])[OX2][#6](=[OX1])',
    'Alkyl_Sulfonate': '[OX2][SX4](=[OX1])(=[OX1])[CX4]',
    'Peroxide': '[OX2][OX2]',
    'Azide': '[NX1]~[NX2]~[NX1]',
    'Hydrazine': '[NX3][NX3]',
    'Isocyanate': '[NX2]=[CX2]=[OX1]',
    'Phosphite_Ester': '[PX3]([OX2])[OX2][OX2]',
    'Sulfonyl_Halide': '[SX4](=[OX1])(=[OX1])[Cl,Br,I,F]',
}

# Compile SMARTS patterns once for performance
_COMPILED_PATTERNS = {name: Chem.MolFromSmarts(smarts) for name, smarts in REACTIVE_SMARTS.items()}

# ═════════════════════════════════════════════════════════════════════
# RDKit FILTER CATALOGS (PAINS & BRENK)
# ═════════════════════════════════════════════════════════════════════
_PAINS_CATALOG: Optional[FilterCatalog] = None
_BRENK_CATALOG: Optional[FilterCatalog] = None

def _get_pains_catalog() -> FilterCatalog:
    global _PAINS_CATALOG
    if _PAINS_CATALOG is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        _PAINS_CATALOG = FilterCatalog(params)
    return _PAINS_CATALOG

def _get_brenk_catalog() -> FilterCatalog:
    global _BRENK_CATALOG
    if _BRENK_CATALOG is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        _BRENK_CATALOG = FilterCatalog(params)
    return _BRENK_CATALOG

# ═════════════════════════════════════════════════════════════════════
# ADMET CALCULATION & FILTERING
# ═════════════════════════════════════════════════════════════════════

def calculate_physchem_properties(mol: Chem.Mol) -> Dict[str, float]:
    """Calculate standard physicochemical properties for a molecule."""
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    
    # QED (Quantitative Estimate of Drug-likeness)
    try:
        from rdkit.Chem import QED
        qed_val = QED.qed(mol)
    except Exception:
        qed_val = 0.5

    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': int(hbd),
        'HBA': int(hba),
        'TPSA': round(tpsa, 2),
        'RotBonds': int(rotb),
        'HeavyAtoms': int(heavy_atoms),
        'QED': round(qed_val, 3)
    }

def check_reactive_groups(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """Check for reactive groups and toxicophores defined in REACTIVE_SMARTS."""
    matches = []
    for name, pattern in _COMPILED_PATTERNS.items():
        if pattern and mol.HasSubstructMatch(pattern):
            matches.append(name)
    return (len(matches) > 0, matches)

def check_pains(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """Check against RDKit official PAINS filters."""
    catalog = _get_pains_catalog()
    if catalog.HasMatch(mol):
        matches = [entry.GetDescription() for entry in catalog.GetMatches(mol)]
        return True, matches
    return False, []

def check_brenk(mol: Chem.Mol) -> Tuple[bool, List[str]]:
    """Check against Brenk unwanted functional groups."""
    catalog = _get_brenk_catalog()
    if catalog.HasMatch(mol):
        matches = [entry.GetDescription() for entry in catalog.GetMatches(mol)]
        return True, matches
    return False, []

def evaluate_admet(
    smiles_or_mol: Any,
    max_mw: float = 500.0,
    max_logp: float = 5.0,
    max_hbd: int = 5,
    max_hba: int = 10,
    max_tpsa: float = 140.0,
    max_rotb: int = 10,
    filter_pains: bool = True,
    filter_reactive: bool = True,
    filter_brenk: bool = False
) -> Dict[str, Any]:
    """
    Evaluates a single SMILES or RDKit Mol object against comprehensive ADMET & drug-likeness rules.
    Returns property dict, pass/fail booleans, and list of violation reasons.
    """
    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
        smiles = smiles_or_mol
    else:
        mol = smiles_or_mol
        smiles = Chem.MolToSmiles(mol) if mol else ""

    if mol is None:
        return {
            'valid': False,
            'passed': False,
            'SMILES': smiles,
            'reasons': ['Invalid SMILES structure'],
            'properties': {},
            'Tox_Alert': 'Invalid'
        }

    props = calculate_physchem_properties(mol)
    reasons = []

    # 1. Lipinski Rule of 5
    lipinski_violations = 0
    if props['MW'] > max_mw:
        reasons.append(f"MW > {max_mw} ({props['MW']})")
        lipinski_violations += 1
    if props['LogP'] > max_logp:
        reasons.append(f"LogP > {max_logp} ({props['LogP']})")
        lipinski_violations += 1
    if props['HBD'] > max_hbd:
        reasons.append(f"HBD > {max_hbd} ({props['HBD']})")
        lipinski_violations += 1
    if props['HBA'] > max_hba:
        reasons.append(f"HBA > {max_hba} ({props['HBA']})")
        lipinski_violations += 1

    # 2. Veber parameters
    if props['TPSA'] > max_tpsa:
        reasons.append(f"TPSA > {max_tpsa} ({props['TPSA']})")
    if props['RotBonds'] > max_rotb:
        reasons.append(f"RotBonds > {max_rotb} ({props['RotBonds']})")

    # 3. Reactive / Toxicophore checks
    has_reactive, reactive_matches = check_reactive_groups(mol)
    if filter_reactive and has_reactive:
        reasons.append(f"Reactive groups: {', '.join(reactive_matches)}")

    # 4. PAINS check
    has_pains, pains_matches = check_pains(mol)
    if filter_pains and has_pains:
        reasons.append(f"PAINS alert: {', '.join(pains_matches[:2])}")

    # 5. Brenk check (optional strict filter)
    has_brenk, brenk_matches = check_brenk(mol)
    if filter_brenk and has_brenk:
        reasons.append(f"Brenk alert: {', '.join(brenk_matches[:2])}")

    passed = len(reasons) == 0
    tox_alert = reactive_matches[0] if has_reactive else ("PAINS" if has_pains else "None")

    return {
        'valid': True,
        'passed': passed,
        'SMILES': smiles,
        'MW': props['MW'],
        'LogP': props['LogP'],
        'HBD': props['HBD'],
        'HBA': props['HBA'],
        'TPSA': props['TPSA'],
        'RotBonds': props['RotBonds'],
        'QED': props['QED'],
        'lipinski_violations': lipinski_violations,
        'reasons': reasons,
        'has_reactive': has_reactive,
        'reactive_alerts': reactive_matches,
        'has_pains': has_pains,
        'pains_alerts': pains_matches,
        'Tox_Alert': tox_alert
    }

def filter_compound_library(
    smiles_list: List[str],
    max_mw: float = 500.0,
    max_logp: float = 5.0,
    max_hbd: int = 5,
    max_hba: int = 10,
    max_tpsa: float = 140.0,
    filter_pains: bool = True,
    filter_reactive: bool = True
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Screens an entire library of SMILES.
    Returns: (passed_records, rejected_records)
    """
    passed = []
    rejected = []

    for sm in smiles_list:
        res = evaluate_admet(
            sm,
            max_mw=max_mw,
            max_logp=max_logp,
            max_hbd=max_hbd,
            max_hba=max_hba,
            max_tpsa=max_tpsa,
            filter_pains=filter_pains,
            filter_reactive=filter_reactive
        )
        if res.get('passed', False):
            passed.append(res)
        else:
            rejected.append(res)

    return passed, rejected
