"""
FAZENA Tools — ADMEFilter Compatibility Wrapper
Wraps the consolidated ADMET engine in modules/admet.py
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import Dict, Any, List
from rdkit import Chem
from modules.admet import (
    evaluate_admet,
    check_pains,
    check_reactive_groups,
    calculate_physchem_properties
)

class ADMEFilter:
    """Legacy compatibility class wrapping the enhanced ADMET engine."""
    
    def __init__(self):
        pass

    def filter(self, smiles: str) -> Dict[str, Any]:
        """Filters a single SMILES string using full Lipinski, PAINS, and reactive checks."""
        res = evaluate_admet(smiles, filter_pains=True, filter_reactive=True)
        if not res.get('valid', False):
            return {
                'passed': False,
                'lipinski_passed': False,
                'pains_passed': False,
                'reasons_for_failure': res.get('reasons', ['Invalid SMILES']),
                'properties': {}
            }
        
        lipinski_ok = (res.get('lipinski_violations', 0) == 0)
        pains_ok = not res.get('has_pains', False)
        
        return {
            'passed': res.get('passed', False),
            'lipinski_passed': lipinski_ok,
            'pains_passed': pains_ok,
            'reasons_for_failure': res.get('reasons', []),
            'properties': {
                'MW': res.get('MW'),
                'LogP': res.get('LogP'),
                'HBD': res.get('HBD'),
                'HBA': res.get('HBA'),
                'TPSA': res.get('TPSA'),
                'RotBonds': res.get('RotBonds'),
                'QED': res.get('QED')
            }
        }

    def check_pains(self, mol: Chem.Mol) -> (bool, List[str]):
        has_pains, alerts = check_pains(mol)
        return (not has_pains, alerts)
