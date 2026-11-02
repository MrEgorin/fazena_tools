from typing import List, Dict
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd

class ADMEFilter:
    def __init__(self):
        self.pains = self.load_pains()

    def load_pains(self) -> pd.DataFrame:
        # Load PAINS data from a CSV file or other source
        # For demonstration purposes, we'll create a mock DataFrame
        return pd.DataFrame({
            'smarts': ['C1=CC(=O)N2C=C(C=C1)C3=CC=CC=C32',  # Example PAINS pattern
                      'C1=CC(=O)NC1'],
            'name': ['PAINS_Aromatic_Aldehyde', 'PAINS_Aromatic_Ketone']
        })

    def filter(self, smiles: str) -> Dict[str, any]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string")

            # Calculate properties
            mw = Descriptors.MolWt(mol)
            logp = rdMolDescriptors.CalcClogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            tpsa = rdMolDescriptors.CalcTPSA(mol)

            # Check Lipinski's Rule of Five
            lipinski_passed = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)

            # Check PAINS patterns
            pains_passed, pains_reasons = self.check_pains(mol)

            # Overall passed status
            passed = lipinski_passed and pains_passed

            return {
                'passed': passed,
                'lipinski_passed': lipinski_passed,
                'pains_passed': pains_passed,
                'reasons_for_failure': pains_reasons if not passed else [],
                'properties': {
                    'MW': mw,
                    'LogP': logp,
                    'HBD': hbd,
                    'HBA': hba,
                    'TPSA': tpsa
                }
            }

        except Exception as e:
            return {
                'passed': False,
                'lipinski_passed': None,
                'pains_passed': None,
                'reasons_for_failure': [str(e)],
                'properties': {}
            }

    def check_pains(self, mol) -> (bool, List[str]):
        pains_reasons = []
        for index, row in self.pains.iterrows():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(row['smarts'])):
                pains_reasons.append(row['name'])
        return len(pains_reasons) == 0, pains_reasons
