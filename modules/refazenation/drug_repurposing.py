"""
FAZENA Tools — DrugRepurposing Compatibility Class
Wraps target prediction engine in modules/target_prediction.py
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

from typing import Dict, Any
from modules.target_prediction import predict_targets

class DrugRepurposing:
    """Programmatic API for biological target prediction and drug repurposing."""
    def __init__(self):
        pass

    def predict(self, smiles: str) -> Dict[str, Any]:
        """Predicts macromolecular targets and repurposing opportunities for input SMILES."""
        df = predict_targets(smiles, top_k=10)
        return {
            'predicted_targets': df.to_dict(orient='records')
        }
