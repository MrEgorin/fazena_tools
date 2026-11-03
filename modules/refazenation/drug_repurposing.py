"""
FAZENA Tools — DrugRepurposing Compatibility & Programmatic API
Provides high-level biological target prediction and drug repurposing pipelines
supporting both Curated Bioactive Panels and ComTarget 26,000+ Proteome Combinatorial Modeling.

Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

from typing import Dict, Any, Optional
from modules.target_prediction import predict_targets
from modules.comtarget_engine import predict_comtarget_targets, get_target_box, evaluate_delta_score

class DrugRepurposing:
    """Programmatic API for biological target prediction and drug repurposing."""
    def __init__(self, default_engine: str = "comtarget"):
        self.default_engine = default_engine

    def predict(self, smiles: str, engine: Optional[str] = None, top_k: int = 15) -> Dict[str, Any]:
        """
        Predicts macromolecular targets and repurposing opportunities for input SMILES.
        
        Args:
            smiles: Small molecule SMILES string.
            engine: 'comtarget' (26,000+ proteome targets) or 'curated' (53 curated targets).
            top_k: Number of candidate targets to return.
            
        Returns:
            Dictionary containing predicted target records, engine metadata, and candidate count.
        """
        chosen_engine = engine or self.default_engine
        if chosen_engine == "comtarget":
            df = predict_comtarget_targets(smiles, top_k=top_k)
        else:
            df = predict_targets(smiles, top_k=top_k)

        return {
            "engine": chosen_engine,
            "query_smiles": smiles,
            "target_count": len(df),
            "predicted_targets": df.to_dict(orient="records")
        }
