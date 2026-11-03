"""
FAZENA Tools — Refazenation Module Package
Fast Target Prediction (SwissTargetPrediction alternative) and automated AutoDock Vina reverse docking.
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

from .drug_repurposing import DrugRepurposing
from .ui import render_page

__all__ = [
    "DrugRepurposing",
    "render_page",
]
