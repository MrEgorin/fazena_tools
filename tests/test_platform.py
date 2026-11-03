"""
FAZENA Tools — Core Platform Test Suite
Verifies module imports, compatibility classes, deep learning co-embeddings, and docking utilities.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

import os
import sys
from pathlib import Path
import unittest

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

class TestFazenaPlatform(unittest.TestCase):
    def test_imports_and_render_pages(self):
        import app
        import modules.config as config
        import modules.input_parser as input_parser
        import modules.admet as admet
        import modules.docking_engine as docking_engine
        import modules.conplex_screener as conplex_screener
        import modules.rna_screener as rna_screener
        import modules.target_prediction as target_prediction
        import modules.ui_components as ui_components
        import modules.visualization as visualization
        import modules.virtual_screening as rna_vs
        import modules.binding_prediction as rna_bsp
        import modules.prozena as prozena
        import modules.refazenation as refazenation
        import modules.rna_hunter as rna_hunter
        import modules.docking as docking
        import modules.admet_filter as admet_filter

        self.assertTrue(callable(prozena.render_page))
        self.assertTrue(callable(refazenation.render_page))
        self.assertTrue(callable(rna_vs.render_page))
        self.assertTrue(callable(rna_bsp.render_page))
        self.assertTrue(callable(rna_hunter.render_page))
        self.assertTrue(callable(app.render_home_dashboard))
        self.assertTrue(callable(app.main))

    def test_compatibility_classes(self):
        from modules.prozena import ProteinScreening
        from modules.refazenation import DrugRepurposing
        from modules.rna_hunter import RNAHunter
        from modules.docking import AutoDocking
        from modules.admet_filter import ADMEFilter

        adme = ADMEFilter()
        filt_res = adme.filter("CC(=O)Oc1ccccc1C(=O)O")
        self.assertTrue(filt_res['passed'])
        self.assertTrue(filt_res['lipinski_passed'])

        rh = RNAHunter()
        rh_res = rh.predict("AUCGAUCG", ["CC(=O)Oc1ccccc1C(=O)O"])
        self.assertIn("CC(=O)Oc1ccccc1C(=O)O", rh_res)

        dr = DrugRepurposing()
        dr_res = dr.predict("CC(=O)Oc1ccccc1C(=O)O")
        self.assertIn('predicted_targets', dr_res)

        ps = ProteinScreening()
        ps_res = ps.predict("MKTVRQERLKSIVRILGARFTNVTGIQKLVCAIVILIVSLITGGLVLTQALQQQQAEVRRLR", ["CC(=O)Oc1ccccc1C(=O)O"])
        self.assertIn("CC(=O)Oc1ccccc1C(=O)O", ps_res)

    def test_conplex_screening(self):
        from modules.conplex_screener import run_conplex_screening
        protein = "MKTVRQERLKSIVRILGARFTNVTGIQKLVCAIVILIVSLITGGLVLTQALQQQQAEVRRLR"
        ligands = [{'name': 'Aspirin', 'SMILES': 'CC(=O)Oc1ccccc1C(=O)O'}]
        df, is_real = run_conplex_screening(protein, ligands)
        self.assertTrue(is_real)
        self.assertEqual(len(df), 1)
        self.assertGreater(df.iloc[0]['pKi'], 0)

    def test_target_prediction(self):
        from modules.target_prediction import predict_targets
        df = predict_targets("CC(=O)Oc1ccccc1C(=O)O", top_k=5)
        self.assertEqual(len(df), 5)
        self.assertIn('Target Name', df.columns)

    def test_rna_screening(self):
        from modules.rna_screener import screen_rna_library, predict_rna_binding_sites
        seq = "AUCGAUCGAUCG"
        ligands = [{'name': 'Aspirin', 'SMILES': 'CC(=O)Oc1ccccc1C(=O)O'}]
        df = screen_rna_library(seq, ligands)
        self.assertEqual(len(df), 1)
        sites = predict_rna_binding_sites(seq, "CC(=O)Oc1ccccc1C(=O)O")
        self.assertEqual(len(sites), len(seq))

    def test_docking_engine_prep(self):
        from modules.docking_engine import calculate_grid_box, prepare_receptor_pdbqt, prepare_ligand_pdbqt
        pdb_sample = """ATOM      1  N   ASP A   1      20.154  30.450  15.200  1.00 20.00           N
ATOM      2  CA  ASP A   1      21.000  31.000  16.000  1.00 20.00           C
ATOM      3  C   ASP A   1      22.000  32.000  17.000  1.00 20.00           C
ATOM      4  O   ASP A   1      23.000  33.000  18.000  1.00 20.00           O
END
"""
        grid = calculate_grid_box(pdb_sample)
        self.assertIn('center_x', grid)
        self.assertIn('size_x', grid)
        rec = prepare_receptor_pdbqt(pdb_sample)
        self.assertIsNotNone(rec)
        lig = prepare_ligand_pdbqt("CC(=O)Oc1ccccc1C(=O)O")
        self.assertIsNotNone(lig)

if __name__ == '__main__':
    unittest.main()
