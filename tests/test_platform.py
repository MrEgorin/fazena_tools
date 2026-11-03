"""
FAZENA Tools — Core Platform Test Suite
Verifies module imports, compatibility classes, deep learning co-embeddings,
AutoDock Vina RNA/protein docking, flexible parameters, and 2D/3D visualizations.
Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
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

    def test_conplex_screening_and_weights(self):
        from modules.conplex_screener import run_conplex_screening, CHECKPOINT_PATH
        import torch

        # Verify real pre-trained weights
        self.assertTrue(CHECKPOINT_PATH.exists(), f"Weights file missing: {CHECKPOINT_PATH}")
        ckpt = torch.load(CHECKPOINT_PATH, map_location='cpu', weights_only=False)
        self.assertIn('drug_projector.0.weight', ckpt)
        self.assertEqual(ckpt['drug_projector.0.weight'].shape, torch.Size([1024, 2048]))
        self.assertIn('target_projector.0.weight', ckpt)
        self.assertEqual(ckpt['target_projector.0.weight'].shape, torch.Size([1024, 1024]))

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

    def test_docking_engine_prep_and_metals(self):
        from modules.docking_engine import calculate_grid_box, prepare_receptor_pdbqt, prepare_ligand_pdbqt, convert_pose_pdbqt_to_pdb
        pdb_sample = """ATOM      1  N   ASP A   1      20.154  30.450  15.200  1.00 20.00           N
ATOM      2  CA  ASP A   1      21.000  31.000  16.000  1.00 20.00           C
ATOM      3  C   ASP A   1      22.000  32.000  17.000  1.00 20.00           C
ATOM      4  O   ASP A   1      23.000  33.000  18.000  1.00 20.00           O
HETATM    5  MG   MG A   2      24.000  34.000  19.000  1.00 20.00          MG
END
"""
        grid = calculate_grid_box(pdb_sample)
        self.assertIn('center_x', grid)
        self.assertIn('size_x', grid)
        rec = prepare_receptor_pdbqt(pdb_sample)
        self.assertIsNotNone(rec)
        self.assertIn("Mg", rec)  # Must be case-sensitive 'Mg' for AutoDock Vina

        lig, mol = prepare_ligand_pdbqt("CC(=O)Oc1ccccc1C(=O)O")
        self.assertIsNotNone(lig)

        # Verify RNA receptor atom typing with polar hydrogens and metals
        rna_sample = """ATOM      1  P     G A   1       8.808   7.877  12.396  1.00 70.01           P
ATOM      2  C2    G A   1       5.553  13.896  18.828  1.00 47.92           C
ATOM      3  N1    G A   1       4.832  13.684  19.983  1.00 43.51           N
ATOM      4  1H5'  G A   1       7.100   9.000  12.000  1.00 50.00
HETATM    5  MN   MN A   2       8.000  10.000  15.000  1.00 40.00          MN
END
"""
        rec_rna = prepare_receptor_pdbqt(rna_sample, is_rna=True)
        self.assertIsNotNone(rec_rna)
        self.assertIn("A ", rec_rna)   # Aromatic base carbon
        self.assertIn("Mn", rec_rna)  # Case-sensitive Manganese
        # Ensure invalid type '1' is never produced
        self.assertNotIn("+0.000 1", rec_rna)

        # Verify pose conversion produces valid HETATM records with element C instead of A
        pose_sample = """ATOM      1  C   UNL     1       5.731   9.986  18.353  1.00  0.00    +0.145 A 
ENDMDL
"""
        converted_pdb = convert_pose_pdbqt_to_pdb(pose_sample)
        self.assertTrue(converted_pdb.startswith("HETATM"))
        self.assertTrue(converted_pdb.splitlines()[0].endswith(" C"))

    def test_flexible_docking_and_rna_execution(self):
        from modules.docking_engine import dock_single_compound, calculate_grid_box

        # Load reference canonical RNA structure
        rna_path = REPO_ROOT / "modules" / "rna_model.pdb"
        with open(rna_path, 'r') as f:
            rna_pdb = f.read()

        default_box = calculate_grid_box(rna_pdb)
        custom_box = {
            'center_x': default_box['center_x'],
            'center_y': default_box['center_y'],
            'center_z': default_box['center_z'],
            'size_x': 25.0,
            'size_y': 25.0,
            'size_z': 25.0
        }

        # Run real docking with flexible parameters
        res = dock_single_compound(
            receptor_pdb=rna_pdb,
            ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
            compound_name="Aspirin_Test",
            grid_box=custom_box,
            exhaustiveness=4,
            n_poses=2,
            energy_range=2.5,
            seed=42,
            is_rna=True
        )

        self.assertEqual(res.get('status'), 'Success', f"RNA docking failed: {res.get('error')}")
        self.assertLess(res.get('affinity', 0.0), 0.0)
        self.assertIn('complex_pdb', res)
        self.assertEqual(res.get('exhaustiveness'), 4)
        self.assertEqual(res.get('seed'), 42)

    def test_2d_and_3d_interaction_visualization(self):
        import tempfile
        from modules.docking.auto_docking import AutoDocking
        from modules.visualization import (
            InteractionPlots,
            generate_2d_interaction_diagram,
            analyze_interactions
        )

        ad = AutoDocking()
        self.assertTrue(hasattr(ad, 'dock_rna'))
        self.assertTrue(hasattr(ad, 'dock_protein'))
        self.assertTrue(hasattr(ad, 'dock'))

        # Test 2D diagram generation from PDBQT
        ip = InteractionPlots()
        with tempfile.NamedTemporaryFile('w', suffix='.pdbqt', delete=False) as f:
            f.write("ATOM      1  C   UNL     1       5.731   9.986  18.353  1.00  0.00    +0.145 A \nENDMDL\n")
            f_path = f.name
        with tempfile.NamedTemporaryFile('w', suffix='.svg', delete=False) as out:
            out_path = out.name

        ip.generate_2d_interaction_diagram(f_path, out_path)
        with open(out_path) as f:
            svg = f.read()
        self.assertIn("<svg", svg)

        # Test pure python contact analysis on sample complex
        sample_complex = """ATOM      1  P     G A   1       8.808   7.877  12.396  1.00 70.01           P
ATOM      2  OP1   G A   1       9.839   8.021  13.435  1.00 68.96           O
HETATM    3  O   UNL     1       8.900   8.100  12.500  1.00  0.00           O
END
"""
        inter = analyze_interactions(sample_complex)
        self.assertGreater(inter.get('total_interactions', 0), 0)
        self.assertIn('Hydrogen Bond', inter.get('types', {}))

        # Test 2D diagram HTML generator with contact badges
        diag_html = generate_2d_interaction_diagram(sample_complex)
        self.assertIn("<svg", diag_html)
        self.assertIn("G1", diag_html)

if __name__ == '__main__':
    unittest.main()
