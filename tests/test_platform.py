"""
FAZENA Tools — Core Platform Test Suite
Verifies module imports, compatibility classes, deep learning co-embeddings,
AutoDock Vina RNA/protein docking, flexible parameters, and 2D/3D visualizations.
Author: Yehor Franchuk (@MrEgorin), FAZENA (fazena.org)
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
        self.assertIn("<img", diag_html)

    def test_vina_gpu_support_and_discovery(self):
        from modules.docking_engine import find_vina_gpu_binary, dock_single_compound
        from modules.config import HARDWARE

        # Check binary discovery function
        gpu_bin = find_vina_gpu_binary()
        self.assertTrue(gpu_bin is None or os.path.exists(gpu_bin))

        # Check hardware detection detects the RTX 5070 Ti or GPU presence
        self.assertIn('has_cuda', HARDWARE)
        self.assertIn('cuda_device_name', HARDWARE)

        # Test single compound docking with engine_mode='vina_gpu' (with automatic graceful fallback)
        pdb_sample = """ATOM      1  N   ASP A   1      20.154  30.450  15.200  1.00 20.00           N
ATOM      2  CA  ASP A   1      21.000  31.000  16.000  1.00 20.00           C
ATOM      3  C   ASP A   1      22.000  32.000  17.000  1.00 20.00           C
ATOM      4  O   ASP A   1      23.000  33.000  18.000  1.00 20.00           O
END
"""
        res = dock_single_compound(
            receptor_pdb=pdb_sample,
            ligand_smiles="CC(=O)Oc1ccccc1C(=O)O",
            exhaustiveness=1,
            n_poses=1,
            engine_mode="vina_gpu",
            gpu_threads=8000
        )
        self.assertIn('status', res)
        if res.get('status') == 'Success':
            self.assertIn('engine_used', res)

    def test_expanded_target_proteome_database(self):
        from modules.target_prediction import (
            TARGET_PROFILES,
            get_all_target_classes,
            get_all_organisms,
            predict_targets
        )

        self.assertGreaterEqual(len(TARGET_PROFILES), 50)
        classes = get_all_target_classes()
        self.assertIn("Kinase", classes)
        self.assertIn("GPCR", classes)
        self.assertIn("Protease", classes)
        self.assertIn("Ion Channel", classes)

        # Test prediction with class filtering
        kinase_hits = predict_targets("Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", top_k=5, target_class="Kinase")
        self.assertFalse(kinase_hits.empty)
        self.assertTrue((kinase_hits['Class'] == 'Kinase').all())

        # Test prediction with organism filtering
        viral_hits = predict_targets("CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C", top_k=5, organism="Viral")
        self.assertFalse(viral_hits.empty)
        self.assertTrue((viral_hits['Organism'] == 'Viral').all())

    def test_comtarget_engine_and_combinatorial_docking(self):
        from modules.comtarget_engine import (
            compute_c3dd_descriptors,
            predict_comtarget_targets,
            get_target_box,
            get_native_score,
            evaluate_delta_score,
            get_all_comtarget_targets_df
        )
        from modules.refazenation import DrugRepurposing

        # 1. Test 3D C3DD Descriptors Calculation
        aspirin_smi = "CC(=O)Oc1ccccc1C(=O)O"
        desc = compute_c3dd_descriptors(aspirin_smi)
        self.assertIsNotNone(desc)
        self.assertEqual(len(desc), 10)
        self.assertGreater(desc[0], 50.0)  # Surface area > 50 A^2
        self.assertGreater(desc[1], 100.0) # Volume > 100 A^3

        # 2. Test 26,000+ Target Prediction Screening
        df_pred = predict_comtarget_targets(aspirin_smi, top_k=5)
        self.assertEqual(len(df_pred), 5)
        self.assertIn("PDB ID", df_pred.columns)
        self.assertIn("ComTarget 3D Similarity", df_pred.columns)
        self.assertIn("Pocket Center (X, Y, Z)", df_pred.columns)

        # 3. Test Grid Box Retrieval
        box_1iep = get_target_box("1IEP")
        self.assertIsNotNone(box_1iep)
        self.assertIn("center_x", box_1iep)
        self.assertIn("size_x", box_1iep)

        # 4. Test Native Baseline Score & Delta Score
        score_1iep = get_native_score("1IEP")
        self.assertIsNotNone(score_1iep)
        self.assertIn("best_score", score_1iep)

        delta_eval = evaluate_delta_score(-14.5, "1IEP")
        self.assertEqual(delta_eval["category"], "HIGH")
        self.assertLessEqual(delta_eval["delta_score"], 0.0)

        # 5. Test Full Database Retrieval
        full_df = get_all_comtarget_targets_df()
        self.assertGreaterEqual(len(full_df), 25000)

        # 6. Test DrugRepurposing programmatic API with ComTarget
        dr = DrugRepurposing(default_engine="comtarget")
        res = dr.predict(aspirin_smi, top_k=3)
        self.assertEqual(res["engine"], "comtarget")
        self.assertEqual(len(res["predicted_targets"]), 3)

if __name__ == '__main__':
    unittest.main()
