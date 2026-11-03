"""
FAZENA Tools — Refazenation Target Prediction & Reverse Docking
Fast Target Prediction (SwissTargetPrediction alternative) and automated AutoDock Vina reverse docking.
Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

import streamlit as st
import pandas as pd
import numpy as np
import time
import math
import os
from rdkit import Chem
from modules.target_prediction import predict_targets
from modules.docking_engine import (
    dock_single_compound,
    prepare_ligand_pdbqt,
    calculate_grid_box,
    analyze_interactions
)
from modules.input_parser import parse_ligand_text, parse_ligand_file, clean_pdb_and_calc_box
from modules.visualization import (
    render_ligand_2d,
    get_base64_mol_image,
    render_complex_3d,
    render_interaction_summary_ui
)
from modules.ui_components import (
    render_header,
    render_section_divider,
    render_footer,
    render_success_alert
)

# Icons
ICON_TARGET = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><path d="M12 8v8M8 12h8"/><circle cx="12" cy="12" r="3"/></svg>'
ICON_REVERSE = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M7 11l5-5 5 5M7 17l5-5 5 5"/></svg>'
ICON_3D = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/><polyline points="3.27 6.96 12 12.01 20.73 6.96"/><line x1="12" y1="22.08" x2="12" y2="12"/></svg>'

def render_page():
    render_header(
        "Refazenation | Target Prediction & Reverse Docking",
        "Predict macromolecular targets and repurposing opportunities using lightweight fingerprint intelligence and multi-target AutoDock Vina reverse docking."
    )

    tab_pred, tab_dock = st.tabs(["🎯 Fast Target Prediction", "🔄 Multi-Target Reverse Docking"])

    # ═════════════════════════════════════════════════════════════════
    # TAB 1: FAST TARGET PREDICTION (SwissTargetPrediction alternative)
    # ═════════════════════════════════════════════════════════════════
    with tab_pred:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        col_l1, col_l2 = st.columns([1, 1])

        with col_l1:
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_TARGET}</span> Query Small Molecule</h3>", unsafe_allow_html=True)
            input_mode = st.radio("Input Mode", ["Paste SMILES", "Upload Molecule (.sdf / .mol)"], horizontal=True, key="pred_mode")
            
            query_smiles = ""
            if "Paste SMILES" in input_mode:
                query_smiles = st.text_input("Query SMILES", value="Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", key="pred_smi")
            else:
                up_file = st.file_uploader("Upload Molecule", type=['sdf', 'mol'], key="pred_file")
                if up_file:
                    parsed = parse_ligand_file(up_file.getvalue(), up_file.name)
                    if parsed:
                        query_smiles = parsed[0]['SMILES']
                        st.info(f"Loaded: **{parsed[0]['name']}**")

        with col_l2:
            st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.85rem; text-transform: uppercase; letter-spacing: 1px; margin-top: 1rem;'>2D Molecular Topology</p>", unsafe_allow_html=True)
            if query_smiles:
                render_ligand_2d(query_smiles, width=220, height=220)

        st.markdown("</div>", unsafe_allow_html=True)

        top_k = st.slider("Number of Predicted Targets to Display", min_value=3, max_value=15, value=8)

        if st.button("🔍 Predict Biological Targets", use_container_width=True):
            if not query_smiles:
                st.error("Please enter a query molecule."); return

            with st.spinner("Screening chemical space against therapeutic target knowledge base..."):
                targets_df = predict_targets(query_smiles, top_k=top_k)
                st.session_state['predicted_targets'] = targets_df
                st.session_state['pred_query_smi'] = query_smiles
                render_success_alert(f"Identified {len(targets_df)} predicted target profiles!")

        if 'predicted_targets' in st.session_state:
            targets_df = st.session_state['predicted_targets']
            st.markdown("### 🏆 Predicted Target Profile Rankings")

            render_section_divider()
            
            # Formatted table
            disp_df = targets_df[['Target Name', 'Class', 'Probability', 'Tanimoto Similarity', 'UniProt ID', 'PDB ID', 'Reference Compound']].copy()
            disp_df['Probability'] = disp_df['Probability'].apply(lambda x: f"{x:.1f}%")
            st.dataframe(disp_df, use_container_width=True)

            render_section_divider()
            st.download_button(
                "📥 Export Target Predictions (.csv)",
                targets_df.to_csv(index=False).encode('utf-8'),
                "target_predictions.csv",
                "text/csv",
                use_container_width=True
            )

    # ═════════════════════════════════════════════════════════════════
    # TAB 2: REVERSE DOCKING (AutoDock Vina)
    # ═════════════════════════════════════════════════════════════════
    with tab_dock:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        col_d1, col_d2 = st.columns([1, 1])

        with col_d1:
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_REVERSE}</span> Ligand for Profiling</h3>", unsafe_allow_html=True)
            dock_smi = st.text_input("Ligand SMILES", value=st.session_state.get('pred_query_smi', "CC(=O)Oc1ccccc1C(=O)O"), key="rev_dock_smi")
            if dock_smi:
                render_ligand_2d(dock_smi, width=180, height=180)

        with col_d2:
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_3D}</span> Target Structures</h3>", unsafe_allow_html=True)
            uploaded_pdbs = st.file_uploader(
                "Upload Multiple Protein PDB Structures (.pdb / .pdbqt)",
                type=['pdb', 'pdbqt'],
                accept_multiple_files=True,
                key="rev_pdbs"
            )
            use_builtin_targets = st.checkbox("Include built-in reference target panel (Kinase 9MXC, 9MHF, RNA model)", value=not bool(uploaded_pdbs), key="rev_builtin")
            if uploaded_pdbs:
                st.info(f"Loaded **{len(uploaded_pdbs)}** custom target structures.")
            elif use_builtin_targets:
                st.info("Loaded **3** reference target structures for reverse docking.")

        st.markdown("</div>", unsafe_allow_html=True)

        render_section_divider()

        with st.expander("⚙️ Reverse Docking Parameters", expanded=False):
            col_rp1, col_rp2, col_rp3 = st.columns(3)
            with col_rp1:
                rev_exhaustiveness = st.slider("Exhaustiveness", min_value=1, max_value=24, value=6, step=1, key="rev_ex")
            with col_rp2:
                rev_n_poses = st.slider("Poses per Target", min_value=1, max_value=10, value=3, step=1, key="rev_nposes")
            with col_rp3:
                rev_seed = st.number_input("Random Seed", min_value=0, max_value=999999, value=42, step=1, key="rev_seed")

        if st.button("⚡ Execute Multi-Target Reverse Docking", use_container_width=True):
            if not dock_smi:
                st.error("Please provide a ligand SMILES."); return
            
            # Prepare targets list
            targets_to_dock = []
            if uploaded_pdbs:
                for f in uploaded_pdbs:
                    targets_to_dock.append((f.name, f.getvalue().decode('utf-8', errors='ignore'), False))
            
            if use_builtin_targets or not targets_to_dock:
                base_dir = os.path.dirname(os.path.dirname(__file__))
                b1 = os.path.join(base_dir, "tools", "AutoRevDock", "examples", "pro", "9MXC_A.pdb")
                b2 = os.path.join(base_dir, "tools", "AutoRevDock", "examples", "pro", "9MHF_A.pdb")
                b3 = os.path.join(base_dir, "rna_model.pdb")
                for path, is_rna in [(b1, False), (b2, False), (b3, True)]:
                    if os.path.exists(path):
                        with open(path, 'r') as fp:
                            targets_to_dock.append((os.path.basename(path), fp.read(), is_rna))

            if not targets_to_dock:
                st.error("Please upload one or more target PDB files or select built-in targets."); return

            dock_results = []
            progress_bar = st.progress(0)
            status_text = st.empty()

            for i, (t_name, pdb_content, is_rna) in enumerate(targets_to_dock):
                status_text.text(f"Docking against {t_name} ({i+1}/{len(targets_to_dock)})...")
                
                res = dock_single_compound(
                    receptor_pdb=pdb_content,
                    ligand_smiles=dock_smi,
                    compound_name=t_name,
                    exhaustiveness=rev_exhaustiveness,
                    n_poses=rev_n_poses,
                    seed=rev_seed,
                    is_rna=is_rna
                )
                
                score = res.get('affinity', 0.0)
                ka = math.exp(-score / (1.987e-3 * 298.15)) if score != 0.0 else 0.0
                logka = round(math.log10(max(ka, 1.0)), 2) if ka > 0 else 0.0

                dock_results.append({
                    'Target Macromolecule': t_name,
                    'Target Type': 'RNA' if is_rna else 'Protein',
                    'Vina Score (kcal/mol)': score,
                    'Estimated Ki': res.get('ki_estimated', 'N/A'),
                    'Ka (M⁻¹)': f"{ka:.2e}",
                    'logKa': logka,
                    'Status': res.get('status'),
                    'complex_pdb': res.get('complex_pdb', '')
                })
                progress_bar.progress((i + 1) / len(targets_to_dock))

            progress_bar.empty()
            status_text.empty()
            dock_results.sort(key=lambda x: x.get('Vina Score (kcal/mol)', 0.0))
            st.session_state['reverse_dock_results'] = dock_results
            render_success_alert(f"Reverse Docking Complete across {len(targets_to_dock)} targets!")

        if 'reverse_dock_results' in st.session_state:
            rev_res = st.session_state['reverse_dock_results']
            st.markdown("### 📊 Target Binding Affinity Profiles")
            
            res_df = pd.DataFrame(rev_res).drop(columns=['complex_pdb'], errors='ignore')
            st.dataframe(res_df, use_container_width=True)

            render_section_divider()

            # 3D inspection
            st.markdown("#### 🔍 3D Complex & 2D Interaction Inspector")
            valid_targets = [r['Target Macromolecule'] for r in rev_res if r.get('Status') == 'Success']
            if valid_targets:
                chosen_target = st.selectbox("Select Target Complex to Inspect", valid_targets)
                chosen_data = next((r for r in rev_res if r['Target Macromolecule'] == chosen_target), None)
                if chosen_data and chosen_data.get('complex_pdb'):
                    col_3d, col_2d = st.columns([3, 2])
                    with col_3d:
                        st.caption(f"3D Complex: Ligand in {chosen_target}")
                        render_complex_3d(chosen_data['complex_pdb'], height=450)
                    with col_2d:
                        st.caption("2D Binding Residue Interactions & Topology")
                        inter_summary = analyze_interactions(chosen_data['complex_pdb'])
                        render_interaction_summary_ui(inter_summary, complex_pdb_str=chosen_data['complex_pdb'])
                        st.download_button(
                            "📥 Download Complex PDB",
                            chosen_data['complex_pdb'],
                            f"{chosen_target}_complex.pdb",
                            "chemical/x-pdb",
                            use_container_width=True
                        )

    render_footer()
