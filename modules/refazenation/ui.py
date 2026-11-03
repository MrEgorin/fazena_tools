"""
FAZENA Tools — Refazenation Target Prediction & Reverse Docking
High-coverage target prediction and multi-target AutoDock Vina / Vina-GPU reverse docking.
Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

import streamlit as st
import pandas as pd
import numpy as np
import time
import math
import os
import urllib.request
from typing import Optional, Dict, List, Any, Tuple
from rdkit import Chem
from modules.target_prediction import (
    predict_targets,
    get_all_target_classes,
    get_all_organisms,
    TARGET_PROFILES
)
from modules.docking_engine import (
    dock_single_compound,
    prepare_ligand_pdbqt,
    calculate_grid_box,
    analyze_interactions,
    find_vina_gpu_binary
)
from modules.input_parser import parse_ligand_text, parse_ligand_file, clean_pdb_and_calc_box
from modules.visualization import (
    render_ligand_2d,
    get_base64_mol_image,
    render_complex_3d,
    render_receptor_3d,
    render_interaction_summary_ui
)
from modules.ui_components import (
    render_header,
    render_section_divider,
    render_footer,
    render_success_alert,
    render_docking_parameter_controls
)
from modules.config import HARDWARE

# Professional SVG Icons (No playful emojis)
ICON_TARGET = '<svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><path d="M12 8v8M8 12h8"/><circle cx="12" cy="12" r="3"/></svg>'
ICON_REVERSE = '<svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M7 11l5-5 5 5M7 17l5-5 5 5"/></svg>'
ICON_3D = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/><polyline points="3.27 6.96 12 12.01 20.73 6.96"/><line x1="12" y1="22.08" x2="12" y2="12"/></svg>'
ICON_DATABASE = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><ellipse cx="12" cy="5" rx="9" ry="3"/><path d="M21 12c0 1.66-4 3-9 3s-9-1.34-9-3"/><path d="M3 5v14c0 1.66 4 3 9 3s9-1.34 9-3V5"/></svg>'


def fetch_pdb_from_rcsb(pdb_id: str) -> Optional[str]:
    """Fetches high-resolution PDB coordinate text directly from RCSB PDB server."""
    clean_id = pdb_id.strip().upper()
    if not clean_id or len(clean_id) != 4:
        return None
    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'FAZENA-Tools/2.0'})
        with urllib.request.urlopen(req, timeout=12) as response:
            content = response.read().decode('utf-8', errors='ignore')
            if "ATOM" in content or "HETATM" in content:
                return content
    except Exception:
        pass
    return None


def render_page():
    render_header(
        "Refazenation | Target Prediction & Reverse Docking",
        "Screen small molecules across high-coverage human and pathogen proteomes to identify polypharmacology, repurposing targets, and execute multi-target AutoDock Vina/Vina-GPU reverse docking."
    )

    tab_pred, tab_dock, tab_db = st.tabs([
        "Target Prediction",
        "Multi-Target Reverse Docking",
        "Proteome Target Database"
    ])

    # Initialize session storage for fetched targets
    if 'fetched_targets' not in st.session_state:
        st.session_state['fetched_targets'] = {}

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
                render_ligand_2d(query_smiles, width=240, height=200)

        st.markdown("</div>", unsafe_allow_html=True)

        # Filters for navigating large proteomes
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        st.markdown("<h4 style='margin-bottom: 0.5rem;'>Proteome Navigation & Filtering</h4>", unsafe_allow_html=True)
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            target_class_choice = st.selectbox("Filter Target Functional Class", get_all_target_classes(), index=0, key="pred_class")
        with col_f2:
            organism_choice = st.selectbox("Filter Organism Proteome", get_all_organisms(), index=0, key="pred_org")
        with col_f3:
            top_k = st.slider("Number of Predicted Targets to Rank", min_value=5, max_value=30, value=12, key="pred_topk")

        min_sim = st.slider("Minimum Tanimoto Similarity Cutoff", min_value=0.0, max_value=0.8, value=0.05, step=0.05, key="pred_minsim")
        st.markdown("</div>", unsafe_allow_html=True)

        if st.button("Predict Biological Targets", use_container_width=True):
            if not query_smiles:
                st.error("Please enter a query molecule."); return

            with st.spinner("Screening chemical space across proteome target knowledge base..."):
                targets_df = predict_targets(
                    query_smiles,
                    top_k=top_k,
                    target_class=target_class_choice,
                    organism=organism_choice,
                    min_similarity=min_sim
                )
                st.session_state['predicted_targets'] = targets_df
                st.session_state['pred_query_smi'] = query_smiles
                render_success_alert(f"Identified {len(targets_df)} predicted target profiles matching criteria.")

        if 'predicted_targets' in st.session_state:
            targets_df = st.session_state['predicted_targets']
            st.markdown("### Predicted Target Profile Rankings")

            render_section_divider()

            if targets_df.empty:
                st.info("No targets met the minimum similarity cutoff. Try lowering the cutoff threshold.")
            else:
                disp_df = targets_df[['Target Name', 'Gene', 'Class', 'Organism', 'Probability', 'Tanimoto Similarity', 'UniProt ID', 'PDB ID', 'Reference Compound']].copy()
                disp_df['Probability'] = disp_df['Probability'].apply(lambda x: f"{x:.1f}%")
                st.dataframe(disp_df, use_container_width=True)

                render_section_divider()
                col_exp1, col_exp2 = st.columns([1, 1])
                with col_exp1:
                    st.download_button(
                        "Export Target Predictions (.csv)",
                        targets_df.to_csv(index=False).encode('utf-8'),
                        "target_predictions.csv",
                        "text/csv",
                        use_container_width=True
                    )
                with col_exp2:
                    selected_pdb = st.selectbox(
                        "Send Predicted Target to Reverse Docking",
                        [f"{r['Target Name']} ({r['PDB ID']})" for _, r in targets_df.iterrows() if r['PDB ID']]
                    )
                    if st.button("Fetch & Stage Selected Target Structure", use_container_width=True):
                        pdb_code = selected_pdb.split("(")[-1].replace(")", "").strip()
                        with st.spinner(f"Fetching PDB {pdb_code} from RCSB..."):
                            fetched_content = fetch_pdb_from_rcsb(pdb_code)
                            if fetched_content:
                                st.session_state['fetched_targets'][f"{pdb_code}.pdb"] = fetched_content
                                render_success_alert(f"Successfully staged target {pdb_code}.pdb for Reverse Docking!")
                            else:
                                st.error(f"Could not retrieve PDB {pdb_code} from RCSB server.")

    # ═════════════════════════════════════════════════════════════════
    # TAB 2: REVERSE DOCKING (AutoDock Vina & Vina-GPU)
    # ═════════════════════════════════════════════════════════════════
    with tab_dock:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        col_d1, col_d2 = st.columns([1, 1])

        with col_d1:
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_REVERSE}</span> Ligand for Profiling</h3>", unsafe_allow_html=True)
            dock_smi = st.text_input("Ligand SMILES", value=st.session_state.get('pred_query_smi', "CC(=O)Oc1ccccc1C(=O)O"), key="rev_dock_smi")
            if dock_smi:
                render_ligand_2d(dock_smi, width=200, height=180)

        with col_d2:
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_3D}</span> Target Macromolecules</h3>", unsafe_allow_html=True)
            
            # Direct PDB ID fetcher
            col_rcsb1, col_rcsb2 = st.columns([2, 1])
            with col_rcsb1:
                custom_pdb_code = st.text_input("Fetch PDB Structure by ID from RCSB", placeholder="e.g. 1M17, 1IEP, 6LU7, 2RH1", key="fetch_pdb_id")
            with col_rcsb2:
                st.markdown("<div style='margin-top: 28px;'></div>", unsafe_allow_html=True)
                if st.button("Fetch PDB", use_container_width=True):
                    if custom_pdb_code:
                        with st.spinner(f"Downloading {custom_pdb_code.upper()}..."):
                            data = fetch_pdb_from_rcsb(custom_pdb_code)
                            if data:
                                st.session_state['fetched_targets'][f"{custom_pdb_code.upper()}.pdb"] = data
                                render_success_alert(f"Fetched {custom_pdb_code.upper()}.pdb!")
                            else:
                                st.error("Failed to fetch PDB. Please verify the 4-character ID.")

            uploaded_pdbs = st.file_uploader(
                "Upload Custom Protein PDB Files (.pdb / .pdbqt)",
                type=['pdb', 'pdbqt'],
                accept_multiple_files=True,
                key="rev_pdbs"
            )
            use_builtin_targets = st.checkbox(
                "Include Reference Panel (Kinases 9MXC, 9MHF, RNA model)",
                value=not bool(uploaded_pdbs or st.session_state.get('fetched_targets')),
                key="rev_builtin"
            )

            total_available = len(uploaded_pdbs or []) + len(st.session_state.get('fetched_targets', {})) + (3 if use_builtin_targets else 0)
            st.caption(f"Total target structures staged for reverse docking: **{total_available}**")

        st.markdown("</div>", unsafe_allow_html=True)

        render_section_divider()

        # Docking Engine & Parameters
        default_box = {'center_x': 0.0, 'center_y': 0.0, 'center_z': 0.0, 'size_x': 25.0, 'size_y': 25.0, 'size_z': 25.0}
        dock_params = render_docking_parameter_controls(default_box, is_rna=False, key_prefix="rev_dock", allow_gpu=True)

        if st.button("Execute Multi-Target Reverse Docking", use_container_width=True):
            if not dock_smi:
                st.error("Please provide a ligand SMILES."); return
            
            # Prepare targets list: (name, pdb_content, is_rna)
            targets_to_dock = []
            
            # 1. Uploaded PDBs
            if uploaded_pdbs:
                for f in uploaded_pdbs:
                    targets_to_dock.append((f.name, f.getvalue().decode('utf-8', errors='ignore'), False))
            
            # 2. Online Fetched PDBs
            for fname, content in st.session_state.get('fetched_targets', {}).items():
                targets_to_dock.append((fname, content, False))

            # 3. Built-in reference targets
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
                st.error("Please provide at least one target structure (upload PDB, fetch PDB ID, or select built-in targets)."); return

            dock_results = []
            progress_bar = st.progress(0)
            status_text = st.empty()

            engine_mode = dock_params.get('engine_mode', 'vina')
            gpu_threads = dock_params.get('gpu_threads', 8000)
            vina_gpu_bin = dock_params.get('vina_gpu_bin')

            for i, (t_name, pdb_content, is_rna) in enumerate(targets_to_dock):
                status_text.text(f"Docking against {t_name} ({i+1}/{len(targets_to_dock)})...")
                
                # Auto-calculate target box center
                t_box = calculate_grid_box(pdb_content)

                res = dock_single_compound(
                    receptor_pdb=pdb_content,
                    ligand_smiles=dock_smi,
                    compound_name=t_name,
                    grid_box=t_box,
                    exhaustiveness=dock_params['exhaustiveness'],
                    n_poses=dock_params['n_poses'],
                    seed=dock_params['seed'],
                    cpu=dock_params['cpu'],
                    sf_name=dock_params['sf_name'],
                    is_rna=is_rna,
                    engine_mode=engine_mode,
                    gpu_threads=gpu_threads,
                    vina_gpu_bin=vina_gpu_bin
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
                    'Engine': res.get('engine_used', 'AutoDock Vina'),
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
            st.markdown("### Target Binding Affinity Profiles")
            
            res_df = pd.DataFrame(rev_res).drop(columns=['complex_pdb'], errors='ignore')
            st.dataframe(res_df, use_container_width=True)

            render_section_divider()

            # 3D inspection
            st.markdown("#### 3D Complex & 2D Interaction Inspector")
            valid_targets = [r['Target Macromolecule'] for r in rev_res if r.get('Status') == 'Success']
            if valid_targets:
                chosen_target = st.selectbox("Select Target Complex to Inspect", valid_targets)
                chosen_data = next((r for r in rev_res if r['Target Macromolecule'] == chosen_target), None)
                if chosen_data and chosen_data.get('complex_pdb'):
                    col_3d, col_2d = st.columns([3, 2])
                    with col_3d:
                        st.caption(f"Interactive 3D Complex: Ligand in {chosen_target}")
                        render_complex_3d(chosen_data['complex_pdb'], height=450)
                    with col_2d:
                        st.caption("2D Binding Residue Interactions & Topology")
                        inter_summary = analyze_interactions(chosen_data['complex_pdb'])
                        render_interaction_summary_ui(inter_summary, complex_pdb_str=chosen_data['complex_pdb'])
                        st.download_button(
                            "Download Complex PDB",
                            chosen_data['complex_pdb'],
                            f"{chosen_target}_complex.pdb",
                            "chemical/x-pdb",
                            use_container_width=True
                        )

    # ═════════════════════════════════════════════════════════════════
    # TAB 3: PROTEOME TARGET DATABASE
    # ═════════════════════════════════════════════════════════════════
    with tab_db:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_DATABASE}</span> Proteome Target Knowledge Base Overview</h3>", unsafe_allow_html=True)
        st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.92rem;'>Curated reference panels across human kinome, GPCRome, nuclear receptors, proteases, ion channels, transporters, and pathogen proteomes.</p>", unsafe_allow_html=True)

        rows = []
        for t in TARGET_PROFILES:
            comp_names = ", ".join([c[0] for c in t.get('compounds', [])])
            rows.append({
                'Target Name': t['target_name'],
                'Gene': t.get('gene', ''),
                'Class': t['class'],
                'Organism': t.get('organism', 'Human'),
                'UniProt ID': t['uniprot'],
                'PDB ID': t['pdb_id'],
                'Reference Chemotypes': comp_names
            })
        full_db_df = pd.DataFrame(rows)

        col_st1, col_st2, col_st3 = st.columns(3)
        col_st1.metric("Total Curated Targets", len(full_db_df))
        col_st2.metric("Target Classes", full_db_df['Class'].nunique())
        col_st3.metric("Organisms", full_db_df['Organism'].nunique())

        st.dataframe(full_db_df, use_container_width=True, height=400)

        st.download_button(
            "Download Full Proteome Target Knowledge Base (.csv)",
            full_db_df.to_csv(index=False).encode('utf-8'),
            "proteome_target_knowledge_base.csv",
            "text/csv",
            use_container_width=True
        )
        st.markdown("</div>", unsafe_allow_html=True)

    render_footer()
