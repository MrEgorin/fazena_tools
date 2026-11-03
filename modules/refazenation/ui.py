"""
FAZENA Tools — Refazenation Target Prediction & Reverse Docking
High-coverage target prediction and multi-target AutoDock Vina / Vina-GPU reverse docking.
Supports both ComTarget Combinatorial 3D Proteome Fishing (26,000+ crystallographic targets)
and Curated Bioactive Reference Chemotype Panels (53 targets).

References:
1. Li et al. ComTarget: Small-Molecule Target Prediction with Combinatorial Modeling.
   Pharmaceuticals 2026, 19, 715.
2. Trott & Olson. AutoDock Vina: Improving the speed and accuracy of docking.
   J. Comput. Chem. 2010, 31, 455-461.

Author: Yehor Franchuk (@MrEgorin), FAZENA (fazena.org)
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
from modules.comtarget_engine import (
    predict_comtarget_targets,
    get_target_box,
    get_native_score,
    evaluate_delta_score,
    get_all_comtarget_targets_df
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

# Professional SVG Icons
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
        "Screen small molecules across high-coverage proteomes (ComTarget 26,000+ crystallographic targets & Curated Clinical Chemotypes) to identify polypharmacology, repurposing targets, and execute multi-target AutoDock Vina/Vina-GPU reverse docking."
    )

    tab_pred, tab_dock, tab_db = st.tabs([
        "Target Prediction",
        "Multi-Target Reverse Docking",
        "Proteome Target Database"
    ])

    # Initialize session storage
    if 'fetched_targets' not in st.session_state:
        st.session_state['fetched_targets'] = {}
    if 'staged_boxes' not in st.session_state:
        st.session_state['staged_boxes'] = {}
    if 'staged_native_scores' not in st.session_state:
        st.session_state['staged_native_scores'] = {}

    # ═════════════════════════════════════════════════════════════════
    # TAB 1: TARGET PREDICTION (ComTarget 26,000+ & Curated Bioactive Panel)
    # ═════════════════════════════════════════════════════════════════
    with tab_pred:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        col_l1, col_l2 = st.columns([1, 1])

        with col_l1:
            st.markdown(
                f"<h3 style='display: flex; align-items: center; gap: 10px;'>"
                f"<span style='color: var(--fz-accent);'>{ICON_TARGET}</span> Query Small Molecule</h3>",
                unsafe_allow_html=True
            )
            input_mode = st.radio("Input Mode", ["Paste SMILES", "Upload Molecule (.sdf / .mol)"], horizontal=True, key="pred_mode")

            query_smiles = ""
            if "Paste SMILES" in input_mode:
                query_smiles = st.text_input(
                    "Query SMILES",
                    value="Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
                    key="pred_smi"
                )
            else:
                up_file = st.file_uploader("Upload Molecule", type=['sdf', 'mol'], key="pred_file")
                if up_file:
                    parsed = parse_ligand_file(up_file.getvalue(), up_file.name)
                    if parsed:
                        query_smiles = parsed[0]['SMILES']
                        st.info(f"Loaded: **{parsed[0]['name']}**")

        with col_l2:
            st.markdown(
                "<p style='color: var(--fz-text-muted); font-size: 0.85rem; text-transform: uppercase; letter-spacing: 1px; margin-top: 1rem;'>2D Molecular Topology</p>",
                unsafe_allow_html=True
            )
            if query_smiles:
                render_ligand_2d(query_smiles, width=240, height=200)

        st.markdown("</div>", unsafe_allow_html=True)

        # Engine Selection & Navigation Filters
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        st.markdown("<h4 style='margin-bottom: 0.5rem;'>Target Prediction Engine & Proteome Scope</h4>", unsafe_allow_html=True)

        engine_choice = st.radio(
            "Select Prediction Engine",
            [
                "ComTarget 3D Combinatorial Modeling (26,000+ Crystallographic Targets)",
                "Curated Bioactive Reference Panel (53 Clinical Targets)"
            ],
            index=0,
            key="pred_engine_choice"
        )

        use_comtarget = "ComTarget" in engine_choice

        if use_comtarget:
            st.caption("ComTarget integrates 3D C3DD shape/surface/Hu moments across 26,271 crystallized macromolecular pockets with pre-computed grid boxes and native ligand affinities.")
            col_cf1, col_cf2 = st.columns([2, 1])
            with col_cf1:
                search_query = st.text_input(
                    "Search / Filter Targets (Optional keyword: Gene, Protein Name, UniProt ID, PDB)",
                    placeholder="e.g. Kinase, EGFR, P00533, 1M17, SARS",
                    key="comtarget_search_filter"
                )
            with col_cf2:
                top_k = st.slider("Number of Predicted Targets to Rank", min_value=5, max_value=50, value=15, key="comtarget_topk")
        else:
            st.caption("Curated bioactive panel covers human kinome, GPCRome, nuclear receptors, proteases, ion channels, transporters, and pathogen proteomes via 2D Morgan fingerprints.")
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
                st.error("Please enter a query molecule.")
                return

            with st.spinner("Screening chemical space across proteome target knowledge base..."):
                t_start = time.time()
                if use_comtarget:
                    targets_df = predict_comtarget_targets(
                        query_smiles,
                        top_k=top_k,
                        search_filter=search_query
                    )
                else:
                    targets_df = predict_targets(
                        query_smiles,
                        top_k=top_k,
                        target_class=target_class_choice,
                        organism=organism_choice,
                        min_similarity=min_sim
                    )
                t_elapsed = time.time() - t_start

                st.session_state['predicted_targets'] = targets_df
                st.session_state['predicted_engine'] = "comtarget" if use_comtarget else "curated"
                st.session_state['pred_query_smi'] = query_smiles
                render_success_alert(f"Identified {len(targets_df)} predicted target profiles in {t_elapsed:.2f} seconds.")

        if 'predicted_targets' in st.session_state:
            targets_df = st.session_state['predicted_targets']
            current_engine = st.session_state.get('predicted_engine', 'comtarget')
            st.markdown("### Predicted Target Profile Rankings")

            render_section_divider()

            if targets_df.empty:
                st.info("No targets met the search criteria. Try broadening your query filter.")
            else:
                if current_engine == "comtarget":
                    disp_cols = [
                        'PDB ID', 'Target Name', 'Gene', 'UniProt ID',
                        'ComTarget 3D Similarity', 'Estimated Probability',
                        'Native Affinity (kcal/mol)', 'Pocket Center (X, Y, Z)', 'Pocket Size (X, Y, Z)'
                    ]
                    available_cols = [c for c in disp_cols if c in targets_df.columns]
                    disp_df = targets_df[available_cols].copy()
                    if 'Estimated Probability' in disp_df.columns:
                        disp_df['Estimated Probability'] = disp_df['Estimated Probability'].apply(lambda x: f"{x:.1f}%")
                    st.dataframe(disp_df, use_container_width=True)
                else:
                    disp_cols = ['Target Name', 'Gene', 'Class', 'Organism', 'Probability', 'Tanimoto Similarity', 'UniProt ID', 'PDB ID', 'Reference Compound']
                    disp_df = targets_df[[c for c in disp_cols if c in targets_df.columns]].copy()
                    if 'Probability' in disp_df.columns:
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
                    target_options = []
                    for _, r in targets_df.iterrows():
                        pdb_code = str(r.get('PDB ID', '')).strip().upper()
                        t_name = str(r.get('Target Name', 'Target'))
                        if pdb_code:
                            target_options.append(f"{t_name[:45]} ({pdb_code})")

                    if target_options:
                        selected_option = st.selectbox("Send Predicted Target to Reverse Docking", target_options)
                        if st.button("Stage Target Structure & Pocket Box for Reverse Docking", use_container_width=True):
                            pdb_code = selected_option.split("(")[-1].replace(")", "").strip().upper()
                            with st.spinner(f"Fetching PDB {pdb_code} from RCSB..."):
                                fetched_content = fetch_pdb_from_rcsb(pdb_code)
                                if fetched_content:
                                    st.session_state['fetched_targets'][f"{pdb_code}.pdb"] = fetched_content

                                    # Pre-configure pocket bounding box from box.data
                                    box_info = get_target_box(pdb_code)
                                    if box_info:
                                        st.session_state['staged_boxes'][f"{pdb_code}.pdb"] = box_info

                                    # Pre-configure native ligand baseline score from score_data.txt
                                    score_info = get_native_score(pdb_code)
                                    if score_info:
                                        st.session_state['staged_native_scores'][f"{pdb_code}.pdb"] = score_info

                                    box_desc = (
                                        f"center ({box_info['center_x']:.1f}, {box_info['center_y']:.1f}, {box_info['center_z']:.1f})"
                                        if box_info else "automatic geometry"
                                    )
                                    render_success_alert(
                                        f"Successfully staged target {pdb_code}.pdb with binding pocket {box_desc} for Reverse Docking!"
                                    )
                                else:
                                    st.error(f"Could not retrieve PDB {pdb_code} from RCSB server.")

    # ═════════════════════════════════════════════════════════════════
    # TAB 2: REVERSE DOCKING (AutoDock Vina & Vina-GPU + Delta Score)
    # ═════════════════════════════════════════════════════════════════
    with tab_dock:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        col_d1, col_d2 = st.columns([1, 1])

        with col_d1:
            st.markdown(
                f"<h3 style='display: flex; align-items: center; gap: 10px;'>"
                f"<span style='color: var(--fz-accent);'>{ICON_REVERSE}</span> Ligand for Profiling</h3>",
                unsafe_allow_html=True
            )
            dock_smi = st.text_input(
                "Ligand SMILES",
                value=st.session_state.get('pred_query_smi', "CC(=O)Oc1ccccc1C(=O)O"),
                key="rev_dock_smi"
            )
            if dock_smi:
                render_ligand_2d(dock_smi, width=200, height=180)

        with col_d2:
            st.markdown(
                f"<h3 style='display: flex; align-items: center; gap: 10px;'>"
                f"<span style='color: var(--fz-accent);'>{ICON_3D}</span> Target Macromolecules</h3>",
                unsafe_allow_html=True
            )

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
                                fname = f"{custom_pdb_code.upper()}.pdb"
                                st.session_state['fetched_targets'][fname] = data
                                box_info = get_target_box(custom_pdb_code.upper())
                                if box_info:
                                    st.session_state['staged_boxes'][fname] = box_info
                                score_info = get_native_score(custom_pdb_code.upper())
                                if score_info:
                                    st.session_state['staged_native_scores'][fname] = score_info
                                render_success_alert(f"Fetched & staged {fname}!")
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
                st.error("Please provide a ligand SMILES.")
                return

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
                st.error("Please provide at least one target structure (upload PDB, fetch PDB ID, or select built-in targets).")
                return

            dock_results = []
            progress_bar = st.progress(0)
            status_text = st.empty()

            engine_mode = dock_params.get('engine_mode', 'vina')
            gpu_threads = dock_params.get('gpu_threads', 8000)
            vina_gpu_bin = dock_params.get('vina_gpu_bin')

            for i, (t_name, pdb_content, is_rna) in enumerate(targets_to_dock):
                status_text.text(f"Docking against {t_name} ({i+1}/{len(targets_to_dock)})...")

                # Check if pocket grid box is pre-configured from ComTarget box.data
                staged_box = st.session_state.get('staged_boxes', {}).get(t_name)
                pdb_clean_code = t_name.split(".")[0].split("_")[0].upper()
                if staged_box is None:
                    staged_box = get_target_box(pdb_clean_code)

                if staged_box:
                    t_box = staged_box
                else:
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

                # Evaluate ComTarget Delta Score
                delta_info = evaluate_delta_score(score, pdb_clean_code)

                dock_results.append({
                    'Target Macromolecule': t_name,
                    'Target Type': 'RNA' if is_rna else 'Protein',
                    'Vina Score (kcal/mol)': score,
                    'Native Reference (kcal/mol)': delta_info['native_score'],
                    'Delta Score (kcal/mol)': delta_info['delta_score'],
                    'Binding Likelihood': delta_info['category'],
                    'Likelihood Description': delta_info['description'],
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
            st.markdown("### Target Binding Affinity Profiles & Delta Scoring")

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
                        st.caption("2D Binding Residue Interactions & ComTarget Profile")
                        st.markdown(
                            f"<div style='background: var(--fz-card-bg); border: 1px solid var(--fz-border); border-radius: 8px; padding: 12px; margin-bottom: 12px;'>"
                            f"<div style='font-size: 0.85rem; color: var(--fz-text-muted);'>ComTarget Delta Score Classification</div>"
                            f"<div style='font-size: 1.25rem; font-weight: bold; color: {chosen_data.get('color', 'var(--fz-accent)')};'>{chosen_data.get('Binding Likelihood', 'N/A')}</div>"
                            f"<div style='font-size: 0.85rem;'>{chosen_data.get('Likelihood Description', '')}</div>"
                            f"<div style='margin-top: 6px; font-size: 0.82rem; color: var(--fz-text-muted);'>"
                            f"Vina Affinity: <b>{chosen_data.get('Vina Score (kcal/mol)')} kcal/mol</b> | Native Co-Crystal: <b>{chosen_data.get('Native Reference (kcal/mol)')} kcal/mol</b> | Δ: <b>{chosen_data.get('Delta Score (kcal/mol)')} kcal/mol</b>"
                            f"</div></div>",
                            unsafe_allow_html=True
                        )
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
        st.markdown(
            f"<h3 style='display: flex; align-items: center; gap: 10px;'>"
            f"<span style='color: var(--fz-accent);'>{ICON_DATABASE}</span> Proteome Target Knowledge Base Overview</h3>",
            unsafe_allow_html=True
        )

        db_view = st.radio(
            "Select Database View",
            [
                "ComTarget High-Coverage Proteome (26,000+ Crystallographic Targets)",
                "Curated Bioactive Reference Panel (53 Targets)"
            ],
            horizontal=True,
            key="db_view_radio"
        )

        if "ComTarget" in db_view:
            st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.92rem;'>Comprehensive repository of 26,271 crystallographic targets with annotated UniProt accessions, gene symbols, pre-computed binding pockets, and native co-crystallized drug baseline affinities.</p>", unsafe_allow_html=True)

            col_ct1, col_ct2, col_ct3 = st.columns(3)
            col_ct1.metric("Total Macromolecular Pockets", "26,271")
            col_ct2.metric("Target Annotations", "29,342")
            col_ct3.metric("Baseline Pocket Grids", "26,339")

            col_s1, col_s2 = st.columns([3, 1])
            with col_s1:
                db_filter = st.text_input(
                    "Search Database by Gene, Protein Name, UniProt, or PDB ID",
                    placeholder="e.g. Kinase, CDK, P00533, 1IEP",
                    key="db_comtarget_filter"
                )

            # Lazy-load full database dataframe
            if 'comtarget_full_df' not in st.session_state:
                with st.spinner("Indexing 26,000+ proteome target database..."):
                    st.session_state['comtarget_full_df'] = get_all_comtarget_targets_df()

            full_ct_df = st.session_state['comtarget_full_df']

            if db_filter.strip():
                q = db_filter.strip().lower()
                filtered_ct_df = full_ct_df[
                    full_ct_df['PDB ID'].str.lower().str.contains(q, na=False) |
                    full_ct_df['Protein Name'].str.lower().str.contains(q, na=False) |
                    full_ct_df['Gene'].str.lower().str.contains(q, na=False) |
                    full_ct_df['UniProt ID'].str.lower().str.contains(q, na=False)
                ]
            else:
                filtered_ct_df = full_ct_df

            st.caption(f"Showing **{min(len(filtered_ct_df), 1000)}** of **{len(filtered_ct_df)}** matching targets:")
            st.dataframe(filtered_ct_df.head(1000), use_container_width=True, height=450)

            st.download_button(
                "Download ComTarget Proteome Knowledge Base (.csv)",
                filtered_ct_df.to_csv(index=False).encode('utf-8'),
                "comtarget_proteome_database.csv",
                "text/csv",
                use_container_width=True
            )
        else:
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
            curated_df = pd.DataFrame(rows)

            col_st1, col_st2, col_st3 = st.columns(3)
            col_st1.metric("Total Curated Targets", len(curated_df))
            col_st2.metric("Target Classes", curated_df['Class'].nunique())
            col_st3.metric("Organisms", curated_df['Organism'].nunique())

            st.dataframe(curated_df, use_container_width=True, height=400)

            st.download_button(
                "Download Curated Reference Knowledge Base (.csv)",
                curated_df.to_csv(index=False).encode('utf-8'),
                "curated_target_knowledge_base.csv",
                "text/csv",
                use_container_width=True
            )

        st.markdown("</div>", unsafe_allow_html=True)

    render_footer()
