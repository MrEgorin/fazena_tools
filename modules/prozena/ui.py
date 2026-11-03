"""
FAZENA Tools — Prozena Protein Screening & GPU Docking
High-throughput protein-ligand screening using ConPLex neural network with AutoDock Vina GPU docking.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

import streamlit as st
import pandas as pd
import numpy as np
import time
import os
from rdkit import Chem
from modules.conplex_screener import run_conplex_screening
from modules.admet import filter_compound_library, evaluate_admet
from modules.input_parser import parse_ligand_text, parse_ligand_file, parse_sequence_input
from modules.docking_engine import dock_top_candidates, analyze_interactions, calculate_grid_box
from modules.visualization import (
    render_ligand_2d,
    get_base64_mol_image,
    render_complex_3d,
    render_receptor_3d,
    render_complex_2d,
    render_interaction_summary_ui
)
from modules.ui_components import (
    render_header,
    render_section_divider,
    render_footer,
    render_success_alert,
    render_docking_parameter_controls
)
from modules.config import HARDWARE, DOCKING_CONFIG

# Icons
ICON_PROZENA = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M20 21v-2a4 4 0 0 0-4-4H8a4 4 0 0 0-4 4v2"/><circle cx="12" cy="7" r="4"/></svg>'
ICON_BIND = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M12 2l3 9h9l-7 5 3 9-8-7-8 7 3-9-7-5h9z"/></svg>'
ICON_EXPORT = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>'
ICON_DOCK = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 2v20M2 12h20"/><circle cx="12" cy="12" r="4"/></svg>'

def render_page():
    render_header(
        "Prozena | Protein Targeted Screening & Docking",
        "Genome-scale protein-ligand screening using ConPLex contrastive co-embedding with automated AutoDock Vina docking."
    )

    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)

    with col1:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_PROZENA}</span> Protein Target</h3>", unsafe_allow_html=True)
        default_protein = "MKTVRQERLKSIVRILGARFTNVTGIQKLVCAIVILIVSLITGGLVLTQALQQQQAEVRRLR"
        protein_input = st.text_area("Protein Sequence (FASTA or Raw)", height=140, value=default_protein)
        clean_protein_seq, prot_name = parse_sequence_input(protein_input, seq_type='protein')
        if clean_protein_seq:
            st.caption(f"Parsed Target Sequence: **{len(clean_protein_seq)} amino acids** ({prot_name})")

        protein_pdb_file = st.file_uploader("Upload 3D Protein Structure (.pdb / .pdbqt) for Docking", type=['pdb', 'pdbqt'])

        default_protein_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tools", "AutoRevDock", "examples", "pro", "9MXC_A.pdb")
        target_preview_content = None
        if protein_pdb_file:
            target_preview_content = protein_pdb_file.getvalue().decode('utf-8', errors='ignore')
            st.caption(f"Loaded custom Protein: **{protein_pdb_file.name}**")
        elif os.path.exists(default_protein_path):
            with open(default_protein_path, 'r') as f:
                target_preview_content = f.read()

        if target_preview_content:
            with st.expander("Target Protein 3D Structure Preview", expanded=False):
                st.caption("3D macromolecular view of target protein architecture:")
                render_receptor_3d(target_preview_content, height=320)

    with col2:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_BIND}</span> Chemical Library</h3>", unsafe_allow_html=True)
        lib_source = st.radio("Library Source", ["Paste SMILES (Multi-line / Library)", "Upload Library File (.sdf, .csv, .tsv, .smi)"], horizontal=True)

        parsed_ligands = []
        if "Paste SMILES" in lib_source:
            default_smiles = (
                "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1 Imatinib\n"
                "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
                "Cn1cnc2c1c(=O)n(C)c(=O)n2C Caffeine\n"
                "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1 Gefitinib\n"
                "CCN1CCN(CC1)C(=O)C2=CC=C(C=C2)NC3=NC=C4C(=N3)N(N=C4)C5CCCC5 Dinaciclib"
            )
            smiles_text = st.text_area("SMILES (one per line, with optional name)", height=140, value=default_smiles)
            parsed_ligands = parse_ligand_text(smiles_text)
        else:
            uploaded_file = st.file_uploader("Upload Compound Library", type=['sdf', 'mol', 'smi', 'smiles', 'csv', 'tsv', 'txt'])
            if uploaded_file:
                parsed_ligands = parse_ligand_file(uploaded_file.getvalue(), uploaded_file.name)

        if parsed_ligands:
            st.info(f"Loaded **{len(parsed_ligands)}** valid molecules in library.")

    st.markdown("</div>", unsafe_allow_html=True)

    render_section_divider()

    # ADMET & Safety Filters
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    with st.expander("Advanced ADMET & Safety Filters (Lipinski, PAINS, Toxicophores)", expanded=True):
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            max_mw = st.number_input("Max Molecular Weight (MW)", value=550.0, step=10.0)
            max_logp = st.number_input("Max LogP", value=5.0, step=0.5)
        with col_f2:
            max_hbd = st.number_input("Max H-Bond Donors", value=5, step=1)
            max_hba = st.number_input("Max H-Bond Acceptors", value=10, step=1)
        with col_f3:
            max_tpsa = st.number_input("Max TPSA (Å²)", value=140.0, step=10.0)
            filter_tox = st.checkbox("Filter Reactive Groups & Toxicophores", value=True)
            filter_pains = st.checkbox("Filter RDKit PAINS Alerts", value=True)
    st.markdown("</div>", unsafe_allow_html=True)

    render_section_divider()

    if st.button("Execute ConPLex Screening Engine", use_container_width=True):
        if not clean_protein_seq or not parsed_ligands:
            st.error("Please provide both protein sequence and chemical library."); return

        with st.spinner("Processing Protein Language & Molecular Co-Embeddings..."):
            results_df, is_real = run_conplex_screening(clean_protein_seq, parsed_ligands)
            st.session_state['prozena_results'] = results_df
            st.session_state['prozena_is_real'] = is_real
            render_success_alert(f"ConPLex Screening Complete! Evaluated {len(results_df)} candidate pairs.")

    if 'prozena_results' in st.session_state:
        results_df = st.session_state['prozena_results']

        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_BIND}</span> Screening Results</h3>", unsafe_allow_html=True)
        
        # Apply filters
        filtered_df = results_df.copy()
        if filter_tox:
            filtered_df = filtered_df[filtered_df['Tox_Alert'] == 'None']
        filtered_df = filtered_df[
            (filtered_df['MW'] <= max_mw) &
            (filtered_df['LogP'] <= max_logp) &
            (filtered_df['HBD'] <= max_hbd) &
            (filtered_df['HBA'] <= max_hba) &
            (filtered_df['TPSA'] <= max_tpsa)
        ].reset_index(drop=True)

        col_m1, col_m2, col_m3, col_m4 = st.columns(4)
        col_m1.metric("Molecules Screened", len(results_df))
        col_m2.metric("Passed ADMET", len(filtered_df))
        top_pki = filtered_df['pKi'].max() if not filtered_df.empty else 0.0
        col_m3.metric("Top pKi", f"{top_pki:.2f}" if top_pki else "N/A")
        mean_pki = filtered_df['pKi'].mean() if not filtered_df.empty else 0.0
        col_m4.metric("Mean pKi", f"{mean_pki:.2f}" if mean_pki else "N/A")

        if not filtered_df.empty:
            filtered_df['Structure'] = filtered_df['SMILES'].apply(get_base64_mol_image)

            html_table = "<div class='results-table-wrap'><table class='results-table'><thead><tr><th>Rank</th><th>Structure</th><th>Name / SMILES</th><th>pKi</th><th>Binding Score</th><th>ADMET Safety</th></tr></thead><tbody>"
            for idx, row in filtered_df.head(25).iterrows():
                bg = "rgba(255,255,255,0.03)" if idx % 2 == 0 else "transparent"
                img_tag = f"<img src='{row['Structure']}' width='110'/>" if row['Structure'] else "N/A"
                admet_str = f"<b>MW:</b> {row['MW']} | <b>LogP:</b> {row['LogP']}<br><b>HBD/HBA:</b> {row['HBD']}/{row['HBA']} | <b>TPSA:</b> {row['TPSA']}"
                html_table += f"<tr style='background-color:{bg}'><td><span class='rank-badge'>{idx+1}</span></td><td>{img_tag}</td><td><b>{row['name']}</b><br><span class='smiles-cell'>{row['SMILES']}</span></td><td class='score-cell' style='color:var(--fz-accent); font-weight:bold; font-size:1.15rem;'>{row['pKi']}</td><td>{row['Binding_Score']:.4f}</td><td style='font-size:0.85rem;'>{admet_str}</td></tr>"
            html_table += "</tbody></table></div>"
            st.markdown(html_table, unsafe_allow_html=True)

            render_section_divider()
            col_ex1, col_ex2 = st.columns(2)
            with col_ex1:
                st.download_button("Export All Results (.csv)", results_df.to_csv(index=False).encode('utf-8'), "prozena_all.csv", "text/csv", use_container_width=True)
            with col_ex2:
                st.download_button("Export Passed Filter Hits (.csv)", filtered_df.drop(columns=['Structure'], errors='ignore').to_csv(index=False).encode('utf-8'), "prozena_passed.csv", "text/csv", use_container_width=True)

            # ═════════════════════════════════════════════════════════════
            # AUTOMATED TOP-5 DOCKING SECTION WITH FULL PARAMETER TUNING
            # ═════════════════════════════════════════════════════════════
            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_DOCK}</span> Automated Top-5 AutoDock Vina Protein Docking</h3>", unsafe_allow_html=True)
            st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.92rem;'>AutoDock Vina performs physics-based docking with full grid exploration and energetic pose ranking for the top hits.</p>", unsafe_allow_html=True)

            default_protein_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tools", "AutoRevDock", "examples", "pro", "9MXC_A.pdb")
            if protein_pdb_file:
                target_pdb_content = protein_pdb_file.getvalue().decode('utf-8', errors='ignore')
                st.info(f"Using uploaded Protein Structure: **{protein_pdb_file.name}**")
            elif os.path.exists(default_protein_path):
                with open(default_protein_path, 'r') as f:
                    target_pdb_content = f.read()
                st.info("Using reference 3D Protein Structure (`9MXC_A.pdb`) for automated docking. (You can also upload a custom .pdb above).")
            else:
                target_pdb_content = None
                st.warning("Please upload a 3D Protein Structure (.pdb) above to run AutoDock Vina docking.")

            dock_params = None
            if target_pdb_content:
                default_box = calculate_grid_box(target_pdb_content)
                dock_params = render_docking_parameter_controls(default_box, is_rna=False, key_prefix="prozena_dock", allow_gpu=True)

                col_btn, col_info = st.columns([1, 2])
                with col_btn:
                    run_dock_btn = st.button("Run AutoDock Vina on Top-5 Compounds", use_container_width=True)
                with col_info:
                    st.caption(f"Search Box Center: `({dock_params['grid_box']['center_x']}, {dock_params['grid_box']['center_y']}, {dock_params['grid_box']['center_z']})` | Box Size: `({dock_params['grid_box']['size_x']} × {dock_params['grid_box']['size_y']} × {dock_params['grid_box']['size_z']} Å)` | Engine: `{dock_params.get('engine_mode', 'vina').upper()}`")

                if run_dock_btn:
                    top_5_candidates = filtered_df.head(5).to_dict(orient='records')
                    with st.spinner("Executing molecular docking for Top 5 candidates..."):
                        dock_results = dock_top_candidates(
                            receptor_pdb=target_pdb_content,
                            candidates=top_5_candidates,
                            top_n=5,
                            grid_box=dock_params['grid_box'],
                            exhaustiveness=dock_params['exhaustiveness'],
                            n_poses=dock_params['n_poses'],
                            energy_range=dock_params['energy_range'],
                            seed=dock_params['seed'],
                            cpu=dock_params['cpu'],
                            sf_name=dock_params['sf_name'],
                            is_rna=False,
                            engine_mode=dock_params.get('engine_mode', 'vina'),
                            gpu_threads=dock_params.get('gpu_threads', 8000),
                            vina_gpu_bin=dock_params.get('vina_gpu_bin')
                        )
                        st.session_state['prozena_dock_results'] = dock_results
                        st.session_state['prozena_dock_params'] = dock_params
                        render_success_alert("AutoDock Vina Docking Complete for Top 5 Hits!")

            if 'prozena_dock_results' in st.session_state:
                dock_res_list = st.session_state['prozena_dock_results']
                saved_params = st.session_state.get('prozena_dock_params', dock_params or {})
                st.markdown("#### Docking Affinity Rankings")
                
                dock_summary = []
                for d in dock_res_list:
                    dock_summary.append({
                        'Compound': d.get('name'),
                        'Affinity (kcal/mol)': d.get('affinity'),
                        'Estimated Ki': d.get('ki_estimated'),
                        'Engine': d.get('engine_used', 'AutoDock Vina'),
                        'Status': d.get('status')
                    })
                st.dataframe(pd.DataFrame(dock_summary), use_container_width=True)

                st.markdown("#### 3D Complex & 2D Interaction Inspector")
                comp_names = [d.get('name') for d in dock_res_list if d.get('status') == 'Success']
                if comp_names:
                    selected_comp = st.selectbox("Select Pose to Inspect in 3D & 2D", comp_names)
                    selected_data = next((d for d in dock_res_list if d.get('name') == selected_comp), None)
                    if selected_data and selected_data.get('complex_pdb'):
                        col_3d, col_2d = st.columns([3, 2])
                        with col_3d:
                            st.caption(f"Interactive 3D View: {selected_comp} in Target Protein")
                            show_box = saved_params.get('show_box', True)
                            render_complex_3d(
                                selected_data['complex_pdb'],
                                height=480,
                                grid_box=saved_params.get('grid_box'),
                                show_box=show_box
                            )
                        with col_2d:
                            st.caption("2D Binding Residue Interactions & Topology")
                            inter_summary = analyze_interactions(selected_data['complex_pdb'])
                            render_interaction_summary_ui(inter_summary, complex_pdb_str=selected_data['complex_pdb'])
                            st.download_button(
                                "Download Complex PDB",
                                selected_data['complex_pdb'],
                                f"{selected_comp}_protein_complex.pdb",
                                "chemical/x-pdb",
                                use_container_width=True
                            )

    render_footer()
