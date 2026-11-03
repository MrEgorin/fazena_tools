"""
FAZENA Tools — RNA Hunter Virtual Screening
Foundry-scale virtual screening against RNA architectures with automated AutoDock Vina docking.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

import streamlit as st
import pandas as pd
import time
import os
from rdkit import Chem
from modules.smartbind_wrapper import run_virtual_screening
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

# SVG Icons
ICON_HELIX = ('<svg width="22" height="22" viewBox="0 0 512 512" fill="none" xmlns="http://www.w3.org/2000/svg">'
    '<path d="M150.815 0L160.901 0C163.181 0.913013 165.648 1.1956 167.927 2.05723C176.627 5.3469 183.973 11.9088 187.816 20.4291C189.266 23.6033 190.218 26.9817 190.64 30.4456C191.333 35.9937 190.903 42.4391 191.01 48.0944C191.254 60.969 190.609 74.692 191.188 87.4759C192.454 112.813 199.326 137.553 211.308 159.913C224.218 184.108 242.749 204.85 265.343 220.392C275.329 227.236 284.55 232.164 294.558 239.707C351.742 282.685 386.923 348.809 390.612 420.248C391.293 432.191 391.007 445.089 391.016 457.143C391.021 463.433 391.35 474.232 390.785 480.139C390.509 483.358 389.782 486.522 388.625 489.538C385.312 498.136 378.697 505.135 370.265 508.951C368.183 509.893 362.052 511.475 361.189 512L351.099 512C349.619 511.283 347.389 510.957 345.78 510.476C335.105 507.335 325.843 498.304 322.616 487.606C320.144 479.409 321.155 469.254 321.01 460.653C320.831 450.072 321.177 439.414 320.957 428.82C320.583 414.857 318.523 400.992 314.82 387.524C304.443 349.851 281.548 316.833 249.903 293.907C243.32 289.15 236.658 285.361 229.889 280.975C222.927 276.462 216.212 271.577 209.774 266.343C177.896 240.396 153.227 206.688 138.137 168.455C128.405 143.821 122.728 117.774 121.326 91.3243C120.469 74.5769 121.473 57.8148 120.991 40.922C120.502 23.7554 125.88 9.87253 142.714 2.58691C145.25 1.48965 148.096 1.08011 150.622 0.0787148L150.815 0Z" fill="var(--fz-accent)"/>'
    '</svg>')
ICON_HTVS = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M3 3h18v18H3zM3 9h18M3 15h18M9 3v18M15 3v18"/><path d="M9 9l6 6M15 9l-6 6" stroke-width="2" stroke="var(--fz-accent)"/></svg>'
ICON_STATS = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M3 20h18M3 16h4v4H3zM9 12h4v8H9zM15 8h4v12h-4zM3 12l6-4 6 4 6-8"/></svg>'
ICON_EXPORT = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>'
ICON_DOCK = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 2v20M2 12h20"/><circle cx="12" cy="12" r="4"/></svg>'

def render_page():
    render_header(
        "RNA Hunter | Virtual Screening",
        "Screen complex chemical libraries against RNA targets using physics-informed bio-models with automated AutoDock Vina docking."
    )

    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 1])

    with col1:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HELIX}</span> Target RNA Architecture</h3>", unsafe_allow_html=True)
        rna_input = st.text_area(
            "Input RNA Sequence (Raw or FASTA)",
            height=140,
            value="GGCCUUAGCGUAGGCUAGGCCAAAA",
            placeholder="e.g. GGCGCUUGCUGGGUAAAA..."
        )
        rna_seq, rna_name = parse_sequence_input(rna_input, seq_type='rna')
        if rna_seq:
            st.caption(f"Parsed RNA Sequence: **{len(rna_seq)} nt** ({rna_name})")

        rna_pdb_file = st.file_uploader("Optional 3D RNA Structure (.pdb) for Automated Docking (default: built-in rna_model.pdb)", type=['pdb', 'pdbqt'])

        default_rna_model_path = os.path.join(os.path.dirname(__file__), "rna_model.pdb")
        target_preview_content = None
        if rna_pdb_file:
            target_preview_content = rna_pdb_file.getvalue().decode('utf-8', errors='ignore')
            st.caption(f"Loaded uploaded RNA: **{rna_pdb_file.name}**")
        elif os.path.exists(default_rna_model_path):
            with open(default_rna_model_path, 'r') as f:
                target_preview_content = f.read()

        if target_preview_content:
            with st.expander("Target RNA 3D Architecture Preview", expanded=False):
                st.caption("3D macromolecular view of target RNA architecture:")
                render_receptor_3d(target_preview_content, height=320)

    with col2:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HTVS}</span> Chemical Library Input</h3>", unsafe_allow_html=True)
        lib_source = st.radio("Library Source", ["Paste SMILES (Multi-line / Library)", "Upload Library File (.sdf, .csv, .tsv, .smi)"], horizontal=True)
        
        parsed_ligands = []
        if "Paste SMILES" in lib_source:
            default_smiles = (
                "COc1ccc2occ(CC(=O)Nc3ccccc3F)c2c1\n"
                "CCN(CC)CCNC(=O)c1c(C)[nH]c(C=C2C(=O)Nc3ccc(F)cc32)c1C\n"
                "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1\n"
                "CC(=O)Oc1ccccc1C(=O)O\n"
                "CN1CCN(CC1)C2=NC3=CC=CC=C3NC4=C2C=C(C=C4)Cl"
            )
            smiles_text = st.text_area("SMILES (one per line, with optional name)", height=140, value=default_smiles)
            parsed_ligands = parse_ligand_text(smiles_text)
        else:
            uploaded_file = st.file_uploader("Upload Chemical Library", type=['sdf', 'mol', 'smi', 'smiles', 'csv', 'tsv', 'txt'])
            if uploaded_file:
                parsed_ligands = parse_ligand_file(uploaded_file.getvalue(), uploaded_file.name)

        if parsed_ligands:
            st.info(f"Loaded **{len(parsed_ligands)}** valid small molecules in chemical library.")

    st.markdown("</div>", unsafe_allow_html=True)

    render_section_divider()

    # ADMET & Safety Filters
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    with st.expander("Advanced ADMET & Safety Filters", expanded=True):
        st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.9rem;'>Compounds violating Lipinski Rule of 5 or containing reactive toxicophores will be flagged or excluded.</p>", unsafe_allow_html=True)
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            max_mw = st.number_input("Max Molecular Weight (MW)", value=550.0, step=10.0)
            max_logp = st.number_input("Max LogP", value=5.0, step=0.5)
        with col_f2:
            max_hbd = st.number_input("Max H-Bond Donors (HBD)", value=5, step=1)
            max_hba = st.number_input("Max H-Bond Acceptors (HBA)", value=10, step=1)
        with col_f3:
            max_tpsa = st.number_input("Max TPSA (Å²)", value=140.0, step=10.0)
            filter_tox = st.checkbox("Filter Toxicophores & Reactive Groups", value=True)
            filter_pains = st.checkbox("Apply Official RDKit PAINS Filters", value=True)
    st.markdown("</div>", unsafe_allow_html=True)

    render_section_divider()

    if st.button("Execute RNA Virtual Screening", use_container_width=True):
        if not rna_seq:
            st.error("Please provide a target RNA sequence."); return
        if not parsed_ligands:
            st.error("Please provide a chemical library."); return

        with st.spinner("Executing Sequence-Structure AI Inference Engine..."):
            progress_bar = st.progress(0)
            for i in range(100):
                time.sleep(0.005)
                progress_bar.progress(i + 1)

            results_df = run_virtual_screening(rna_seq, parsed_ligands)
            progress_bar.empty()
            st.session_state['rna_vs_results'] = results_df
            st.session_state['rna_vs_seq'] = rna_seq
            render_success_alert(f"Screening Complete! Evaluated {len(results_df)} compounds.")

    # Results Section
    if 'rna_vs_results' in st.session_state:
        results_df = st.session_state['rna_vs_results']

        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HTVS}</span> Screening Hit Table</h3>", unsafe_allow_html=True)
        
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

        col_m1, col_m2, col_m3 = st.columns(3)
        col_m1.metric("Total Screened", len(results_df))
        col_m2.metric("Passed ADMET", len(filtered_df))
        top_score = filtered_df['Binding_Score'].max() if not filtered_df.empty else 0.0
        col_m2_val = f"{top_score:.3f}" if top_score else "N/A"
        col_m3.metric("Top Binding Score", col_m2_val)

        if not filtered_df.empty:
            filtered_df['Structure'] = filtered_df['SMILES'].apply(get_base64_mol_image)
            
            html_table = "<div class='results-table-wrap'><table class='results-table'><thead><tr><th>Rank</th><th>Structure</th><th>Name / SMILES</th><th>Binding Score</th><th>ADMET Safety</th></tr></thead><tbody>"
            for idx, row in filtered_df.head(25).iterrows():
                bg = "rgba(255,255,255,0.03)" if idx % 2 == 0 else "transparent"
                img_tag = f"<img src='{row['Structure']}' width='110'/>" if row['Structure'] else "N/A"
                admet_str = f"<b>MW:</b> {row['MW']} | <b>LogP:</b> {row['LogP']}<br><b>HBD/HBA:</b> {row['HBD']}/{row['HBA']} | <b>TPSA:</b> {row['TPSA']}"
                html_table += f"<tr style='background-color:{bg}'><td><span class='rank-badge'>{idx+1}</span></td><td>{img_tag}</td><td><b>{row['name']}</b><br><span class='smiles-cell'>{row['SMILES']}</span></td><td class='score-cell' style='font-size:1.1rem; color:var(--fz-accent);'>{row['Binding_Score']:.4f}</td><td style='font-size:0.85rem;'>{admet_str}</td></tr>"
            html_table += "</tbody></table></div>"
            st.markdown(html_table, unsafe_allow_html=True)

            # Export
            render_section_divider()
            col_ex1, col_ex2 = st.columns(2)
            with col_ex1:
                st.download_button("Export All Results (.csv)", results_df.to_csv(index=False).encode('utf-8'), "rna_screen_all.csv", "text/csv", use_container_width=True)
            with col_ex2:
                st.download_button("Export Passed Filter Hits (.csv)", filtered_df.drop(columns=['Structure'], errors='ignore').to_csv(index=False).encode('utf-8'), "rna_screen_passed.csv", "text/csv", use_container_width=True)

            # ═════════════════════════════════════════════════════════════
            # AUTOMATED TOP-5 DOCKING SECTION WITH FULL PARAMETER TUNING
            # ═════════════════════════════════════════════════════════════
            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_DOCK}</span> Automated Top-5 AutoDock Vina RNA Docking</h3>", unsafe_allow_html=True)
            st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.92rem;'>Perform full 3D physics-based docking with AutoDock Vina for the top 5 screening candidates with flexible grid box and sampling parameters.</p>", unsafe_allow_html=True)

            default_rna_model_path = os.path.join(os.path.dirname(__file__), "rna_model.pdb")
            if rna_pdb_file:
                target_pdb_content = rna_pdb_file.getvalue().decode('utf-8', errors='ignore')
                st.info(f"Using uploaded 3D RNA Structure: **{rna_pdb_file.name}**")
            elif os.path.exists(default_rna_model_path):
                with open(default_rna_model_path, 'r') as f:
                    target_pdb_content = f.read()
                st.info("Using reference 3D RNA Target (`rna_model.pdb`) for automated docking. (You can also upload a custom .pdb above).")
            else:
                target_pdb_content = None
                st.warning("Please upload a 3D RNA Structure (.pdb) above to run full AutoDock Vina docking.")

            dock_params = None
            if target_pdb_content:
                default_box = calculate_grid_box(target_pdb_content)
                dock_params = render_docking_parameter_controls(default_box, is_rna=True, key_prefix="rna_dock", allow_gpu=True)

                col_btn, col_info = st.columns([1, 2])
                with col_btn:
                    run_dock_btn = st.button("Run AutoDock Vina on Top-5 Compounds", use_container_width=True)
                with col_info:
                    st.caption(f"Search Box Center: `({dock_params['grid_box']['center_x']}, {dock_params['grid_box']['center_y']}, {dock_params['grid_box']['center_z']})` | Box Size: `({dock_params['grid_box']['size_x']} × {dock_params['grid_box']['size_y']} × {dock_params['grid_box']['size_z']} Å)` | Engine: `{dock_params.get('engine_mode', 'vina').upper()}`")

                if run_dock_btn:
                    top_5_candidates = filtered_df.head(5).to_dict(orient='records')
                    with st.spinner("Executing AutoDock Vina molecular docking for Top 5 candidates..."):
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
                            is_rna=True,
                            engine_mode=dock_params.get('engine_mode', 'vina'),
                            gpu_threads=dock_params.get('gpu_threads', 8000),
                            vina_gpu_bin=dock_params.get('vina_gpu_bin')
                        )
                        st.session_state['rna_dock_results'] = dock_results
                        st.session_state['rna_dock_params'] = dock_params
                        render_success_alert("AutoDock Vina Docking Complete for Top 5 RNA Hits!")

            if 'rna_dock_results' in st.session_state:
                dock_res_list = st.session_state['rna_dock_results']
                saved_params = st.session_state.get('rna_dock_params', dock_params or {})
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

                # Interactive 3D Viewer & 2D Contacts for selected pose
                st.markdown("#### 3D Complex & 2D Interaction Inspector")
                comp_names = [d.get('name') for d in dock_res_list if d.get('status') == 'Success']
                if comp_names:
                    selected_comp = st.selectbox("Select Docked Pose to Inspect in 3D & 2D", comp_names)
                    selected_data = next((d for d in dock_res_list if d.get('name') == selected_comp), None)
                    if selected_data and selected_data.get('complex_pdb'):
                        col_3d, col_2d = st.columns([3, 2])
                        with col_3d:
                            st.caption(f"Interactive 3D View: {selected_comp} in Target RNA Structure")
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
                                f"{selected_comp}_rna_complex.pdb",
                                "chemical/x-pdb",
                                use_container_width=True
                            )

    render_footer()
