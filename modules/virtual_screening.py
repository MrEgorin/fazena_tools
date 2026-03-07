import streamlit as st
import pandas as pd
import time
import os
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Draw
from io import BytesIO
import base64
from .smartbind_wrapper import run_virtual_screening

# Basic Toxicophores / Reactive groups SMARTS patterns
TOXIC_SMARTS = {
    'Epoxide': '[OX3]1[CX4][CX4]1',
    'Quinone': 'O=C1[#6]~[#6]C(=O)[#6]~[#6]1',
    'Nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',
    'Michael_Acceptor': 'C=C-C(=O)',
    'Aliphatic_Halide': '[CX4][F,Cl,Br,I]',
}

def check_toxicity(mol):
    """Returns True if any toxicophores are found, False otherwise."""
    if not mol: return False
    for name, smarts in TOXIC_SMARTS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True, name
    return False, "None"

def get_admet_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    has_tox, tox_alert = check_toxicity(mol)
    
    return {
        'MW': round(mw, 2),
        'LogP': round(logp, 2),
        'HBD': hbd,
        'HBA': hba,
        'TPSA': round(tpsa, 2),
        'Tox_Alert': tox_alert
    }

def filter_top_candidates(df, max_mw, max_logp, max_hbd, max_hba, max_tpsa, filter_tox, target_count=100, min_score=0.8):
    """
    Iterates through the DataFrame (sorted by Binding_Score desc), 
    calculates ADMET, and collects up to target_count molecules 
    that pass the filters and have score >= min_score.
    """
    df = df.sort_values(by='Binding_Score', ascending=False).reset_index(drop=True)
    
    passed_records = []
    
    for _, row in df.iterrows():
        # Stop condition 1: reached target count
        if len(passed_records) >= target_count:
            break
            
        # Stop condition 2: score dropped below min_score
        if row['Binding_Score'] < min_score:
            break
            
        smiles = row['SMILES']
        props = get_admet_properties(smiles)
        
        if not props:
            continue
            
        # Check Filters
        passed = True
        if props['MW'] > max_mw: passed = False
        if props['LogP'] > max_logp: passed = False
        if props['HBD'] > max_hbd: passed = False
        if props['HBA'] > max_hba: passed = False
        if props['TPSA'] > max_tpsa: passed = False
        if filter_tox and props['Tox_Alert'] != 'None': passed = False
        
        if passed:
            record = row.to_dict()
            record.update(props)
            passed_records.append(record)
            
    return pd.DataFrame(passed_records)

def get_base64_img(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        img = Draw.MolToImage(mol, size=(250, 250))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    except:
        return ""

def render_page():
    st.markdown("## 🔍 Virtual Target Screening")
    st.markdown("<p style='color: var(--text-muted);'>Screen a library of small molecules against an RNA target to predict binding affinities.</p>", unsafe_allow_html=True)
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.markdown("### 1️⃣ RNA Target")
        rna_seq = st.text_area("Input RNA Sequence", height=150, placeholder="e.g., GGCGCUUGCUGGGUAAAA...")
        
    with col2:
        st.markdown("### 2️⃣ Ligand Library")
        library_source = st.radio("Select Library Source", ["Paste SMILES", "Upload File (.smi, .csv, .sdf)"])
        
        ligands_input = []
        if library_source == "Paste SMILES":
            smiles_text = st.text_area("Input SMILES (one per line)", height=150, placeholder="CCO\nC1CCCCC1")
            if smiles_text:
                ligands_input = [s.strip() for s in smiles_text.split() if s.strip()]
        elif library_source == "Upload File (.smi, .csv, .sdf)":
            uploaded_file = st.file_uploader("Upload Ligand Library", type=['smi', 'csv', 'txt', 'sdf'])
            if uploaded_file:
                filename = uploaded_file.name.lower()
                if filename.endswith('.sdf'):
                    # Use RDKit to parse SDF from binary buffer
                    supplier = Chem.ForwardSDMolSupplier(uploaded_file)
                    for mol in supplier:
                        if mol is not None:
                            try:
                                smi = Chem.MolToSmiles(mol)
                                ligands_input.append(smi)
                            except:
                                pass
                else:
                    # Basic reading for .smi / .csv
                    content = uploaded_file.getvalue().decode('utf-8').splitlines()
                    ligands_input = [line.strip().split()[0] for line in content if line.strip()]
                    
                st.info(f"Loaded {len(ligands_input)} ligands.")
        else:
            st.warning("Please choose a library source.")
    st.markdown("</div>", unsafe_allow_html=True)
    st.markdown("---")
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    with st.expander("⚙️ Advanced ADMET & Toxicity Filters", expanded=False):
        st.markdown("Molecules passing these filters with a Binding Score $\\ge$ 0.8 will be prioritized in the Top 100.")
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            max_mw = st.number_input("Max Molecular Weight", value=500.0)
            max_logp = st.number_input("Max LogP", value=5.0)
        with col_f2:
            max_hbd = st.number_input("Max H-Bond Donors", value=5)
            max_hba = st.number_input("Max H-Bond Acceptors", value=10)
        with col_f3:
            max_tpsa = st.number_input("Max TPSA", value=140.0)
            filter_tox = st.checkbox("Filter known Toxicophores", value=True, help="Removes epoxides, quinones, pan-assay interference (PAINS) proxies, and other reactive groups.")
    st.markdown("</div>", unsafe_allow_html=True)
            
    st.markdown("---")
    
    if st.button("🚀 Run RNAHunter Virtual Screening", use_container_width=True):
        if not rna_seq:
            st.error("Please provide an RNA sequence.")
            return
        if not ligands_input:
            st.error("Please provide at least one ligand.")
            return
            
        with st.spinner("Initializing models and running contrastive screening..."):
            # Progress bar simulation for UX
            progress_bar = st.progress(0)
            for i in range(100):
                time.sleep(0.01)
                progress_bar.progress(i + 1)
                
            # Run wrapper
            results_df = run_virtual_screening(rna_seq, ligands_input)
            
            progress_bar.empty()
            st.success("Screening Complete!")
            
            # Display results
            st.markdown("### 📊 Top Binding Candidates")
            
            # Show a metric
            col_a, col_b, col_c = st.columns(3)
            col_a.metric("Total Ligands Screened", len(results_df))
            top_score = results_df['Binding_Score'].max()
            col_b.metric("Highest Affinity Score", f"{top_score:.2f}")
            col_c.metric("Computation Time", "~1.5s (Demo)")
            
            # Generate 2D images for Top 100 candidates
            with st.spinner("Applying ADMET Filters and Generating Visualizations..."):
                display_df = filter_top_candidates(
                    results_df, 
                    max_mw=max_mw, 
                    max_logp=max_logp, 
                    max_hbd=max_hbd, 
                    max_hba=max_hba, 
                    max_tpsa=max_tpsa, 
                    filter_tox=filter_tox, 
                    target_count=100, 
                    min_score=0.8
                )
                
                if len(display_df) > 0:
                    display_df['Structure'] = display_df['SMILES'].apply(get_base64_img)
                else:
                    st.warning("No molecules passed the filters with a score >= 0.8.")
            
            if len(display_df) > 0:
                # Rearrange columns
                display_df = display_df[['Structure', 'SMILES', 'Binding_Score', 'MW', 'LogP', 'HBD', 'HBA', 'TPSA']]
            
            # To bypass Streamlit's restrictive row height on st.dataframe for images, 
            # we will render a stylized HTML table for the Top 100.
            
            # Create interactive HTML table
            html_table = "<div style='max-height: 600px; overflow-y: auto; border-radius: 12px; border: 1px solid var(--glass-border);'>"
            html_table += "<table style='width: 100%; border-collapse: collapse; text-align: left; color: var(--text-main);'>"
            html_table += "<thead style='background-color: var(--secondary-bg); position: sticky; top: 0; z-index: 1;'>"
            html_table += "<tr><th style='padding: 12px; border-bottom: 1px solid var(--glass-border);'>Rank</th>"
            html_table += "<th style='padding: 12px; border-bottom: 1px solid var(--glass-border);'>2D Structure</th>"
            html_table += "<th style='padding: 12px; border-bottom: 1px solid var(--glass-border);'>SMILES</th>"
            html_table += "<th style='padding: 12px; border-bottom: 1px solid var(--glass-border);'>Binding Score</th>"
            html_table += "<th style='padding: 12px; border-bottom: 1px solid var(--glass-border);'>ADMET Properties</th></tr></thead><tbody>"
            
            for idx, row in display_df.iterrows():
                bg_color = "var(--primary-bg)" if idx % 2 == 0 else "var(--secondary-bg)"
                # Add highlighting for top 3
                rank = idx + 1
                rank_str = f"🥇 {rank}" if rank == 1 else f"🥈 {rank}" if rank == 2 else f"🥉 {rank}" if rank == 3 else f"#{rank}"
                
                img_tag = f"<img src='{row['Structure']}' width='200' style='border-radius: 8px;'/>" if row['Structure'] else "N/A"
                
                html_table += f"<tr style='background-color: {bg_color}; border-bottom: 1px solid var(--glass-border);'>"
                html_table += f"<td style='padding: 12px; font-weight: bold;'>{rank_str}</td>"
                html_table += f"<td style='padding: 12px;'>{img_tag}</td>"
                html_table += f"<td style='padding: 12px; font-family: monospace; word-break: break-all; max-width: 200px;'>{row['SMILES']}</td>"
                html_table += f"<td style='padding: 12px; font-weight: bold; color: var(--accent-orange);'>{row['Binding_Score']:.4f}</td>"
                
                admet_str = f"MW: {row.get('MW', 'N/A')}<br>LogP: {row.get('LogP', 'N/A')}<br>HBD/HBA: {row.get('HBD', 'N/A')}/{row.get('HBA', 'N/A')}<br>TPSA: {row.get('TPSA', 'N/A')}"
                html_table += f"<td style='padding: 12px; font-size: 0.9em; color: var(--text-muted);'>{admet_str}</td>"
                html_table += "</tr>"
                
            html_table += "</tbody></table></div>"
            
            if len(display_df) > 0:
                st.markdown(f"**Found {len(display_df)} candidates** meeting ADMET filters and binding confidence $\ge$ 0.8.")
            st.markdown(html_table, unsafe_allow_html=True)

            st.markdown("---")
            st.markdown("### 📈 Screening Statistics & Distribution")
            
            # Statistics
            stats_a, stats_b, stats_c, stats_d = st.columns(4)
            stats_a.metric("Mean Score", f"{results_df['Binding_Score'].mean():.4f}")
            stats_b.metric("Median Score", f"{results_df['Binding_Score'].median():.4f}")
            stats_c.metric("Min Score", f"{results_df['Binding_Score'].min():.4f}")
            stats_d.metric("Max Score", f"{results_df['Binding_Score'].max():.4f}")
            
            # Histogram
            import numpy as np
            counts, bins = np.histogram(results_df['Binding_Score'], bins=20)
            hist_df = pd.DataFrame({'Frequency (# of Ligands)': counts}, index=np.round(bins[:-1], 3))
            st.bar_chart(hist_df, use_container_width=True)
            
            st.markdown("---")
            st.markdown("### 📥 Export Results")
            col_export1, col_export2 = st.columns(2)
            
            csv_all = results_df.to_csv(index=False).encode('utf-8')
            # Drop 'Structure' for CSV exports
            csv_export_df = display_df.drop(columns=['Structure'], errors='ignore')
            csv_top100 = csv_export_df.to_csv(index=False).encode('utf-8')
            
            with col_export1:
                st.download_button(
                    label="Download ALL Results (.csv)",
                    data=csv_all,
                    file_name='rnahunter_vs_all_results.csv',
                    mime='text/csv',
                    use_container_width=True
                )
            with col_export2:
                st.download_button(
                    label=f"Download Filtered TOP {len(display_df)} (.csv)",
                    data=csv_top100,
                    file_name='rnahunter_vs_filtered_results.csv',
                    mime='text/csv',
                    use_container_width=True
                )
