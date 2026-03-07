import streamlit as st
import time
import pandas as pd
from rdkit import Chem
from .smartbind_wrapper import run_binding_site_prediction
from .visualization import highlight_rna_sequence, render_ligand_2d

def render_page():
    st.markdown("## 🎯 Binding Site Prediction")
    st.markdown("<p style='color: var(--text-muted);'>Predict per-residue binding probabilities for a specific RNA and ligand pair.</p>", unsafe_allow_html=True)
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        rna_seq = st.text_area("RNA Sequence", height=100)
    with col2:
        input_type = st.radio("Ligand Input Type", ["Paste SMILES", "Upload SDF"])
        ligand_smiles = ""
        
        if input_type == "Paste SMILES":
            ligand_smiles = st.text_input("Ligand SMILES", "CC(C)(C)c1ccc2occ(CC(=O)Nc3ccccc3F)c2c1")
        else:
            uploaded_file = st.file_uploader("Upload Ligand (.sdf)", type=['sdf'])
            if uploaded_file:
                supplier = Chem.ForwardSDMolSupplier(uploaded_file)
                for mol in supplier:
                    if mol is not None:
                        try:
                            ligand_smiles = Chem.MolToSmiles(mol)
                            break # Just take the first valid molecule for binding site
                        except:
                            pass
                st.info(f"Parsed SMILES: {ligand_smiles}")
                
        if ligand_smiles:
            st.markdown("**Ligand 2D Structure**")
            render_ligand_2d(ligand_smiles)
    st.markdown("</div>", unsafe_allow_html=True)

    if st.button("🔮 Predict Binding Sites", use_container_width=True):
        if not rna_seq or not ligand_smiles:
            st.error("Please provide both RNA Sequence and Ligand SMILES.")
            return
            
        with st.spinner("Analyzing sequence with RNA Foundation Model..."):
            time.sleep(1) # simulate loading
            
            probabilities = run_binding_site_prediction(rna_seq, ligand_smiles)
            
            st.success("Prediction Complete!")
            st.markdown("### 🧬 Binding Site Heatmap")
            st.markdown("Residues with high probabilities are highlighted. Hover over characters to see exact scores.")
            
            html_seq = highlight_rna_sequence(rna_seq, probabilities)
            st.markdown(f"<div class='rna-sequence-container'>{html_seq}</div>", unsafe_allow_html=True)
            
            st.markdown("---")
            st.markdown("### 📋 Tabular Results")
            
            # Create a dataframe: Position (Index), Nucleotide, Probability
            df_results = pd.DataFrame({
                "Position": range(1, len(rna_seq) + 1),
                "Nucleotide": list(rna_seq),
                "Probability": probabilities
            })
            
            st.dataframe(
                df_results.style.background_gradient(cmap='viridis', subset=['Probability']),
                use_container_width=True,
                height=300
            )

            # Export > 0.6 Probabilities
            st.markdown("---")
            st.markdown("### 📥 Export High-Probability Sites")
            st.markdown("Download a `.csv` containing only the residues with a binding probability > 0.6.")
            
            df_filtered = df_results[df_results['Probability'] > 0.6]
            
            if not df_filtered.empty:
                csv_data = df_filtered.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label=f"Download {len(df_filtered)} Sites (.csv)",
                    data=csv_data,
                    file_name="rnahunter_binding_sites_gt_06.csv",
                    mime="text/csv",
                    use_container_width=True
                )
            else:
                st.info("No binding sites found with probability > 0.6.")
