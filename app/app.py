import streamlit as st
from modules.admet_filter.filter_admet import ADMEFilter
from modules.prozena.protein_screening import ProteinScreening
from modules.rna_hunter.rna_binding_prediction import RNAHunter
from modules.refazenation.drug_repurposing import DrugRepurposing
from modules.docking.auto_docking import AutoDocking
from modules.visualization.interaction_plots import InteractionPlots
import os

# Inject custom CSS for styling
def inject_custom_css():
    with open("app/style.css") as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

inject_custom_css()

# Sidebar Navigation / Product Selector
st.sidebar.title("Fazena Tools")
product_selector = st.sidebar.selectbox(
    "Select a Tool",
    ["Prozena", "RNA-Hunter", "Refazenation"]
)

# Initialize modules
adme_filter = ADMEFilter()
protein_screening = ProteinScreening(model_path="path/to/protein_model.pkl")
rna_hunter = RNAHunter(model_path="path/to/rna_model.pkl")
drug_repurposing = DrugRepurposing()
auto_docking = AutoDocking(vina_path="path/to/vina")
interaction_plots = InteractionPlots()

# Input & Parameter Controls
st.title(f"{product_selector} Tool")

if product_selector == "Prozena":
    st.subheader("Protein Screening & GPU Docking")
    compound_library_file = st.file_uploader("Upload Compound Library (SMILES/SDF)", type=["smi", "sdf"])
    protein_target_file = st.file_uploader("Upload Protein Target (PDB/PDBQT)", type=["pdb", "pdbqt"])
    admet_filter_params = {
        'max_mw': st.slider('Max Molecular Weight', min_value=100, max_value=2000, value=500),
        'max_logp': st.slider('Max LogP', min_value=-10, max_value=10, value=5),
        'max_hbd': st.slider('Max H-bond Donors', min_value=0, max_value=20, value=5),
        'max_hba': st.slider('Max H-bond Acceptors', min_value=0, max_value=30, value=10),
        'max_tpsa': st.slider('Max TPSA', min_value=0, max_value=400, value=100)
    }
    grid_box_params = {
        'center_x': st.number_input('Grid Box Center X', value=0.0),
        'center_y': st.number_input('Grid Box Center Y', value=0.0),
        'center_z': st.number_input('Grid Box Center Z', value=0.0),
        'size_x': st.number_input('Grid Box Size X', value=10.0),
        'size_y': st.number_input('Grid Box Size Y', value=10.0),
        'size_z': st.number_input('Grid Box Size Z', value=10.0),
        'exhaustiveness': st.slider('Exhaustiveness', min_value=1, max_value=32, value=8)
    }
    output_dir = st.text_input('Output Directory', value="output/prozena")

elif product_selector == "RNA-Hunter":
    st.subheader("RNA Screening & Docking")
    compound_library_file = st.file_uploader("Upload Compound Library (SMILES/SDF)", type=["smi", "sdf"])
    rna_target_file = st.file_uploader("Upload RNA Target (PDB/PDBQT)", type=["pdb", "pdbqt"])
    admet_filter_params = {
        'max_mw': st.slider('Max Molecular Weight', min_value=100, max_value=2000, value=500),
        'max_logp': st.slider('Max LogP', min_value=-10, max_value=10, value=5),
        'max_hbd': st.slider('Max H-bond Donors', min_value=0, max_value=20, value=5),
        'max_hba': st.slider('Max H-bond Acceptors', min_value=0, max_value=30, value=10),
        'max_tpsa': st.slider('Max TPSA', min_value=0, max_value=400, value=100)
    }
    grid_box_params = {
        'center_x': st.number_input('Grid Box Center X', value=0.0),
        'center_y': st.number_input('Grid Box Center Y', value=0.0),
        'center_z': st.number_input('Grid Box Center Z', value=0.0),
        'size_x': st.number_input('Grid Box Size X', value=10.0),
        'size_y': st.number_input('Grid Box Size Y', value=10.0),
        'size_z': st.number_input('Grid Box Size Z', value=10.0),
        'exhaustiveness': st.slider('Exhaustiveness', min_value=1, max_value=32, value=8)
    }
    output_dir = st.text_input('Output Directory', value="output/rna_hunter")

elif product_selector == "Refazenation":
    st.subheader("Drug Repurposing & Target Prediction")
    compound_library_file = st.file_uploader("Upload Compound Library (SMILES/SDF)", type=["smi", "sdf"])
    drug_library_file = st.file_uploader("Upload Drug Library (CSV)", type=["csv"])
    output_dir = st.text_input('Output Directory', value="output/refazenation")

# Execution & Reporting
if st.button("Run Pre-Filtering (ADMET)"):
    if compound_library_file is None:
        st.error("Please upload a compound library.")
    else:
        # Placeholder for pre-filtering logic
        st.success("Pre-filtering completed.")

if st.button("Run Screening"):
    if compound_library_file is None or protein_target_file is None:
        st.error("Please upload both compound and target files.")
    else:
        # Placeholder for screening logic
        st.success("Screening completed.")

if st.button("Run Automated Top-5 Docking"):
    if compound_library_file is None or protein_target_file is None:
        st.error("Please upload both compound and target files.")
    else:
        # Placeholder for docking logic
        st.success("Docking completed.")

if st.button("Generate Visualizations"):
    if compound_library_file is None or protein_target_file is None:
        st.error("Please upload both compound and target files.")
    else:
        # Placeholder for visualization logic
        st.success("Visualizations generated.")

# Display results in structured tables, interactive charts, 2D interaction plots
st.subheader("Results")
# Placeholder for displaying results
