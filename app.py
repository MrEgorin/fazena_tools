import streamlit as st
from modules.prozena.protein_screening import ProteinScreening
from modules.rna_hunter.rna_binding_prediction import RNAHunter
from modules.refazenation.drug_repurposing import DrugRepurposing
from modules.docking.auto_docking import AutoDocking
from modules.admet_filter.filter_admet import ADMEFilter

def main():
    st.title("Fazena Tools - Modular Drug Discovery Platform")

    product = st.selectbox(
        "Select Product",
        ["RNA Hunter", "Refazenation", "Prozena"],
        index=0
    )

    if product == "RNA Hunter":
        page = st.radio(
            "Sections",
            ["Dashboard", "Virtual Screening", "Binding Site Search"],
            label_visibility="visible",
        )
    elif product == "Refazenation":
        page = st.radio(
            "Sections",
            ["Reverse Docking"],
            label_visibility="visible",
        )
    else:  # Prozena
        page = st.radio(
            "Sections",
            ["HTS Screening"],
            label_visibility="visible",
        )

    if product == "RNA Hunter":
        if page == "Dashboard":
            st.markdown(
                "<h1 class='hero-title'>RNA Hunter</h1>",
                unsafe_allow_html=True
            )
            st.markdown(
                "<p style='color: var(--fz-text-muted); font-size: 1.1rem; "
                "margin-bottom: 2.5rem; line-height: 1.7;'>"
                "Advanced Structure-Agnostic RNA Ligand Discovery</p>",
                unsafe_allow_html=True
            )

        elif page == "Virtual Screening":
            st.header("RNA Virtual Screening")
            rna_sequence = st.text_input("Enter RNA Sequence", value="AUCG")
            smiles_file = st.file_uploader("Upload SMILES File", type=["smi"])
            model_path = st.text_input("Model Path", value="path/to/model.pth")
            output_file = st.text_input("Output File", value="output.csv")

            if st.button("Run Virtual Screening"):
                rna_hunter = RNAHunter(model_path)
                results = rna_hunter.predict(rna_sequence)
                st.write(f"Prediction Results: {results}")

        elif page == "Binding Site Search":
            st.header("RNA Binding Site Search")
            # Add binding site search logic here

    elif product == "Refazenation":
        if page == "Reverse Docking":
            st.header("Drug Repurposing and Fast Target Prediction")
            smiles = st.text_input("Enter SMILES", value="CCO")
            if st.button("Run Drug Repurposing"):
                drug_repurposer = DrugRepurposer()
                results = drug_repurposer.predict(smiles)
                st.write(f"Prediction Results: {results}")

    elif product == "Prozena":
        if page == "HTS Screening":
            st.header("Protein Targeted Screening & Docking")
            protein_seq = st.text_input("Enter Protein Sequence", value="MKTVRQERLKSIVRILGARFTNVTGIQKLVCAIVI")
            ligands_file = st.file_uploader("Upload Ligands File", type=["smi"])
            model_path = st.text_input("Model Path", value="path/to/model.pth")
            output_file = st.text_input("Output File", value="output.csv")

            if st.button("Run Protein Screening"):
                protein_screener = ProteinScreening(model_path)
                results = protein_screener.predict(protein_seq)
                st.write(f"Prediction Results: {results}")

if __name__ == "__main__":
    main()
