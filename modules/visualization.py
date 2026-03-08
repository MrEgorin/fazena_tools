import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def render_ligand_2d(smiles: str):
    """
    Renders a 2D image of a ligand given its SMILES string within a Streamlit app.
    Styled to match the Fazena dark design language.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(300, 300))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            
            st.markdown(
                f"<div style='text-align: center; background: var(--fz-bg-primary); "
                f"border: 1px solid var(--fz-border); border-radius: var(--fz-radius-card); "
                f"padding: 1rem; display: inline-block;'>"
                f"<img src='data:image/png;base64,{img_str}' "
                f"style='border-radius: 8px; max-width: 280px;'/>"
                f"</div>",
                unsafe_allow_html=True
            )
        else:
            st.error(f"Invalid SMILES string: {smiles}")
    except Exception as e:
        st.error(f"Failed to render ligand: {e}")

def highlight_rna_sequence(sequence: str, probabilities: list) -> str:
    """
    Takes an RNA sequence and a list of probabilities (0.0 to 1.0) of the same length.
    Returns an HTML string with each nucleotide highlighted using the Fazena orange accent.
    High probability → strong orange. Low probability → subtle.
    """
    if len(sequence) != len(probabilities):
        return (
            f"Error: Sequence length ({len(sequence)}) does not match "
            f"probabilities length ({len(probabilities)})."
        )
    
    html_output = ""
    for char, prob in zip(sequence, probabilities):
        # Use orange accent color with varying alpha based on probability
        alpha = max(0.05, prob)
        
        tooltip = f"{char}: {prob:.3f}"
        
        # Bolder text for high-probability residues
        font_weight = "bold" if prob > 0.5 else "normal"
        # Brighter text for higher scores
        color = "#ffffff" if prob > 0.3 else "rgba(255, 255, 255, 0.5)"
        
        # Use orange-based background matching the accent
        bg_color = f"rgba(218, 117, 45, {alpha})"
        
        span = (
            f"<span title='{tooltip}' style='"
            f"background-color: {bg_color}; color: {color}; "
            f"font-weight: {font_weight}; padding: 2px 1px; "
            f"border-radius: 2px; transition: background-color 0.2s;'>"
            f"{char}</span>"
        )
        html_output += span
        
    return html_output
