import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import base64
from io import BytesIO

def render_ligand_2d(smiles: str):
    """
    Renders a 2D image of a ligand given its SMILES string within a Streamlit app.
    Provides fallback to text if RDKit fails or SMILES is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate image in memory
            img = Draw.MolToImage(mol, size=(300, 300))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            
            # Display image centered
            st.markdown(
                f"<div style='text-align: center;'><img src='data:image/png;base64,{img_str}' style='border-radius: 8px;'/></div>",
                unsafe_allow_html=True
            )
        else:
            st.error(f"Invalid SMILES string: {smiles}")
    except Exception as e:
        st.error(f"Failed to render ligand: {e}")

def highlight_rna_sequence(sequence: str, probabilities: list) -> str:
    """
    Takes an RNA sequence and a list of probabilities (0.0 to 1.0) of the same length.
    Returns an HTML string where each nucleotide is wrapped in a span with a background
    color indicating its binding probability. Let's use a heatmap from transparent to deep accent (cyan).
    """
    if len(sequence) != len(probabilities):
        return f"Error: Sequence length ({len(sequence)}) does not match probabilities length ({len(probabilities)})."
    
    html_output = ""
    for char, prob in zip(sequence, probabilities):
        # Calculate RGBA color based on probability. Let's use a green/cyan hue: rgba(0, 242, 254, alpha)
        # Ensure minimum visibility
        alpha = max(0.05, prob)
        
        # Tooltip for exact score
        tooltip = f"{char}: {prob:.3f}"
        
        # Style
        # Bolder text for scores > 0.5
        font_weight = "bold" if prob > 0.5 else "normal"
        color = "#ffffff" if prob > 0.3 else "#a0aabf"
        
        span = f"<span title='{tooltip}' style='background-color: rgba(0, 242, 254, {alpha}); color: {color}; font-weight: {font_weight}; padding: 2px 1px; border-radius: 2px; transition: background-color 0.2s;'>{char}</span>"
        html_output += span
        
    return html_output
