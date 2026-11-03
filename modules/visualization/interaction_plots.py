"""
FAZENA Tools — Molecular & Interaction Visualization Utilities
Provides 2D molecule rendering, RNA sequence heatmaps, 3D py3Dmol complex viewer,
2D protein-ligand interaction diagrams, and InteractionPlots class.
Author: Yehor Franchuk, MrEgorin, FAZENA
"""

import os
import base64
from io import BytesIO
from typing import List, Dict, Any, Optional

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import py3Dmol

try:
    from modules.docking_engine import analyze_interactions
except ImportError:
    analyze_interactions = None


def render_ligand_2d(smiles: str, width: int = 260, height: int = 260):
    """
    Renders a 2D image of a ligand given its SMILES string within a Streamlit app.
    Styled to match the Fazena dark design language.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(width, height))
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            
            st.markdown(
                f"<div style='text-align: center; background: var(--fz-bg-primary, #1D1D1D); "
                f"border: 1px solid var(--fz-border, #333); border-radius: var(--fz-radius-card, 12px); "
                f"padding: 1rem; display: inline-block;'>"
                f"<img src='data:image/png;base64,{img_str}' "
                f"style='border-radius: 8px; max-width: {width}px;'/>"
                f"</div>",
                unsafe_allow_html=True
            )
        else:
            st.error(f"Invalid SMILES string: {smiles}")
    except Exception as e:
        st.error(f"Failed to render ligand: {e}")


def get_base64_mol_image(smiles: str, size: tuple = (160, 160)) -> str:
    """Returns base64 encoded PNG image data URL for table embedding."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        img = Draw.MolToImage(mol, size=size)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return f"data:image/png;base64,{base64.b64encode(buffered.getvalue()).decode()}"
    except Exception:
        return ""


def highlight_rna_sequence(sequence: str, probabilities: list) -> str:
    """
    Takes an RNA sequence and a list of probabilities (0.0 to 1.0) of the same length.
    Returns an HTML string with each nucleotide highlighted using the Fazena orange accent.
    """
    if len(sequence) != len(probabilities):
        return (
            f"Error: Sequence length ({len(sequence)}) does not match "
            f"probabilities length ({len(probabilities)})."
        )
    
    html_output = ""
    for char, prob in zip(sequence, probabilities):
        alpha = max(0.06, prob)
        tooltip = f"{char}: {prob:.3f}"
        font_weight = "bold" if prob > 0.5 else "normal"
        color = "#ffffff" if prob > 0.3 else "rgba(255, 255, 255, 0.55)"
        bg_color = f"rgba(218, 117, 45, {alpha})"
        
        span = (
            f"<span title='{tooltip}' style='"
            f"background-color: {bg_color}; color: {color}; "
            f"font-weight: {font_weight}; padding: 3px 2px; "
            f"border-radius: 3px; margin: 1px; transition: background-color 0.2s; "
            f"display: inline-block; min-width: 14px; text-align: center;'>"
            f"{char}</span>"
        )
        html_output += span
        
    return html_output


def render_complex_3d(
    pdb_str: str,
    height: int = 500,
    receptor_style: str = "cartoon",
    ligand_style: str = "stick"
):
    """
    Renders an interactive 3D viewer for a receptor-ligand complex using py3Dmol.
    Shows the receptor macromolecule as a cartoon and docked ligand as distinct sticks.
    """
    if not pdb_str or not pdb_str.strip():
        st.warning("No 3D PDB structure data available to render.")
        return

    viewer = py3Dmol.view(height=height, width="100%")
    viewer.addModel(pdb_str, 'pdb')

    # Receptor styling (cartoon representation with color gradient)
    viewer.setStyle(
        {'hetflag': False},
        {'cartoon': {'color': 'spectrum', 'opacity': 0.85}}
    )

    # Docked ligand styling (bright orange sticks with element-colored atoms)
    viewer.setStyle(
        {'hetflag': True},
        {
            'stick': {
                'colorscheme': 'orangeCarbon',
                'radius': 0.28
            },
            'sphere': {
                'scale': 0.25,
                'opacity': 0.6
            }
        }
    )

    # Dark background matching Fazena theme (#1D1D1D)
    viewer.setBackgroundColor('#1D1D1D')
    viewer.zoomTo({'hetflag': True})
    
    html = viewer._make_html()
    components.html(html, height=height + 20)


def render_interaction_summary_ui(interaction_summary: Dict[str, Any]):
    """Renders 2D interaction type badges and residue details in Streamlit."""
    total = interaction_summary.get('total_interactions', 0)
    types = interaction_summary.get('types', {})
    residues = interaction_summary.get('interacting_residues', [])

    if total == 0:
        st.info("No close contact interactions detected under default cutoff distance.")
        return

    st.markdown(f"**Total Contacts Detected:** <span style='color: var(--fz-accent, #DA752D); font-weight: bold;'>{total}</span>", unsafe_allow_html=True)
    
    # Interaction type chips
    type_chips = []
    for itype, count in types.items():
        type_chips.append(
            f"<span style='background: rgba(218,117,45,0.15); border: 1px solid rgba(218,117,45,0.3); "
            f"border-radius: 12px; padding: 4px 10px; margin-right: 8px; font-size: 0.85rem; color: #FEFEFE;'>"
            f"<b>{itype}</b>: {count}</span>"
        )
    st.markdown("<div style='margin-bottom: 12px;'>" + "".join(type_chips) + "</div>", unsafe_allow_html=True)

    # Residue breakdown table
    if residues:
        res_df = pd.DataFrame(residues)
        res_df.columns = ['Target Residue', 'Interaction Type']
        st.dataframe(res_df, use_container_width=True, height=200)


class InteractionPlots:
    """Compatibility wrapper for 2D and 3D interaction visualizers."""
    def __init__(self):
        pass

    def plot_2d_interaction(self, complex_pdb_str: str):
        if analyze_interactions:
            summary = analyze_interactions(complex_pdb_str)
            render_interaction_summary_ui(summary)
        else:
            st.info("Interaction analysis engine not available.")

    def plot_3d_interaction(self, complex_pdb_str: str):
        render_complex_3d(complex_pdb_str)

    def generate_2d_interaction_diagram(self, pdbqt_file: str, output_file: str) -> None:
        """Generate a 2D interaction diagram from a PDBQT file."""
        try:
            mol = Chem.MolFromPDBFile(pdbqt_file)
            if mol is None:
                raise ValueError("Invalid PDBQT file")
            Draw.MolToFile(mol, output_file, size=(400, 400), useSVG=True)
        except Exception as e:
            raise RuntimeError(f"Failed to generate 2D interaction diagram: {e}")
