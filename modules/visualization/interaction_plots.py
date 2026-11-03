"""
FAZENA Tools — Molecular & Interaction Visualization Utilities
Provides 2D molecule rendering, 2D protein/RNA-ligand interaction diagrams,
interactive 3D py3Dmol macromolecule & complex viewer with search box preview,
and InteractionPlots compatibility class.
Author: Yehor Franchuk (@MrEgorin), FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

import os
import base64
import html
import json
import numpy as np
from io import BytesIO
from typing import List, Dict, Any, Optional, Tuple

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

import re

try:
    from modules.docking_engine import analyze_interactions
except ImportError:
    analyze_interactions = None


def _mol_to_clean_svg_b64(mol: Chem.Mol, width: int = 280, height: int = 240) -> str:
    """Renders an RDKit Mol to a base64 encoded SVG data URI with clean XML and high contrast."""
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.clearBackground = True
        opts.setBackgroundColour((1.0, 1.0, 1.0, 1.0))
        opts.bondLineWidth = 2
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_raw = drawer.GetDrawingText()
        svg_clean = re.sub(r'<\?xml[^>]*\?>', '', svg_raw).strip()
        b64 = base64.b64encode(svg_clean.encode('utf-8')).decode('ascii')
        return f"data:image/svg+xml;base64,{b64}"
    except Exception:
        return ""


def render_ligand_2d(smiles: str, width: int = 260, height: int = 240):
    """
    Renders a crisp 2D chemical depiction of a ligand given its SMILES string.
    Uses base64 data URI image embedding inside the Fazena dark blue-green card,
    preventing Streamlit markdown DOM sanitization from stripping raw SVG elements.
    """
    try:
        if not smiles or not smiles.strip():
            st.caption("No ligand SMILES provided.")
            return

        mol = Chem.MolFromSmiles(smiles.strip())
        if not mol:
            st.error(f"Invalid SMILES string: {smiles}")
            return

        AllChem.Compute2DCoords(mol)
        data_uri = _mol_to_clean_svg_b64(mol, width=width, height=height)

        # Fallback to PNG if SVG data URI generation failed
        if not data_uri:
            data_uri = get_base64_mol_image(smiles, size=(width, height))

        if data_uri:
            st.markdown(
                f"<div style='text-align: center; background: #1E293B; "
                f"border: 1px solid rgba(14, 165, 233, 0.25); border-radius: 12px; "
                f"padding: 12px; display: inline-block; box-shadow: 0 4px 14px rgba(0,0,0,0.25); max-width: 100%;'>"
                f"<div style='background: #FFFFFF; border-radius: 8px; padding: 8px; display: inline-block;'>"
                f"<img src='{data_uri}' alt='2D Molecule' style='width: {width}px; max-width: 100%; height: auto; display: block;'/>"
                f"</div>"
                f"</div>",
                unsafe_allow_html=True
            )
        else:
            st.warning("Could not generate 2D image for molecule.")
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
    Returns an HTML string with each nucleotide highlighted using the Fazena Emerald-Cyan gradient accent.
    """
    if len(sequence) != len(probabilities):
        return (
            f"Error: Sequence length ({len(sequence)}) does not match "
            f"probabilities length ({len(probabilities)})."
        )

    html_output = ""
    for char, prob in zip(sequence, probabilities):
        alpha = max(0.08, float(prob))
        tooltip = f"{char}: {prob:.3f}"
        font_weight = "bold" if prob > 0.5 else "normal"
        color = "#ffffff" if prob > 0.3 else "rgba(255, 255, 255, 0.65)"
        bg_color = f"rgba(16, 185, 129, {alpha})" if prob > 0.5 else f"rgba(14, 165, 233, {alpha})"

        span = (
            f"<span title='{tooltip}' style='"
            f"background-color: {bg_color}; color: {color}; "
            f"font-weight: {font_weight}; padding: 3px 2px; "
            f"border-radius: 4px; margin: 1px; transition: background-color 0.2s; "
            f"display: inline-block; min-width: 15px; text-align: center; border: 1px solid rgba(255,255,255,0.06);'>"
            f"{char}</span>"
        )
        html_output += span

    return html_output


def render_complex_3d(
    pdb_str: str,
    height: int = 500,
    grid_box: Optional[Dict[str, float]] = None,
    show_box: bool = False,
    receptor_style: str = "cartoon",
    ligand_style: str = "stick",
    is_receptor_only: bool = False
):
    """
    Renders an interactive 3D molecular viewer for macromolecule or complex PDB structures using 3Dmol.js.
    Robustly handles both RNA and protein polymers, applies high-contrast styling for ligands,
    safely checks selection counts before zooming to prevent camera drop-out bugs,
    and optionally renders the 3D grid bounding box wireframe.
    """
    if not pdb_str or not pdb_str.strip():
        st.warning("No 3D PDB structure data available to render.")
        return

    # Prepare escaped PDB string for inline JS
    clean_lines = [l for l in pdb_str.splitlines() if not l.startswith(('CONECT', 'MASTER'))]
    clean_pdb = "\n".join(clean_lines)
    escaped_pdb = clean_pdb.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')

    # Grid box parameters
    box_js = ""
    if show_box and grid_box:
        cx = float(grid_box.get('center_x', 0.0))
        cy = float(grid_box.get('center_y', 0.0))
        cz = float(grid_box.get('center_z', 0.0))
        sx = float(grid_box.get('size_x', 20.0))
        sy = float(grid_box.get('size_y', 20.0))
        sz = float(grid_box.get('size_z', 20.0))
        box_js = f"""
        viewer.addBox({{
            center: {{x: {cx}, y: {cy}, z: {cz}}},
            dimensions: {{w: {sx}, h: {sy}, d: {sz}}},
            color: '#10b981',
            wireframe: true
        }});
        """

    div_id = f"viewer_3d_{abs(hash(pdb_str[:120])) % 1000000}"
    rec_only_js = "true" if is_receptor_only else "false"

    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8"/>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js"></script>
      <style>
        body, html {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; background: #0B1329; }}
        #{div_id} {{ width: 100%; height: {height}px; position: relative; border-radius: 12px; }}
        .badge-ctrl {{ position: absolute; bottom: 8px; left: 8px; z-index: 10; font-family: sans-serif; font-size: 11px; background: rgba(15, 23, 42, 0.85); color: #94A3B8; padding: 4px 10px; border-radius: 6px; border: 1px solid rgba(14,165,233,0.3); font-weight: 600; letter-spacing: 0.5px; }}
        .action-bar {{ position: absolute; top: 8px; right: 8px; z-index: 10; display: flex; gap: 6px; }}
        .btn-act {{ font-family: sans-serif; font-size: 11px; font-weight: 600; background: rgba(15, 23, 42, 0.85); color: #38BDF8; border: 1px solid rgba(14,165,233,0.35); padding: 3px 8px; border-radius: 5px; cursor: pointer; transition: all 0.2s; }}
        .btn-act:hover {{ background: #0ea5e9; color: #ffffff; border-color: #0ea5e9; }}
      </style>
    </head>
    <body>
      <div id="{div_id}">
        <div class="badge-ctrl">FAZENA 3D Interactive View</div>
        <div class="action-bar">
          <button class="btn-act" onclick="window._resetCam_{div_id}()">Center View</button>
          <button class="btn-act" onclick="window._toggleSpin_{div_id}()">Spin</button>
        </div>
      </div>
      <script>
        (function() {{
          var element = document.getElementById("{div_id}");
          var config = {{ backgroundColor: "#0B1329" }};
          var viewer = $3Dmol.createViewer(element, config);

          var pdbData = `{escaped_pdb}`;
          viewer.addModel(pdbData, "pdb");

          // 1. Polymer receptor: cartoon styling + subtle sticks
          viewer.setStyle({{ hetflag: false }}, {{
              cartoon: {{ color: 'spectrum', opacity: 0.92 }},
              stick: {{ radius: 0.12, opacity: 0.35 }}
          }});

          // 2. Nucleic acids (RNA/DNA): cartoon ribbon + sticks for bases to ensure visibility on all formats
          viewer.setStyle({{ nucleic: true }}, {{
              cartoon: {{ color: '#0ea5e9', opacity: 0.90 }},
              stick: {{ radius: 0.16, colorscheme: 'cyanCarbon' }}
          }});

          // 3. Ligand styling: high-contrast emerald sticks + atomic spheres
          viewer.setStyle({{ hetflag: true }}, {{
              stick: {{ colorscheme: 'greenCarbon', radius: 0.28 }},
              sphere: {{ scale: 0.24, opacity: 0.85 }}
          }});
          viewer.setStyle({{ resname: ['UNL', 'LIG', 'MOL', 'DRG', 'UNK'] }}, {{
              stick: {{ colorscheme: 'greenCarbon', radius: 0.28 }},
              sphere: {{ scale: 0.24, opacity: 0.85 }}
          }});

          {box_js}

          // Safe Camera Zoom: verify atoms exist before zooming to selection
          var isReceptorOnly = {rec_only_js};
          function doZoom() {{
              if (!isReceptorOnly) {{
                  var ligAtoms = viewer.selectedAtoms({{ resname: ['UNL', 'LIG', 'MOL', 'DRG', 'UNK'] }});
                  if (ligAtoms && ligAtoms.length > 0) {{
                      viewer.zoomTo({{ resname: ['UNL', 'LIG', 'MOL', 'DRG', 'UNK'] }});
                      return;
                  }}
              }}
              viewer.zoomTo();
          }}
          doZoom();
          viewer.render();

          // Controls
          var spinning = false;
          window._resetCam_{div_id} = function() {{
              viewer.zoomTo();
              viewer.render();
          }};
          window._toggleSpin_{div_id} = function() {{
              spinning = !spinning;
              viewer.spin(spinning ? "y" : false, 1.0);
          }};
        }})();
      </script>
    </body>
    </html>
    """
    components.html(html_code, height=height + 15)


def render_receptor_3d(
    pdb_str: str,
    height: int = 420,
    grid_box: Optional[Dict[str, float]] = None,
    show_box: bool = False
):
    """
    Renders a standalone macromolecule (protein or RNA) with optional grid box wireframe.
    """
    render_complex_3d(pdb_str, height=height, grid_box=grid_box, show_box=show_box, is_receptor_only=True)


# ═════════════════════════════════════════════════════════════════════
# 2D INTERACTION DIAGRAM GENERATION
# ═════════════════════════════════════════════════════════════════════

def generate_2d_interaction_diagram(
    complex_or_pdbqt: str,
    output_file: Optional[str] = None
) -> str:
    """
    Generates a 2D interaction diagram from a complex PDB or ligand PDBQT.
    Extracts the ligand molecule, calculates key binding residue contacts,
    and returns a clean, fully-rendered high-contrast base64 data URI image
    with contact residue annotations that cannot be stripped by markdown parsers.
    """
    content = ""
    if os.path.isfile(complex_or_pdbqt):
        with open(complex_or_pdbqt, 'r') as f:
            content = f.read()
    else:
        content = complex_or_pdbqt

    # 1. Parse ligand
    mol = None
    if "UNL" in content or "HETATM" in content or "ATOM" in content:
        # Extract HETATM records and clean AutoDock types (e.g. 'A' -> 'C')
        het_lines = []
        for l in content.splitlines():
            if l.startswith('HETATM') or ('ATOM' in l and ('UNL' in l or 'LIG' in l or 'MOL' in l or 'DRG' in l)):
                # If column 76-78 contains AutoDock type 'A' or 'OA' or 'NA'
                if len(l) >= 78:
                    c_type = l[76:78].strip()
                    if c_type in ('A', 'C'):
                        l = f"{l[:76]} C{l[78:]}"
                    elif c_type in ('OA', 'O'):
                        l = f"{l[:76]} O{l[78:]}"
                    elif c_type in ('NA', 'N'):
                        l = f"{l[:76]} N{l[78:]}"
                het_lines.append(l)

        if het_lines:
            try:
                clean_het = "\n".join(het_lines) + "\nEND\n"
                mol = Chem.MolFromPDBBlock(clean_het, removeHs=False)
            except Exception:
                mol = None

    if mol is None and ("ROOT" in content or "BRANCH" in content):
        try:
            from meeko import PDBQTMolecule, RDKitMolCreate
            pdbqt_mol = PDBQTMolecule(content)
            mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
            if mols:
                mol = mols[0]
        except Exception:
            mol = None

    if mol is None:
        # Fallback to SMILES in REMARK or default molecule
        for line in content.splitlines():
            if line.startswith("REMARK SMILES"):
                smi = line.split("REMARK SMILES")[-1].strip()
                mol = Chem.MolFromSmiles(smi)
                break

    if mol is None:
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")

    Chem.rdDepictor.Compute2DCoords(mol)

    # 2. Extract contacts
    contacts = []
    if analyze_interactions:
        inter_summary = analyze_interactions(content)
        contacts = inter_summary.get('interacting_residues', [])[:10]

    # 3. Generate clean high-contrast base64 SVG data URI
    data_uri = _mol_to_clean_svg_b64(mol, width=420, height=280)
    if not data_uri:
        # PNG fallback
        buffered = BytesIO()
        img = Draw.MolToImage(mol, size=(420, 280))
        img.save(buffered, format="PNG")
        data_uri = f"data:image/png;base64,{base64.b64encode(buffered.getvalue()).decode()}"

    # Wrap in styled Fazena container with contact badges
    badge_chips = []
    color_map = {
        'Hydrogen Bond': '#0ea5e9',
        'Pi-Stacking': '#10b981',
        'Electrostatic': '#a855f7',
        'Hydrophobic': '#38bdf8',
        'Van der Waals': '#34d399'
    }

    for c in contacts:
        res = c.get('residue', 'Res')
        itype = c.get('type', 'Contact')
        dist = c.get('distance', 0.0)
        c_color = color_map.get(itype, '#0ea5e9')
        badge_chips.append(
            f"<div style='display: inline-flex; align-items: center; gap: 6px; background: rgba(15,23,42,0.85); "
            f"border: 1px solid {c_color}; border-radius: 8px; padding: 4px 9px; margin: 3px; font-size: 0.76rem; color: #F8FAFC;'>"
            f"<span style='font-weight: 700; color: {c_color};'>{res}</span>"
            f"<span style='color: #94A3B8;'>({dist} Å)</span>"
            f"<span style='background: rgba(255,255,255,0.08); border-radius: 4px; padding: 1px 5px; font-size: 0.68rem;'>{itype}</span>"
            f"</div>"
        )

    badges_html = "".join(badge_chips) if badge_chips else "<div style='color: #94A3B8; font-size: 0.8rem;'>No close contacts detected under cutoff.</div>"

    final_svg = f"""
    <!-- <svg> depiction -->
    <div style='background: #1E293B; border: 1px solid rgba(14,165,233,0.25); border-radius: 12px; padding: 14px; text-align: center;'>
        <div style='margin-bottom: 8px; font-size: 0.82rem; font-weight: 700; color: #10b981; letter-spacing: 1px;'>
            2D LIGAND DEPICTION & BINDING POCKET CONTACTS
        </div>
        <div style='display: flex; justify-content: center; background: #0B1329; border-radius: 8px; padding: 10px;'>
            <div style='background: #FFFFFF; border-radius: 8px; padding: 8px; display: inline-block;'>
                <img src='{data_uri}' alt='2D Interaction Topology' style='max-width: 100%; height: auto; display: block;'/>
            </div>
        </div>
        <div style='margin-top: 12px; text-align: left;'>
            <div style='font-size: 0.74rem; color: #94A3B8; text-transform: uppercase; font-weight: 600; letter-spacing: 0.5px; margin-bottom: 6px;'>Key Contacting Residues:</div>
            {badges_html}
        </div>
    </div>
    """

    if output_file:
        try:
            drawer = rdMolDraw2D.MolDraw2DSVG(440, 340)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            with open(output_file, 'w') as f:
                f.write(drawer.GetDrawingText())
        except Exception:
            pass

    return final_svg

    if output_file:
        with open(output_file, 'w') as f:
            f.write(base_svg)

    return final_svg


def render_complex_2d(complex_pdb_str: str):
    """Renders the 2D interaction diagram directly in Streamlit."""
    diagram_html = generate_2d_interaction_diagram(complex_pdb_str)
    st.markdown(diagram_html, unsafe_allow_html=True)


def render_interaction_summary_ui(interaction_summary: Dict[str, Any], complex_pdb_str: Optional[str] = None):
    """Renders 2D interaction depiction, badges, and residue details in Streamlit."""
    total = interaction_summary.get('total_interactions', 0)
    types = interaction_summary.get('types', {})
    residues = interaction_summary.get('interacting_residues', [])

    if complex_pdb_str:
        render_complex_2d(complex_pdb_str)
        st.markdown("<br>", unsafe_allow_html=True)

    if total == 0:
        st.info("No close contact interactions detected under default cutoff distance.")
        return

    st.markdown(f"**Total Contacts Detected:** <span style='color: #10b981; font-weight: bold;'>{total}</span>", unsafe_allow_html=True)

    color_map = {
        'Hydrogen Bond': ('#0ea5e9', 'rgba(14, 165, 233, 0.15)'),
        'Pi-Stacking': ('#10b981', 'rgba(16, 185, 129, 0.15)'),
        'Electrostatic': ('#a855f7', 'rgba(168, 85, 247, 0.15)'),
        'Hydrophobic': ('#38bdf8', 'rgba(56, 189, 248, 0.15)'),
        'Van der Waals': ('#34d399', 'rgba(52, 211, 153, 0.15)')
    }

    type_chips = []
    for itype, count in types.items():
        border_c, bg_c = color_map.get(itype, ('#0ea5e9', 'rgba(14,165,233,0.15)'))
        type_chips.append(
            f"<span style='background: {bg_c}; border: 1px solid {border_c}; "
            f"border-radius: 12px; padding: 4px 10px; margin-right: 8px; font-size: 0.82rem; color: #F8FAFC;'>"
            f"<b>{itype}</b>: {count}</span>"
        )
    st.markdown("<div style='margin-bottom: 12px; margin-top: 6px;'>" + "".join(type_chips) + "</div>", unsafe_allow_html=True)

    if residues:
        res_df = pd.DataFrame(residues)
        col_names = {'residue': 'Residue', 'chain': 'Chain', 'type': 'Interaction Type', 'distance': 'Distance (Å)', 'rec_atom': 'Receptor Atom', 'lig_atom': 'Ligand Atom'}
        res_df = res_df.rename(columns=col_names)
        st.dataframe(res_df, use_container_width=True, height=220)


class InteractionPlots:
    """Compatibility wrapper for 2D and 3D interaction visualizers."""
    def __init__(self):
        pass

    def plot_2d_interaction(self, complex_pdb_str: str):
        if analyze_interactions:
            summary = analyze_interactions(complex_pdb_str)
            render_interaction_summary_ui(summary, complex_pdb_str=complex_pdb_str)
        else:
            st.info("Interaction analysis engine not available.")

    def plot_3d_interaction(self, complex_pdb_str: str):
        render_complex_3d(complex_pdb_str)

    def generate_2d_interaction_diagram(self, pdbqt_file: str, output_file: str) -> None:
        generate_2d_interaction_diagram(pdbqt_file, output_file=output_file)
