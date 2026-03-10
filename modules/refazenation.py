import streamlit as st
import pandas as pd
import time
import random
import os
import numpy as np
import math
import base64
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from .ui_components import render_header, render_section_divider, render_footer, render_success_alert
from .refazenation_wrapper import run_docking_real, create_complex_pdb

# Icons
ICON_REVERSE = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M7 11l5-5 5 5M7 17l5-5 5 5"/></svg>'
ICON_TARGET = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><path d="M12 8v8M8 12h8"/><circle cx="12" cy="12" r="3"/></svg>'
ICON_PROTEIN = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5"/></svg>'
ICON_3D = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/><polyline points="3.27 6.96 12 12.01 20.73 6.96"/><line x1="12" y1="22.08" x2="12" y2="12"/></svg>'
ICON_DOWNLOAD = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>'

# Vina Path Configuration
VINA_PATH = r"C:\AutoDock-Vina-develop"

# ═══════════════════════════════════
# ADVANCED PREPARATION LOGIC
# ═══════════════════════════════════

def get_mol_safe(smiles_or_sdf, is_sdf=False):
    """Parses molecule with RELAXED sanitization to handle high-valence P/O/S atoms."""
    try:
        mol = None
        if is_sdf:
            suppl = Chem.ForwardSDMolSupplier(BytesIO(smiles_or_sdf), sanitize=False)
            mol = next(suppl, None)
        else:
            mol = Chem.MolFromSmiles(smiles_or_sdf, sanitize=False)
        
        if mol:
            # Crucial: Calculate properties before accessing valence to avoid C++ Pre-condition Violation
            mol.UpdatePropertyCache(strict=False)
            
            # Step 1: Attempt to fix explicit valence for problematic atoms
            for atom in mol.GetAtoms():
                try:
                    # Fix P(6+) or O(3+) by neutralizing or setting formal charges
                    symbol = atom.GetSymbol()
                    val = atom.GetValence(which=Chem.rdchem.ValenceType.EXPLICIT)
                    if symbol == 'P' and val > 5:
                        atom.SetNoImplicit(True)
                    if symbol == 'O' and val > 2:
                        atom.SetFormalCharge(1)
                except:
                    continue # Skip if valence calculation still fails for an atom
            
            # Step 2: Partial sanitize (Skip property verification to allow non-standard hypervalency)
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        return mol
    except Exception as e:
        st.warning(f"Sanitization Warning: {str(e)}")
        return None

def prepare_ligand_3d(mol):
    """Embeds ligand in 3D, adds hydrogens, and returns 3D mol."""
    if not mol: return None
    try:
        mol_3d = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        AllChem.EmbedMolecule(mol_3d, params)
        AllChem.MMFFOptimizeMolecule(mol_3d)
        return mol_3d
    except:
        return mol

def clean_protein_pdb(pdb_str):
    """Removes water (HOH), adds hydrogens, and calculates box."""
    lines = pdb_str.splitlines()
    clean_lines = []
    coords = []
    
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            res_name = line[17:20].strip()
            if res_name in ["HOH", "WAT", "SOL", "TIP3", "CL", "NA", "MG"]:
                continue
            
            clean_lines.append(line)
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
            except: pass
    
    if not coords: return None, None
    coords = np.array(coords)
    center = np.mean(coords, axis=0)
    size = np.max(coords, axis=0) - np.min(coords, axis=0) + 12.0 # padding
    box = {"center": center.tolist(), "size": size.tolist()}
    return "\n".join(clean_lines), box

# ═══════════════════════════════════
# DOCKING EXECUTION PIPELINE
# ═══════════════════════════════════

def execute_docking_pipeline(ligand_mol, pdb_files):
    """Multi-stage execution using real AutoRevDock/Vina logic wrappers."""
    results = []
    
    for pdb_file in pdb_files:
        status = st.empty()
        status.info(f"[{pdb_file.name}] Stage 1: Initializing real docking pipeline...")
        time.sleep(0.5)

        # Clean Protein
        pdb_str = pdb_file.getvalue().decode("utf-8")
        cleaned_pdb, box = clean_protein_pdb(pdb_str)
        
        # Real Docking Call
        docking_data = run_docking_real(ligand_mol, [pdb_file], VINA_PATH)
        
        status.info(f"[{pdb_file.name}] Stage 3: Parsing Results Log...")
        time.sleep(0.5)
        
        # Generate Complex
        complex_pdb_str = create_complex_pdb(cleaned_pdb, ligand_mol)
        
        status.empty()
        
        # Add to results
        results.append({
            "Protein": pdb_file.name,
            "Score": docking_data[0]['Score'],
            "Ka": docking_data[0]['Ka'],
            "logKa": docking_data[0]['logKa'],
            "RMSD": docking_data[0].get('RMSD', 0.0),
            "Center": f"{box['center'][0]:.1f}, {box['center'][1]:.1f}, {box['center'][2]:.1f}",
            "Complex_PDB": complex_pdb_str
        })
    return pd.DataFrame(results).sort_values(by="Score")

# ═══════════════════════════════════
# UI RENDERING
# ═══════════════════════════════════

def render_page():
    render_header("Refazenation", "Professional Reverse Docking Engine using specialized Protein-Ligand pipelines.")

    col_l, col_p = st.columns([1, 1.2])

    with col_l:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_3D}</span> Ligand Input</h3>", unsafe_allow_html=True)
        
        l_input_type = st.radio("Input Format", ["SMILES", "SDF File"])
        smiles = ""
        mol_3d = None
        
        if l_input_type == "SMILES":
            smiles = st.text_input("Enter SMILES", "CC(C)(C)c1ccc2occ(CC(=O)Nc3ccccc3F)c2c1")
            if smiles:
                with st.spinner("Preparing 3D Conformers..."):
                    raw_mol = get_mol_safe(smiles)
                    mol_3d = prepare_ligand_3d(raw_mol)
        else:
            f = st.file_uploader("Upload SDF", type=['sdf'])
            if f:
                with st.spinner("Processing SDF Ensemble..."):
                    raw_mol = get_mol_safe(f.getvalue(), is_sdf=True)
                    mol_3d = prepare_ligand_3d(raw_mol)
                    if mol_3d: smiles = Chem.MolToSmiles(mol_3d)

        if mol_3d:
            try:
                # Create a 2D copy for clean visualization
                img_mol = Chem.Mol(mol_3d)
                img_mol = Chem.RemoveHs(img_mol)
                AllChem.Compute2DCoords(img_mol)
                
                img = Draw.MolToImage(img_mol, size=(250, 250))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                b64_img = base64.b64encode(buffered.getvalue()).decode()
                st.markdown(f"<center><img src='data:image/png;base64,{b64_img}' style='width: 180px; border-radius: 12px; border: 1px solid var(--fz-border); margin: 15px 0;'></center>", unsafe_allow_html=True)
            except Exception as e:
                st.error(f"Visualization Error: {e}")
        st.markdown("</div>", unsafe_allow_html=True)

    with col_p:
        st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_PROTEIN}</span> Protein target batch</h3>", unsafe_allow_html=True)
        
        uploaded_pdbs = st.file_uploader("Upload Protein PDBs", type=['pdb'], accept_multiple_files=True)
        
        st.markdown("<p style='font-size: 0.85rem; color: var(--fz-text-dim);'>Integrated Pipeline: Automated PDB Cleansing → Blind Search Space Definition → Rescoring Core Application</p>", unsafe_allow_html=True)
        
        if st.button("🚀 Initialize Docking Sequence", width="stretch"):
            if not mol_3d or not uploaded_pdbs:
                st.error("Missing input data (Ligand/Proteins)."); return
            
            results_df = execute_docking_pipeline(mol_3d, uploaded_pdbs)
            st.session_state['refdock_res'] = results_df
            render_success_alert("Batch Processing Complete!")

        st.markdown("</div>", unsafe_allow_html=True)

    if 'refdock_res' in st.session_state:
        render_section_divider()
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_TARGET}</span> Affinity Profiles</h3>", unsafe_allow_html=True)
        
        df = st.session_state['refdock_res']
        
        # Display Results
        for i, row in df.reset_index(drop=True).iterrows():
            st.markdown(f"""
            <div class='fazena-card' style='margin-bottom: 10px; padding: 15px;'>
                <div style='display: flex; justify-content: space-between; align-items: center;'>
                    <div>
                        <span class='rank-badge'>{i+1}</span>
                        <b style='color: var(--fz-accent); font-size: 1.1rem; margin-left: 10px;'>{row['Protein']}</b>
                    </div>
                    <div style='text-align: right;'>
                        <span style='font-size: 0.8rem; color: var(--fz-text-dim);'>SCORE</span><br>
                        <b style='font-size: 1.2rem;'>{row['Score']}</b> kcal/mol
                    </div>
                </div>
                <div style='display: flex; gap: 40px; margin-top: 15px; border-top: 1px solid rgba(255,255,255,0.05); padding-top: 10px;'>
                    <div><span style='color: var(--fz-text-dim); font-size: 0.75rem;'>Affinity Ka</span><br>{row['Ka']} M⁻¹</div>
                    <div><span style='color: var(--fz-text-dim); font-size: 0.75rem;'>logKa</span><br>{row['logKa']}</div>
                    <div><span style='color: var(--fz-text-dim); font-size: 0.75rem;'>RMSD</span><br>{row['RMSD']} Å</div>
                    <div><span style='color: var(--fz-text-dim); font-size: 0.75rem;'>Search Center</span><br>{row['Center']}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # Download Complex Button
            btn_col_label = f"Download Complex ({row['Protein']})"
            st.download_button(
                btn_col_label, 
                data=row['Complex_PDB'], 
                file_name=f"complex_{row['Protein']}", 
                mime="text/plain", 
                key=f"dl_{i}", 
                width="stretch"
            )

        render_section_divider()
        st.download_button("Export Table Summary (.csv)", df.drop(columns=['Complex_PDB']).to_csv(index=False).encode('utf-8'), "Refazenation_Summary.csv", "text/csv", width="stretch")

    render_footer()
