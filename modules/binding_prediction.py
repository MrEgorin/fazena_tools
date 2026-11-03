import streamlit as st
import time
import pandas as pd
from rdkit import Chem
from .smartbind_wrapper import run_binding_site_prediction, is_engine_active
from .visualization import highlight_rna_sequence, render_ligand_2d
from .ui_components import render_header, render_section_divider, render_footer, render_success_alert

# Personalized SVG Icons
ICON_HELIX = ('<svg width="22" height="22" viewBox="0 0 512 512" fill="none" xmlns="http://www.w3.org/2000/svg">'
    '<path d="M150.815 0L160.901 0C163.181 0.913013 165.648 1.1956 167.927 2.05723C176.627 5.3469 183.973 11.9088 187.816 20.4291C189.266 23.6033 190.218 26.9817 190.64 30.4456C191.333 35.9937 190.903 42.4391 191.01 48.0944C191.254 60.969 190.609 74.692 191.188 87.4759C192.454 112.813 199.326 137.553 211.308 159.913C224.218 184.108 242.749 204.85 265.343 220.392C275.329 227.236 284.55 232.164 294.558 239.707C351.742 282.685 386.923 348.809 390.612 420.248C391.293 432.191 391.007 445.089 391.016 457.143C391.021 463.433 391.35 474.232 390.785 480.139C390.509 483.358 389.782 486.522 388.625 489.538C385.312 498.136 378.697 505.135 370.265 508.951C368.183 509.893 362.052 511.475 361.189 512L351.099 512C349.619 511.283 347.389 510.957 345.78 510.476C335.105 507.335 325.843 498.304 322.616 487.606C320.144 479.409 321.155 469.254 321.01 460.653C320.831 450.072 321.177 439.414 320.957 428.82C320.583 414.857 318.523 400.992 314.82 387.524C304.443 349.851 281.548 316.833 249.903 293.907C243.32 289.15 236.658 285.361 229.889 280.975C222.927 276.462 216.212 271.577 209.774 266.343C177.896 240.396 153.227 206.688 138.137 168.455C128.405 143.821 122.728 117.774 121.326 91.3243C120.469 74.5769 121.473 57.8148 120.991 40.922C120.502 23.7554 125.88 9.87253 142.714 2.58691C145.25 1.48965 148.096 1.08011 150.622 0.0787148L150.815 0Z" fill="var(--fz-accent)"/>'
    '<path d="M221.111 31.0136L312.021 30.9886L341.696 30.9785C349.101 30.9769 359.908 29.6548 365.757 34.6551C368.835 37.2954 370.697 41.0807 370.91 45.13C371.149 49.6923 369.587 53.1805 366.613 56.5187C364.759 58.3778 361.563 60.3503 358.981 60.6721C353.891 61.3067 348.097 61.0929 342.935 61.0904L316.26 61.0755L220.973 61.0834C221.293 51.5109 220.527 40.3575 221.111 31.0136Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M221.72 91.1214C221.936 91.1056 222.152 91.0922 222.368 91.0812C231.114 90.6785 242.219 90.9882 251.088 90.9862L305.673 90.9804L339.501 90.983C345.298 90.984 351.533 90.8305 357.306 91.2127C360.713 91.4383 363.732 92.7428 366.255 95.0991C369.135 97.8071 370.816 101.553 370.925 105.505C371.06 110.062 369.382 113.514 366.352 116.781C361.184 121.234 357.245 121.182 350.715 121.092C348.116 121.056 345.366 121.083 342.755 121.079L314.453 121.076L259.052 121.082C248.6 121.086 237.589 121.289 227.181 120.999C226.543 118.116 225.645 115.204 224.953 112.304C223.572 106.608 222.554 100.829 221.906 95.0037C221.792 93.9653 221.348 91.958 221.72 91.1214Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M241.224 151.116C249.114 150.741 258.987 150.991 266.998 150.988L313.441 150.98L342.164 150.982C347.194 150.982 352.655 150.828 357.639 151.244C360.879 151.515 363.855 152.86 366.255 155.098C369.169 157.835 370.854 161.633 370.928 165.63C371.01 170.289 369.28 173.751 366.112 177.005C365.816 177.251 365.516 177.493 365.213 177.73C363.177 179.315 360.743 180.586 358.167 180.797C352.441 181.266 346.329 181.089 340.568 181.089L308.791 181.08L280.79 181.094C276.355 181.107 268.224 181.37 264.138 180.852C260.807 177.894 240.28 154.411 241.224 151.116Z" fill="var(--fz-text)" opacity="0.8"/>'
    '<path d="M154.737 451.091C160.718 450.8 167.799 450.986 173.882 450.986L207.636 450.983C235.189 450.982 263.497 450.642 290.997 451.054L291.012 481.04C286.723 481.108 282.357 481.074 278.062 481.08L195.232 481.085L169.447 481.095C162.741 481.098 152.641 482.154 147.143 478.133C143.87 475.745 141.705 472.134 141.141 468.123C139.92 459.096 145.993 452.365 154.737 451.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '<path d="M154.735 391.091C163.595 390.74 173.937 390.983 182.881 390.985L234.142 390.986L234.142 390.986L267 390.982C271.757 390.982 280.269 390.739 284.701 391.183C286.626 397.819 288.039 403.531 289.146 410.38C289.695 413.772 289.956 417.937 290.584 421.054L277.125 421.078L195.269 421.087L169.391 421.094C162.79 421.096 152.703 422.12 147.268 418.224C143.978 415.877 141.779 412.298 141.172 408.303C139.832 399.221 145.91 392.377 154.735 391.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '<path d="M154.729 331.091C160.486 330.857 166.841 330.985 172.666 330.983L203.319 330.983L230.267 330.987C233.097 330.987 246.882 330.46 248.603 332.024C255.001 337.838 260.828 345.369 265.828 352.444C267.724 355.126 269.404 358.525 271.561 360.935C267.033 361.047 262.504 361.094 257.975 361.077L189.283 361.088L167.552 361.096C161.449 361.098 152.728 361.943 147.585 358.444C144.22 356.168 141.923 352.624 141.219 348.623C139.65 339.487 145.812 332.398 154.729 331.091Z" fill="var(--fz-accent)" opacity="0.6"/>'
    '</svg>')
ICON_DOCKING = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="9" stroke-dasharray="2 2"/><path d="M12 2v20M2 12h20"/><circle cx="12" cy="12" r="3" fill="var(--fz-accent)"/><path d="M16 16l4 4M8 8L4 4" stroke-width="2"/></svg>'
ICON_MAP = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M3 6l6-3 6 3 6-3v15l-6 3-6-3-6 3V6zM9 3v15M15 6v15"/></svg>'
ICON_TABLE = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect width="18" height="18" x="3" y="3" rx="2"/><path d="M3 9h18M3 15h18M12 3v18"/></svg>'
ICON_EXPORT = '<svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>'

def render_page():
    render_header("RNA Hunter | Binding Site Search", "Map atomic-scale interaction probabilities across RNA sequences using cross-modal inference.")
    
    # Removed engine status as per user request
    
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_HELIX}</span> Target Sequence</h3>", unsafe_allow_html=True)
        rna_seq = st.text_area("RNA Sequence", height=100, placeholder="GGCGCUUGCUGGGUAAAA...")
    with col2:
        st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_DOCKING}</span> Ligand Probe</h3>", unsafe_allow_html=True)
        input_type = st.radio("Ligand Input Type", ["Paste SMILES", "Upload SDF"])
        ligand_smiles = ""
        if input_type == "Paste SMILES": ligand_smiles = st.text_input("Ligand SMILES", "CC(C)(C)c1ccc2occ(CC(=O)Nc3ccccc3F)c2c1")
        elif uploaded_file := st.file_uploader("Upload Ligand (.sdf)", type=['sdf']):
            if mol := next(Chem.ForwardSDMolSupplier(uploaded_file), None):
                ligand_smiles = Chem.MolToSmiles(mol); st.info(f"Loaded: {ligand_smiles}")
        
        if ligand_smiles:
            st.markdown("<p style='color: var(--fz-text-muted); font-size: 0.85rem; text-transform: uppercase; letter-spacing: 1px; margin-top: 1rem;'>Ligand 2D Topology</p>", unsafe_allow_html=True)
            render_ligand_2d(ligand_smiles)
    st.markdown("</div>", unsafe_allow_html=True)

    if st.button("Generate Binding Map", width="stretch"):
        if not rna_seq or not ligand_smiles:
            st.error("Please provide both RNA Sequence and Ligand Probe."); return
        with st.spinner("Decoding Interaction interfaces..."):
            time.sleep(1); probabilities = run_binding_site_prediction(rna_seq, ligand_smiles)
            render_success_alert("Prediction Complete!")
            
            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_MAP}</span> Interaction Probability Heatmap</h3>", unsafe_allow_html=True)
            st.markdown(f"<div class='rna-sequence-container'>{highlight_rna_sequence(rna_seq, probabilities)}</div>", unsafe_allow_html=True)
            
            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_TABLE}</span> Quantized Probabilities</h3>", unsafe_allow_html=True)
            df_results = pd.DataFrame({"Position": range(1, len(rna_seq) + 1), "Nucleotide": list(rna_seq), "Probability": probabilities})
            st.dataframe(df_results.style.background_gradient(cmap='YlOrRd', subset=['Probability']), width="stretch", height=300)

            render_section_divider()
            st.markdown(f"<h3 style='display: flex; align-items: center; gap: 10px;'><span style='color: var(--fz-accent);'>{ICON_EXPORT}</span> Export Analytics</h3>", unsafe_allow_html=True)
            df_filtered = df_results[df_results['Probability'] > 0.6]
            if not df_filtered.empty:
                st.download_button(f"Download {len(df_filtered)} Hotspots (.csv)", df_filtered.to_csv(index=False).encode('utf-8'), "hotspots.csv", "text/csv", width="stretch")
            else: st.info("No binding hotspots detected above 0.6 threshold.")
    render_footer()