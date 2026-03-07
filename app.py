import streamlit as st
import os
import sys

# Ensure modules directory is resolvable
sys.path.append(os.path.join(os.path.dirname(__file__), 'modules'))

from ui_components import inject_custom_css, render_header

st.set_page_config(
    page_title="RNA Hunter", 
    page_icon="RNAHunterlogo.jpg",
    layout="wide", 
    initial_sidebar_state="expanded"
)

# Apply global premium CSS
inject_custom_css()

# Sidebar Navigation
with st.sidebar:
    st.image("RNAHunterlogo.jpg", width=120)
    st.markdown("<h2 style='color: #ea7414; font-size: 1.5rem; margin-bottom: 0px;'>RNA Hunter</h2>", unsafe_allow_html=True)
    st.markdown("<p style='color: rgba(255, 255, 255, 0.6); font-size: 0.8em; margin-top: 0px; letter-spacing: 1px;'>POWERED BY RNA FOUNDATION MODELS</p>", unsafe_allow_html=True)
    st.markdown("---")
    
    page = st.radio(
        "Navigation",
        ["Home Dashboard", "Virtual Screening", "Binding Site Prediction"],
        label_visibility="collapsed"
    )
    
    st.markdown("---")
    st.markdown("<p style='color: #8b9bb4; font-size: 0.8em;'>Developer: FAZENA<br>Version: 2.1 (RNAHunter Engine)</p>", unsafe_allow_html=True)

# Routing
if page == "Home Dashboard":
    render_header("Welcome to RNA Hunter", "Advanced Structure-Agnostic RNA Ligand Discovery")
    
    st.markdown("""
    ### 🚀 Next-Gen RNA Targeting
    This revamped version of RNA Hunter leverages a foundation model approach.
    
    **Capabilities:**
    * **Virtual Screening:** Screen thousands of small molecules against an RNA sequence without needing custom 3D structures.
    * **Binding Site Prediction:** Identify exactly which residues are most likely to interact with a specific drug candidate.
    * **Ligand-Specific Decoy Enhancement:** Validated against HARIBOSS & RNAmigos datasets for robust screening.
    """)
        
    st.markdown("---")
    st.markdown("### 📊 System Status")
    
    # Check for SMARTBind weights placeholder
    weights_path = os.path.join(os.path.dirname(__file__), "SMARTBind_weight")
    if os.path.exists(weights_path) and len(os.listdir(weights_path)) > 0:
        st.success("✅ Model Weights Found.")
    else:
        st.warning("⚠️ Model Weights Missing. The app will use mocked data until the weights are downloaded from Zenodo to `SMARTBind_weight/`.")
        
    st.markdown("<br><p style='text-align: center; color: rgba(255, 255, 255, 0.4); font-size: 0.8em; letter-spacing: 2px;'>DEVELOPED BY <a href='https://www.fazena.org' target='_blank' style='color: #ea7414; text-decoration: none; font-weight: bold;'>FAZENA</a></p>", unsafe_allow_html=True)
        
elif page == "Virtual Screening":
    # Refrain from directly importing until we create the virtual_screening.py
    try:
        from modules import virtual_screening
        virtual_screening.render_page()
    except Exception as e:
        st.error(f"Module not found or error loading: {e}")

elif page == "Binding Site Prediction":
    try:
        from modules import binding_prediction
        binding_prediction.render_page()
    except Exception as e:
        st.error(f"Module not found or error loading: {e}")
