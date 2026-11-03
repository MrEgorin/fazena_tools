"""
FAZENA Tools — Modular Drug Discovery Platform
Unified Streamlit Application for RNA Hunter, Prozena, and Refazenation.
Author: Yehor Franchuk, MrEgorin, FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

import streamlit as st

# Configure page layout and metadata
st.set_page_config(
    page_title="FAZENA Tools | Modular Drug Discovery Platform",
    layout="wide",
    initial_sidebar_state="expanded"
)

from modules.ui_components import inject_custom_css, render_footer, render_header
from modules.config import HARDWARE, DOCKING_CONFIG
import modules.virtual_screening as rna_vs
import modules.binding_prediction as rna_bsp
import modules.prozena as prozena_page
import modules.refazenation as refazenation_page

# Inject custom FAZENA design system CSS (Blue-Green fazena.org theme)
inject_custom_css()

def render_home_dashboard():
    """Renders the central platform dashboard overview."""
    st.markdown("""
    <div style='text-align: center; padding: 2.5rem 1rem 1.5rem;'>
        <div style='display: inline-flex; align-items: center; gap: 8px; background: rgba(16, 185, 129, 0.12); border: 1px solid rgba(16, 185, 129, 0.3); border-radius: 20px; padding: 4px 14px; margin-bottom: 1rem;'>
            <span style='color: #10b981; font-size: 0.8rem; font-weight: 700; letter-spacing: 1.5px;'>FAZENA BIOCOMPUTATIONAL PLATFORM</span>
        </div>
        <h1 class='hero-title' style='font-size: 2.8rem; margin-bottom: 0.6rem; background: linear-gradient(135deg, #10b981 0%, #0ea5e9 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
            FAZENA TOOLS
        </h1>
        <p style='color: var(--fz-text-muted); font-size: 1.15rem; max-width: 780px; margin: 0 auto 1.5rem; line-height: 1.7;'>
            Next-generation modular bio-computational platform for RNA and protein drug discovery,
            powered by contrastive machine learning (ConPLex) and AutoDock Vina physics engines.
        </p>
        <div style='display: inline-flex; flex-wrap: wrap; justify-content: center; gap: 10px; margin-bottom: 1.5rem;'>
            <span style='background: rgba(14, 165, 233, 0.1); border: 1px solid rgba(14, 165, 233, 0.25); color: #38bdf8; font-size: 0.78rem; padding: 3px 12px; border-radius: 12px; font-weight: 600;'>
                Offline Local Compute & Physics
            </span>
            <span style='background: rgba(16, 185, 129, 0.1); border: 1px solid rgba(16, 185, 129, 0.25); color: #34d399; font-size: 0.78rem; padding: 3px 12px; border-radius: 12px; font-weight: 600;'>
                Zero External LLM Dependencies
            </span>
            <span style='background: rgba(168, 85, 247, 0.1); border: 1px solid rgba(168, 85, 247, 0.25); color: #c084fc; font-size: 0.78rem; padding: 3px 12px; border-radius: 12px; font-weight: 600;'>
                Institute of Cell Biology NASU
            </span>
        </div>
    </div>
    """, unsafe_allow_html=True)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
        <div class='capability-card' style='height: 100%; border-top: 3px solid #10b981;'>
            <div style='font-size: 0.8rem; font-weight: 700; color: #10b981; letter-spacing: 1px; margin-bottom: 6px; text-transform: uppercase;'>MODULE 01</div>
            <div class='cap-title' style='color: #10b981;'>RNA Hunter</div>
            <div class='cap-desc'>
                Structure-agnostic virtual screening and nucleotide-resolution binding site discovery 
                against target RNA architectures with physics-grounded pharmacophores and AutoDock Vina docking.
            </div>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class='capability-card' style='height: 100%; border-top: 3px solid #0ea5e9;'>
            <div style='font-size: 0.8rem; font-weight: 700; color: #0ea5e9; letter-spacing: 1px; margin-bottom: 6px; text-transform: uppercase;'>MODULE 02</div>
            <div class='cap-title' style='color: #0ea5e9;'>Prozena</div>
            <div class='cap-desc'>
                High-throughput screening against protein targets using the ConPLex contrastive co-embedding model,
                integrated with Lipinski/PAINS ADMET filtering and automated AutoDock Vina/Vina-GPU docking.
            </div>
        </div>
        """, unsafe_allow_html=True)

    with col3:
        st.markdown("""
        <div class='capability-card' style='height: 100%; border-top: 3px solid #38bdf8;'>
            <div style='font-size: 0.8rem; font-weight: 700; color: #38bdf8; letter-spacing: 1px; margin-bottom: 6px; text-transform: uppercase;'>MODULE 03</div>
            <div class='cap-title' style='color: #38bdf8;'>Refazenation</div>
            <div class='cap-desc'>
                Target prediction across high-coverage human and pathogen proteomes and multi-target reverse docking
                across therapeutic classes to identify polypharmacology and repurposing candidates.
            </div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br><hr class='section-divider'><br>", unsafe_allow_html=True)

    # System Status & Hardware Card
    st.markdown("<div class='fazena-card'>", unsafe_allow_html=True)
    st.markdown("### System & Hardware Diagnostics")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Platform OS", HARDWARE.get('os', 'Linux'))
    c2.metric("Python Runtime", HARDWARE.get('python_version', '3.10'))
    dock_status = "Active (Python API)" if DOCKING_CONFIG.get('has_vina_python') else "Available"
    c3.metric("AutoDock Vina", dock_status)
    if HARDWARE.get('has_cuda') and HARDWARE.get('recommended_device') == 'cuda':
        cap_str = f" [{HARDWARE.get('compute_capability')}]" if HARDWARE.get('compute_capability') else ""
        cuda_stat = f"{HARDWARE.get('cuda_device_name')} ({HARDWARE.get('gpu_vram_gb', 0)} GB){cap_str}"
    elif HARDWARE.get('has_cuda'):
        cuda_stat = f"{HARDWARE.get('cuda_device_name')} (CPU Mode)"
    else:
        cuda_stat = "CPU Mode"
    c4.metric("Hardware Compute", cuda_stat)
    st.markdown("</div>", unsafe_allow_html=True)

    render_footer()

def main():
    # ═════════════════════════════════════════════════════════════════
    # SIDEBAR NAVIGATION
    # ═════════════════════════════════════════════════════════════
    with st.sidebar:
        st.markdown("""
        <div style='padding: 1rem 0 0.5rem; text-align: left;'>
            <div style='font-family: var(--fz-font-heading); font-size: 1.45rem; letter-spacing: 1.5px; background: linear-gradient(135deg, #10b981 0%, #0ea5e9 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                FAZENA <span style='color: var(--fz-accent);'>TOOLS</span>
            </div>
            <div style='font-size: 0.72rem; color: #10b981; font-weight: 600; text-transform: uppercase; letter-spacing: 1px;'>
                Biocomputational Drug Foundry
            </div>
        </div>
        <hr style='border: none; border-bottom: 1px solid var(--fz-border); margin: 0.5rem 0 1rem 0;'>
        """, unsafe_allow_html=True)

        product = st.selectbox(
            "Select Module",
            ["RNA Hunter", "Prozena", "Refazenation", "Platform Overview"],
            index=0
        )

        st.markdown("<br>", unsafe_allow_html=True)

        if product == "RNA Hunter":
            page = st.radio(
                "RNA Hunter Navigation",
                ["Virtual Target Screening", "Binding Site Search"],
                label_visibility="visible"
            )
        elif product == "Prozena":
            page = st.radio(
                "Prozena Navigation",
                ["HTS Protein Screening & Docking"],
                label_visibility="visible"
            )
        elif product == "Refazenation":
            page = st.radio(
                "Refazenation Navigation",
                ["Target Prediction & Reverse Docking"],
                label_visibility="visible"
            )
        else:
            page = "Overview"

        st.markdown("""
        <div style='margin-top: 3.5rem; padding-top: 1.5rem; border-top: 1px solid var(--fz-border); font-size: 0.74rem; color: var(--fz-text-dim); line-height: 1.6;'>
            <b style="color: #10b981;">FAZENA Discovery Suite</b><br>
            Developed by <a href="https://www.fazena.org" target="_blank" style="color: var(--fz-accent); text-decoration: none; font-weight: 600;">FAZENA</a><br>
            <span>Lead: <b>Yehor Franchuk (MrEgorin)</b></span><br>
            <span style="opacity: 0.9; color: #a855f7;">Institute of Cell Biology NASU</span>
        </div>
        """, unsafe_allow_html=True)

    # ═════════════════════════════════════════════════════════════════
    # PAGE ROUTER
    # ═════════════════════════════════════════════════════════════
    if product == "RNA Hunter":
        if page == "Virtual Target Screening":
            rna_vs.render_page()
        elif page == "Binding Site Search":
            rna_bsp.render_page()

    elif product == "Prozena":
        prozena_page.render_page()

    elif product == "Refazenation":
        refazenation_page.render_page()

    else:
        render_home_dashboard()

if __name__ == "__main__":
    main()
