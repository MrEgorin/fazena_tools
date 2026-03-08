import streamlit as st
import os
import sys

# Ensure modules directory is resolvable
sys.path.append(os.path.join(os.path.dirname(__file__), 'modules'))

from ui_components import (
    inject_custom_css,
    render_header,
    render_metric_card,
    render_section_divider,
    render_capability_card,
    render_footer,
)

st.set_page_config(
    page_title="RNA Hunter | FAZENA",
    page_icon="RNAHunterlogo.svg",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Apply global premium CSS
inject_custom_css()

# ═══════════════════════════════════
# SIDEBAR
# ═══════════════════════════════════
with st.sidebar:
    st.markdown("<div class='sidebar-logo-area'>", unsafe_allow_html=True)
    st.image("RNAHunterlogo.svg", width=100)
    st.markdown(
        "<h2 style='font-family: var(--fz-font-heading); color: var(--fz-text); "
        "font-size: 1.3rem; margin-bottom: 2px; margin-top: 8px; border: none; "
        "padding-left: 0 !important;'>RNA Hunter</h2>",
        unsafe_allow_html=True
    )
    st.markdown(
        "<p style='color: var(--fz-text-dim); font-size: 0.7rem; margin-top: 0; "
        "letter-spacing: 2px; font-family: var(--fz-font-body);'>"
        "POWERED BY RNA FOUNDATION MODELS</p>",
        unsafe_allow_html=True
    )
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("---")

    page = st.radio(
        "Navigation",
        ["Home Dashboard", "Virtual Screening", "Binding Site Prediction"],
        label_visibility="collapsed",
    )

    st.markdown("---")

    # Version and developer info
    st.markdown(
        "<div style='padding: 0.5rem 0;'>"
        "<span class='version-badge'>v2.1</span>"
        "<p style='color: var(--fz-text-dim); font-size: 0.75rem; margin-top: 10px; "
        "font-family: var(--fz-font-body); line-height: 1.6;'>"
        "Developer: <span style='color: var(--fz-accent);'>FAZENA</span><br>"
        "RNAHunter Engine</p>"
        "</div>",
        unsafe_allow_html=True
    )

# ═══════════════════════════════════
# PERSONALIZED BIOTECH ICON DEFINITIONS
# ═══════════════════════════════════
# Discovery / Rocket -> Stylized "In Silico"
ICON_DISCOVERY = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5"/></svg>'

# Virtual Screening -> High Throughput (HTVS)
ICON_SCREENING = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M3 3h18v18H3zM3 9h18M3 15h18M9 3v18M15 3v18"/></svg>'

# Binding Site -> Targeted CADD mapping
ICON_BINDING = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><path d="M12 8v8M8 12h8"/><circle cx="12" cy="12" r="3" fill="currentColor"/></svg>'

# Decoy / Shield -> Multi-layer protection
ICON_PROTECTION = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M12 22s8-4 8-10V5l-8-3-8 3v7c0 6 8 10 8 10zM12 8v8M8 12h8"/></svg>'

# System Status / Engine -> Nucleus/Core
ICON_ENGINE = '<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="4"/><circle cx="12" cy="12" r="9"/><path d="M12 3v2M12 19v2M3 12h2M19 12h2"/></svg>'

# ═══════════════════════════════════
# ROUTING
# ═══════════════════════════════════
if page == "Home Dashboard":
    # Hero Header
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

    # Next-Gen Section
    st.markdown(
        f"<div class='fazena-card'>"
        f"<h3 style='margin-top: 0; display: flex; align-items: center; gap: 12px;'>"
        f"<span style='color: var(--fz-accent);'>{ICON_DISCOVERY}</span> Next-Gen RNA Targeting</h3>"
        f"<p style='color: var(--fz-text-muted); line-height: 1.8; font-size: 0.95rem;'>"
        f"Leveraging foundational bio-models to enable precise structure-agnostic discovery "
        f"of RNA-binding candidates across high-dimensional chemical spaces."
        f"</p></div>",
        unsafe_allow_html=True
    )

    # Capability Cards
    col1, col2, col3 = st.columns(3)
    with col1:
        render_capability_card(
            ICON_SCREENING,
            "Virtual Screening",
            "Screen massive chemical libraries against RNA targets using high-throughput foundational inference."
        )
    with col2:
        render_capability_card(
            ICON_BINDING,
            "Binding Site Prediction",
            "Identify high-probability interaction interfaces with per-residue resolution."
        )
    with col3:
        render_capability_card(
            ICON_PROTECTION,
            "Decoy Enhancement",
            "Validated against complex datasets to ensure robust target-specific selectivity."
        )

    render_section_divider()

    # System Status
    st.markdown(
        f"<h2 style='display: flex; align-items: center; gap: 10px; font-size: 1.1rem;'>"
        f"System Status</h2>",
        unsafe_allow_html=True
    )

    status_col1, status_col2 = st.columns(2)

    weights_path = os.path.join(os.path.dirname(__file__), "SMARTBind_weight")
    weights_found = os.path.exists(weights_path) and len(os.listdir(weights_path)) > 0

    with status_col1:
        if weights_found:
            st.markdown(
                f"<div class='fazena-card' style='border-left: 3px solid #62b1a8;'>"
                f"<span class='status-badge success'>● Active Engine</span>"
                f"<p style='color: var(--fz-text-muted); margin-top: 0.75rem; font-size: 0.9rem;'>"
                f"Computational weights loaded. Inference pipelines operational.</p>"
                f"</div>",
                unsafe_allow_html=True
            )
        else:
            st.markdown(
                f"<div class='fazena-card' style='border-left: 3px solid var(--fz-accent);'>"
                f"<span class='status-badge warning'>◆ Synthetic Simulation</span>"
                f"<p style='color: var(--fz-text-muted); margin-top: 0.75rem; font-size: 0.9rem;'>"
                f"Model weights offline. Running in synthetic mode via <code style='color: var(--fz-accent); "
                f"font-family: var(--fz-font-mono);'>SMARTBind</code> fallback.</p>"
                f"</div>",
                unsafe_allow_html=True
            )

    with status_col2:
        render_metric_card(
            ICON_ENGINE,
            "RNAHunter",
            "Core Engine"
        )

    render_footer()

elif page == "Virtual Screening":
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
