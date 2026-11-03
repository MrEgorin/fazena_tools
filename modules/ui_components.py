import streamlit as st

def inject_custom_css():
    st.markdown("""
    <style>
    /* ═══════════════════════════════════════════════════════════
       FAZENA Design System — RNA Hunter Beta
       Ported from fazena.org (Next.js) to Streamlit
       ═══════════════════════════════════════════════════════════ */

    /* ── Google Fonts ── */
    @import url('https://fonts.googleapis.com/css2?family=Audiowide&family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

    /* ── Design Tokens (fazena.org Blue-Green Theme) ── */
    :root {
        --fz-bg-primary:    #0B1329;
        --fz-bg-secondary:  #0F172A;
        --fz-bg-card:       #1E293B;
        --fz-bg-card-hover: #24344D;
        --fz-accent:        #0ea5e9;
        --fz-accent-emerald:#10b981;
        --fz-accent-dim:    #0284c7;
        --fz-accent-glow:   rgba(14, 165, 233, 0.28);
        --fz-accent-glow-s: rgba(16, 185, 129, 0.18);
        --fz-text:          #F8FAFC;
        --fz-text-muted:    rgba(248, 250, 252, 0.70);
        --fz-text-dim:      rgba(248, 250, 252, 0.45);
        --fz-border:        rgba(14, 165, 233, 0.18);
        --fz-border-hover:  rgba(16, 185, 129, 0.35);
        --fz-radius-card:   16px;
        --fz-radius-btn:    8px;
        --fz-radius-input:  8px;
        --fz-font-heading:  'Audiowide', sans-serif;
        --fz-font-body:     'Inter', sans-serif;
        --fz-font-mono:     'JetBrains Mono', monospace;
        --fz-gradient:      linear-gradient(135deg, #10b981 0%, #0ea5e9 100%);
    }

    /* ── Text Selection ── */
    ::selection {
        background-color: var(--fz-accent);
        color: #ffffff;
    }
    ::-moz-selection {
        background-color: var(--fz-accent);
        color: #ffffff;
    }

    /* ── App Background ── */
    .stApp {
        background: var(--fz-bg-primary) !important;
        color: var(--fz-text);
        font-family: var(--fz-font-body);
    }

    .stApp > header {
        background: transparent !important;
    }

    /* ── Scrollbars ── */
    ::-webkit-scrollbar {
        width: 6px;
        height: 6px;
    }
    ::-webkit-scrollbar-track {
        background: var(--fz-bg-primary);
    }
    ::-webkit-scrollbar-thumb {
        background: rgba(255, 255, 255, 0.08);
        border-radius: 3px;
    }
    ::-webkit-scrollbar-thumb:hover {
        background: rgba(255, 255, 255, 0.15);
    }

    /* ═══════════════════════════════════
       SIDEBAR
       ═══════════════════════════════════ */
    section[data-testid="stSidebar"] {
        background: var(--fz-bg-primary) !important;
        border-right: 1px solid var(--fz-border) !important;
    }

    section[data-testid="stSidebar"] .stMarkdown p {
        color: var(--fz-text-muted);
        font-family: var(--fz-font-body);
    }

    section[data-testid="stSidebar"] hr {
        border-color: var(--fz-border) !important;
    }

    /* Sidebar logo area */
    .sidebar-logo-area {
        padding: 0.5rem 0 0.5rem 0;
        margin-bottom: 0.5rem;
    }

    /* ═══════════════════════════════════
       TYPOGRAPHY
       ═══════════════════════════════════ */
    h1, h2, h3, h4, h5, h6 {
        font-family: var(--fz-font-heading) !important;
        letter-spacing: 1.5px !important;
        text-transform: uppercase;
        color: var(--fz-text) !important;
    }

    h1 {
        font-size: 2rem !important;
        border-left: 4px solid var(--fz-accent);
        padding-left: 1rem !important;
        margin-bottom: 0.5rem !important;
    }

    h2 {
        font-size: 1.4rem !important;
        margin-bottom: 1rem !important;
    }

    h3 {
        font-size: 1.1rem !important;
    }

    p, li, span {
        font-family: var(--fz-font-body);
    }

    /* ═══════════════════════════════════
       FAZENA CARDS — Glassmorphism
       ═══════════════════════════════════ */
    .fazena-card {
        background: var(--fz-bg-card);
        border: 1px solid var(--fz-border);
        border-radius: var(--fz-radius-card);
        padding: 2rem;
        backdrop-filter: blur(12px);
        -webkit-backdrop-filter: blur(12px);
        margin-bottom: 1.5rem;
        transition: all 0.35s cubic-bezier(0.165, 0.84, 0.44, 1);
        position: relative;
        overflow: hidden;
    }

    .fazena-card::before {
        content: '';
        position: absolute;
        inset: -1px;
        border-radius: var(--fz-radius-card);
        background: linear-gradient(135deg, var(--fz-accent-glow-s), transparent 60%);
        opacity: 0;
        transition: opacity 0.35s ease;
        pointer-events: none;
        z-index: 0;
    }

    .fazena-card:hover {
        border-color: var(--fz-border-hover);
        background: var(--fz-bg-card-hover);
        box-shadow: 0 20px 50px rgba(0, 0, 0, 0.4);
    }

    .fazena-card:hover::before {
        opacity: 1;
    }

    /* Card with orange glow on hover */
    .fazena-card-glow:hover {
        border-color: var(--fz-accent);
        box-shadow: 0 0 30px var(--fz-accent-glow);
    }

    /* ═══════════════════════════════════
       METRIC CARDS
       ═══════════════════════════════════ */
    .metric-card {
        background: var(--fz-bg-card);
        border: 1px solid var(--fz-border);
        border-radius: var(--fz-radius-card);
        padding: 1.5rem;
        text-align: center;
        transition: all 0.3s ease;
    }

    .metric-card:hover {
        border-color: var(--fz-border-hover);
        transform: translateY(-2px);
        box-shadow: 0 10px 30px rgba(0, 0, 0, 0.3);
    }

    .metric-card .metric-icon {
        font-size: 1.8rem;
        margin-bottom: 0.5rem;
    }

    .metric-card .metric-value {
        font-family: var(--fz-font-heading);
        font-size: 1.6rem;
        color: var(--fz-accent);
        margin-bottom: 0.25rem;
    }

    .metric-card .metric-label {
        font-family: var(--fz-font-body);
        font-size: 0.8rem;
        color: var(--fz-text-muted);
        letter-spacing: 1px;
        text-transform: uppercase;
    }

    /* ═══════════════════════════════════
       BUTTONS — Fazena Blue-Green
       ═══════════════════════════════════ */
    .stButton > button {
        background: linear-gradient(135deg, rgba(16, 185, 129, 0.15) 0%, rgba(14, 165, 233, 0.15) 100%) !important;
        color: var(--fz-text) !important;
        font-family: var(--fz-font-heading) !important;
        border: 1px solid var(--fz-accent-emerald) !important;
        border-radius: var(--fz-radius-btn) !important;
        padding: 0.65rem 2rem !important;
        transition: all 0.3s cubic-bezier(0.165, 0.84, 0.44, 1) !important;
        text-transform: uppercase;
        letter-spacing: 2px;
        font-size: 0.82rem !important;
        position: relative;
        overflow: hidden;
        box-shadow: 0 4px 15px rgba(16, 185, 129, 0.15) !important;
    }

    .stButton > button:hover {
        background: linear-gradient(135deg, #10b981 0%, #0ea5e9 100%) !important;
        color: #ffffff !important;
        box-shadow: 0 6px 20px rgba(14, 165, 233, 0.35), inset 0 0 15px rgba(255, 255, 255, 0.2) !important;
        transform: translateY(-2px);
        border-color: #0ea5e9 !important;
    }

    .stButton > button:active {
        transform: translateY(0);
    }

    /* Download Buttons */
    .stDownloadButton > button {
        background: rgba(14, 165, 233, 0.1) !important;
        color: var(--fz-text) !important;
        font-family: var(--fz-font-heading) !important;
        border: 1px solid rgba(14, 165, 233, 0.5) !important;
        border-radius: var(--fz-radius-btn) !important;
        padding: 0.65rem 2rem !important;
        transition: all 0.3s cubic-bezier(0.165, 0.84, 0.44, 1) !important;
        text-transform: uppercase;
        letter-spacing: 1.5px;
        font-size: 0.78rem !important;
    }

    .stDownloadButton > button:hover {
        background: linear-gradient(135deg, #0ea5e9 0%, #10b981 100%) !important;
        color: #ffffff !important;
        border-color: #10b981 !important;
        box-shadow: 0 4px 18px rgba(14, 165, 233, 0.3) !important;
        transform: translateY(-2px);
    }

    /* ═══════════════════════════════════
       INPUT FIELDS
       ═══════════════════════════════════ */
    .stTextInput > div > div > input,
    .stTextArea > div > textarea,
    .stNumberInput > div > div > input {
        background-color: var(--fz-bg-primary) !important;
        color: var(--fz-text) !important;
        border: 1px solid var(--fz-border) !important;
        border-radius: var(--fz-radius-input) !important;
        padding: 10px 14px !important;
        font-family: var(--fz-font-body) !important;
        transition: border-color 0.3s ease, box-shadow 0.3s ease !important;
    }

    .stTextInput > div > div > input:focus,
    .stTextArea > div > textarea:focus,
    .stNumberInput > div > div > input:focus {
        border-color: var(--fz-accent) !important;
        box-shadow: 0 0 12px var(--fz-accent-glow) !important;
        outline: none !important;
    }

    /* Input labels */
    .stTextInput label,
    .stTextArea label,
    .stNumberInput label,
    .stSelectbox label,
    .stMultiSelect label,
    .stFileUploader label,
    .stCheckbox label,
    .stRadio label {
        color: var(--fz-text-muted) !important;
        font-family: var(--fz-font-body) !important;
        font-size: 0.85rem !important;
        letter-spacing: 0.5px !important;
    }

    /* ═══════════════════════════════════
       SELECTBOX & MULTISELECT
       ═══════════════════════════════════ */
    .stSelectbox > div > div,
    .stMultiSelect > div > div {
        background-color: var(--fz-bg-primary) !important;
        border: 1px solid var(--fz-border) !important;
        border-radius: var(--fz-radius-input) !important;
        color: var(--fz-text) !important;
    }

    /* ═══════════════════════════════════
       RADIO BUTTONS
       ═══════════════════════════════════ */
    .stRadio > div {
        background: transparent;
    }

    .stRadio > div > label {
        color: var(--fz-text-muted) !important;
        transition: color 0.2s ease;
    }

    .stRadio > div > label:hover {
        color: var(--fz-text) !important;
    }

    /* ═══════════════════════════════════
       CHECKBOX
       ═══════════════════════════════════ */
    .stCheckbox > label > span {
        color: var(--fz-text-muted) !important;
    }

    /* ═══════════════════════════════════
       FILE UPLOADER
       ═══════════════════════════════════ */
    .stFileUploader > div {
        background-color: var(--fz-bg-card) !important;
        border: 1px dashed var(--fz-border-hover) !important;
        border-radius: var(--fz-radius-card) !important;
        transition: border-color 0.3s ease !important;
    }

    .stFileUploader > div:hover {
        border-color: var(--fz-accent) !important;
    }

    /* ═══════════════════════════════════
       TABS
       ═══════════════════════════════════ */
    .stTabs [data-baseweb="tab-list"] {
        background-color: transparent !important;
        gap: 20px;
        border-bottom: 1px solid var(--fz-border) !important;
    }

    .stTabs [data-baseweb="tab"] {
        font-family: var(--fz-font-heading) !important;
        color: var(--fz-text-muted) !important;
        background-color: transparent !important;
        font-size: 0.8rem;
        letter-spacing: 1px;
        text-transform: uppercase;
        padding: 0.5rem 1rem;
        transition: color 0.3s ease;
    }

    .stTabs [data-baseweb="tab"]:hover {
        color: var(--fz-text) !important;
    }

    .stTabs [aria-selected="true"] {
        color: var(--fz-accent) !important;
        border-bottom-color: var(--fz-accent) !important;
    }

    /* ═══════════════════════════════════
       EXPANDER
       ═══════════════════════════════════ */
    .streamlit-expanderHeader {
        background-color: var(--fz-bg-card) !important;
        border: 1px solid var(--fz-border) !important;
        border-radius: var(--fz-radius-input) !important;
        color: var(--fz-text-muted) !important;
        font-family: var(--fz-font-body) !important;
        transition: border-color 0.3s ease !important;
    }

    .streamlit-expanderHeader:hover {
        border-color: var(--fz-border-hover) !important;
        color: var(--fz-text) !important;
    }

    details {
        background-color: var(--fz-bg-card) !important;
        border: 1px solid var(--fz-border) !important;
        border-radius: var(--fz-radius-card) !important;
    }

    details summary {
        color: var(--fz-text-muted) !important;
        font-family: var(--fz-font-body) !important;
    }

    details summary:hover {
        color: var(--fz-text) !important;
    }

    /* ═══════════════════════════════════
       METRICS (st.metric)
       ═══════════════════════════════════ */
    div[data-testid="stMetric"] {
        background: var(--fz-bg-card);
        border: 1px solid var(--fz-border);
        border-radius: var(--fz-radius-card);
        padding: 1rem 1.25rem;
        transition: all 0.3s ease;
    }

    div[data-testid="stMetric"]:hover {
        border-color: var(--fz-border-hover);
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.3);
    }

    div[data-testid="stMetric"] label {
        color: var(--fz-text-muted) !important;
        font-family: var(--fz-font-body) !important;
        font-size: 0.8rem !important;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
        color: var(--fz-accent) !important;
        font-family: var(--fz-font-heading) !important;
        font-size: 1.5rem !important;
    }

    /* ═══════════════════════════════════
       PROGRESS BAR
       ═══════════════════════════════════ */
    .stProgress > div > div > div > div {
        background-color: var(--fz-accent) !important;
        background-image: linear-gradient(90deg, var(--fz-accent), #f7931a) !important;
    }

    .stProgress > div > div {
        background-color: rgba(255, 255, 255, 0.05) !important;
    }

    /* ═══════════════════════════════════
       DATAFRAME / TABLE
       ═══════════════════════════════════ */
    .stDataFrame {
        border: 1px solid var(--fz-border) !important;
        border-radius: var(--fz-radius-input) !important;
        background: var(--fz-bg-card) !important;
        overflow: hidden;
    }

    /* ═══════════════════════════════════
       ALERTS (st.success, st.error, etc.)
       ═══════════════════════════════════ */
    .stAlert {
        border-radius: var(--fz-radius-input) !important;
        font-family: var(--fz-font-body) !important;
    }

    div[data-testid="stNotification"] {
        border-radius: var(--fz-radius-input) !important;
    }

    /* ═══════════════════════════════════
       SPINNER
       ═══════════════════════════════════ */
    .stSpinner > div {
        border-top-color: var(--fz-accent) !important;
    }

    /* ═══════════════════════════════════
       HR / DIVIDER
       ═══════════════════════════════════ */
    hr {
        border-color: var(--fz-border) !important;
    }

    /* ═══════════════════════════════════
       RNA SEQUENCE CONTAINER
       ═══════════════════════════════════ */
    .rna-sequence-container {
        font-family: var(--fz-font-mono);
        background: var(--fz-bg-primary);
        padding: 1.5rem;
        border-radius: var(--fz-radius-input);
        border: 1px solid var(--fz-border);
        color: var(--fz-accent);
        word-break: break-all;
        line-height: 1.8;
        font-size: 1rem;
        margin-top: 1rem;
        box-shadow: inset 0 0 15px rgba(0, 0, 0, 0.6);
    }

    /* ═══════════════════════════════════
       SECTION DIVIDER
       ═══════════════════════════════════ */
    .section-divider {
        border: none;
        height: 1px;
        background: linear-gradient(90deg, transparent, var(--fz-border-hover), transparent);
        margin: 2rem 0;
    }

    /* ═══════════════════════════════════
       STATUS BADGE
       ═══════════════════════════════════ */
    .status-badge {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        padding: 4px 12px;
        border-radius: 100px;
        font-size: 0.75rem;
        font-family: var(--fz-font-body);
        letter-spacing: 0.5px;
        text-transform: uppercase;
    }

    .status-badge.success {
        background: rgba(34, 197, 94, 0.1);
        color: #22c55e;
        border: 1px solid rgba(34, 197, 94, 0.2);
    }

    .status-badge.warning {
        background: rgba(234, 116, 20, 0.1);
        color: var(--fz-accent);
        border: 1px solid rgba(234, 116, 20, 0.2);
    }

    .status-badge.error {
        background: rgba(239, 68, 68, 0.1);
        color: #ef4444;
        border: 1px solid rgba(239, 68, 68, 0.2);
    }

    /* ═══════════════════════════════════
       SIDEBAR NAV ITEMS
       ═══════════════════════════════════ */
    .sidebar-nav-item {
        display: flex;
        align-items: center;
        gap: 10px;
        padding: 10px 14px;
        border-radius: 10px;
        color: var(--fz-text-muted);
        font-family: var(--fz-font-body);
        font-size: 0.9rem;
        text-decoration: none;
        transition: all 0.25s ease;
        cursor: pointer;
        margin-bottom: 2px;
    }

    .sidebar-nav-item:hover {
        color: var(--fz-text);
        background: rgba(255, 255, 255, 0.03);
    }

    .sidebar-nav-item.active {
        color: var(--fz-text);
        background: rgba(255, 255, 255, 0.04);
    }

    .sidebar-nav-item.active::before {
        content: '';
        position: absolute;
        left: 0;
        top: 50%;
        transform: translateY(-50%);
        width: 3px;
        height: 60%;
        background: var(--fz-accent);
        border-radius: 0 2px 2px 0;
    }

    /* ═══════════════════════════════════
       FOOTER BRANDING
       ═══════════════════════════════════ */
    .fazena-footer-brand {
        text-align: center;
        padding: 2rem 0 1rem;
        color: var(--fz-text-dim);
        font-size: 0.75rem;
        letter-spacing: 2px;
        text-transform: uppercase;
    }

    .fazena-footer-brand a {
        color: var(--fz-accent);
        text-decoration: none;
        font-weight: 600;
        transition: color 0.3s ease;
    }

    .fazena-footer-brand a:hover {
        color: #f7931a;
    }

    /* ═══════════════════════════════════
       HERO TITLE (gradient)
       ═══════════════════════════════════ */
    .hero-title {
        font-family: var(--fz-font-heading);
        font-size: 2.2rem;
        letter-spacing: 2px;
        text-transform: uppercase;
        background: linear-gradient(135deg, #ffffff 0%, var(--fz-accent) 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        margin-bottom: 0.25rem;
        border: none;
        padding: 0 !important;
    }

    /* ═══════════════════════════════════
       CAPABILITY CARD
       ═══════════════════════════════════ */
    .capability-card {
        background: var(--fz-bg-card);
        border: 1px solid var(--fz-border);
        border-radius: var(--fz-radius-card);
        padding: 1.5rem;
        transition: all 0.35s cubic-bezier(0.165, 0.84, 0.44, 1);
        position: relative;
        overflow: hidden;
        height: 100%;
    }

    .capability-card::before {
        content: '';
        position: absolute;
        inset: -1px;
        border-radius: var(--fz-radius-card);
        background: linear-gradient(135deg, var(--fz-accent-glow-s), transparent 60%);
        opacity: 0;
        transition: opacity 0.35s ease;
        pointer-events: none;
    }

    .capability-card:hover {
        border-color: var(--fz-border-hover);
        background: var(--fz-bg-card-hover);
        box-shadow: 0 20px 50px rgba(0, 0, 0, 0.4);
        transform: translateY(-2px);
    }

    .capability-card:hover::before {
        opacity: 1;
    }

    .capability-card .cap-icon {
        display: flex;
        align-items: center;
        justify-content: center;
        width: 48px;
        height: 48px;
        border-radius: 12px;
        background: var(--fz-accent-glow-s);
        border: 1px solid rgba(234, 116, 20, 0.2);
        font-size: 1.4rem;
        margin-bottom: 1rem;
        box-shadow: 0 0 15px rgba(247, 147, 26, 0.1);
        transition: transform 0.3s ease;
    }

    .capability-card:hover .cap-icon {
        transform: scale(1.1);
    }

    .capability-card .cap-title {
        font-family: var(--fz-font-heading);
        font-size: 0.9rem;
        color: var(--fz-text);
        margin-bottom: 0.5rem;
        letter-spacing: 1px;
        text-transform: uppercase;
    }

    .capability-card .cap-desc {
        font-family: var(--fz-font-body);
        font-size: 0.85rem;
        color: var(--fz-text-muted);
        line-height: 1.6;
    }

    /* ═══════════════════════════════════
       ANIMATIONS
       ═══════════════════════════════════ */
    @keyframes fadeInUp {
        from {
            opacity: 0;
            transform: translateY(20px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }

    @keyframes scanLine {
        0% { transform: translateX(-100%); }
        100% { transform: translateX(100%); }
    }

    .scan-line {
        position: relative;
        overflow: hidden;
    }

    .scan-line::after {
        content: '';
        position: absolute;
        bottom: 0;
        left: 0;
        width: 100%;
        height: 2px;
        background: linear-gradient(90deg, transparent, var(--fz-accent), transparent);
        animation: scanLine 3s ease-in-out infinite;
    }

    .fade-in-up {
        animation: fadeInUp 0.6s ease-out forwards;
    }

    /* Version badge */
    .version-badge {
        display: inline-flex;
        align-items: center;
        gap: 4px;
        padding: 3px 10px;
        border-radius: 100px;
        font-size: 0.65rem;
        font-family: var(--fz-font-mono);
        letter-spacing: 0.5px;
        background: rgba(234, 116, 20, 0.08);
        color: var(--fz-accent);
        border: 1px solid rgba(234, 116, 20, 0.15);
    }

    /* Alert styles */
    .fz-alert {
        padding: 0.75rem 1.25rem;
        border-radius: var(--fz-radius-input);
        font-size: 0.9rem;
        margin: 1rem 0;
        display: flex;
        align-items: center;
        gap: 10px;
    }
    .fz-alert-success {
        background: rgba(218, 117, 45, 0.1);
        color: var(--fz-accent);
        border: 1px solid rgba(218, 117, 45, 0.2);
    }
    .rank-badge {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 28px;
        height: 28px;
        border-radius: 50%;
        background: rgba(218, 117, 45, 0.1);
        color: var(--fz-accent);
        border: 1px solid rgba(218, 117, 45, 0.3);
        font-family: var(--fz-font-mono);
        font-size: 0.75rem;
        font-weight: 600;
    }
    </style>
    """, unsafe_allow_html=True)


def render_header(title, subtitle):
    """Render a premium section header matching the main website SectionTitle pattern."""
    st.markdown(f"<h1>{title}</h1>", unsafe_allow_html=True)
    st.markdown(
        f"<p style='color: var(--fz-text-muted); font-size: 1.05rem; "
        f"margin-bottom: 2rem; line-height: 1.7;'>{subtitle}</p>",
        unsafe_allow_html=True
    )


def render_metric_card(icon, value, label):
    """Render a premium metric display card."""
    st.markdown(f"""
    <div class="metric-card">
        <div class="metric-icon">{icon}</div>
        <div class="metric-value">{value}</div>
        <div class="metric-label">{label}</div>
    </div>
    """, unsafe_allow_html=True)


def render_section_divider():
    """Render a styled section divider with gradient fade."""
    st.markdown('<hr class="section-divider">', unsafe_allow_html=True)


def render_capability_card(icon, title, description):
    """Render a capability/feature card matching SingleFeature from main site."""
    st.markdown(f"""
    <div class="capability-card">
        <div class="cap-icon">{icon}</div>
        <div class="cap-title">{title}</div>
        <div class="cap-desc">{description}</div>
    </div>
    """, unsafe_allow_html=True)


def render_footer():
    """Render the Fazena footer branding with creator, team, and affiliation attribution."""
    st.markdown("""
    <div style='margin-top: 3.5rem; padding-top: 2rem; border-top: 1px solid var(--fz-border); text-align: center;'>
        <div style='display: inline-flex; flex-wrap: wrap; justify-content: center; align-items: center; gap: 12px; margin-bottom: 0.8rem;'>
            <a href="https://www.fazena.org" target="_blank" style="text-decoration: none;">
                <span style="background: rgba(16, 185, 129, 0.12); border: 1px solid rgba(16, 185, 129, 0.35); color: #10b981; padding: 5px 12px; border-radius: 20px; font-size: 0.82rem; font-weight: 600; display: inline-flex; align-items: center; gap: 6px;">
                    🧬 FAZENA Team
                </span>
            </a>
            <span style="background: rgba(14, 165, 233, 0.12); border: 1px solid rgba(14, 165, 233, 0.35); color: #38bdf8; padding: 5px 12px; border-radius: 20px; font-size: 0.82rem; font-weight: 500; display: inline-flex; align-items: center; gap: 6px;">
                👨‍💻 Yehor Franchuk (MrEgorin)
            </span>
            <span style="background: rgba(168, 85, 247, 0.12); border: 1px solid rgba(168, 85, 247, 0.35); color: #c084fc; padding: 5px 12px; border-radius: 20px; font-size: 0.82rem; font-weight: 500; display: inline-flex; align-items: center; gap: 6px;">
                🏛️ Institute of Cell Biology NASU
            </span>
        </div>
        <div style="color: var(--fz-text-dim); font-size: 0.76rem; letter-spacing: 0.5px; line-height: 1.6;">
            Next-Generation Modular Drug Discovery Foundry &bull; <a href="https://www.fazena.org" target="_blank" style="color: var(--fz-accent); text-decoration: none;">fazena.org</a><br>
            All AI inferences and physics-based AutoDock Vina calculations execute strictly on local models and trained weights.
        </div>
    </div>
    """, unsafe_allow_html=True)


def render_docking_parameter_controls(default_grid_box: dict, is_rna: bool = False, key_prefix: str = "dock", allow_gpu: bool = True) -> dict:
    """
    Renders flexible docking parameter controls with an intuitive interface.
    Allows configuring box center (x, y, z), box dimensions (size_x, size_y, size_z),
    docking engine (AutoDock Vina CPU vs Vina-GPU CUDA acceleration),
    search exhaustiveness, number of poses, energy range, random seed, and scoring function.
    """
    from modules.config import HARDWARE

    with st.expander("AutoDock Vina Engine & Search Box Configuration", expanded=False):
        st.markdown(
            f"<p style='color: var(--fz-text-muted); font-size: 0.88rem; margin-bottom: 1rem;'>"
            f"Configure search box coordinates and physics simulation settings for <b>{'RNA Target' if is_rna else 'Protein Target'}</b>."
            f"</p>",
            unsafe_allow_html=True
        )

        # Hardware Engine Selection
        has_cuda = HARDWARE.get('has_cuda', False)
        cuda_name = HARDWARE.get('cuda_device_name', 'NVIDIA GPU')
        vram_gb = HARDWARE.get('gpu_vram_gb', 0.0)

        engine_mode = "vina"
        gpu_threads = 8000
        vina_gpu_bin = ""

        if allow_gpu and has_cuda:
            st.markdown(
                f"<div style='background: rgba(16, 185, 129, 0.08); border: 1px solid rgba(16, 185, 129, 0.25); "
                f"border-radius: 8px; padding: 6px 12px; margin-bottom: 12px; font-size: 0.82rem; color: #34d399; display: flex; align-items: center; justify-content: space-between;'>"
                f"<span>Hardware Detected: <b>{cuda_name}</b> ({vram_gb} GB VRAM)</span>"
                f"<span style='background: rgba(16, 185, 129, 0.2); padding: 2px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600;'>CUDA ACTIVE</span>"
                f"</div>",
                unsafe_allow_html=True
            )
            engine_choice = st.radio(
                "Docking Acceleration Engine",
                [
                    "AutoDock Vina-GPU (NVIDIA CUDA Acceleration)",
                    "AutoDock Vina (CPU Multithreading)"
                ],
                index=0,
                key=f"{key_prefix}_engine_choice"
            )
            if "Vina-GPU" in engine_choice:
                engine_mode = "vina_gpu"
                col_g1, col_g2 = st.columns([1, 1])
                with col_g1:
                    gpu_threads = st.slider(
                        "GPU Workgroup Threads (--thread)",
                        min_value=1000,
                        max_value=16000,
                        value=8000,
                        step=1000,
                        help="Number of GPU threads allocated for parallel Monte Carlo search iterations on RTX 5070 Ti.",
                        key=f"{key_prefix}_gthreads"
                    )
                with col_g2:
                    vina_gpu_bin = st.text_input(
                        "Custom Vina-GPU Binary (Auto-detected if empty)",
                        value="",
                        placeholder="e.g. /usr/local/bin/AutoDock-Vina-GPU-2-1",
                        key=f"{key_prefix}_gbin"
                    )
            else:
                engine_mode = "vina"

        col_c1, col_c2, col_c3 = st.columns(3)
        with col_c1:
            center_x = st.number_input("Center X (Å)", value=float(default_grid_box.get('center_x', 0.0)), step=1.0, format="%.2f", key=f"{key_prefix}_cx")
        with col_c2:
            center_y = st.number_input("Center Y (Å)", value=float(default_grid_box.get('center_y', 0.0)), step=1.0, format="%.2f", key=f"{key_prefix}_cy")
        with col_c3:
            center_z = st.number_input("Center Z (Å)", value=float(default_grid_box.get('center_z', 0.0)), step=1.0, format="%.2f", key=f"{key_prefix}_cz")

        col_s1, col_s2, col_s3 = st.columns(3)
        with col_s1:
            size_x = st.slider("Size X (Å)", min_value=10.0, max_value=80.0, value=float(default_grid_box.get('size_x', 24.0)), step=1.0, key=f"{key_prefix}_sx")
        with col_s2:
            size_y = st.slider("Size Y (Å)", min_value=10.0, max_value=80.0, value=float(default_grid_box.get('size_y', 24.0)), step=1.0, key=f"{key_prefix}_sy")
        with col_s3:
            size_z = st.slider("Size Z (Å)", min_value=10.0, max_value=80.0, value=float(default_grid_box.get('size_z', 24.0)), step=1.0, key=f"{key_prefix}_sz")

        st.markdown("<hr style='border: none; border-top: 1px solid var(--fz-border); margin: 0.8rem 0;'>", unsafe_allow_html=True)

        col_p1, col_p2, col_p3 = st.columns(3)
        with col_p1:
            exhaustiveness = st.slider("Exhaustiveness", min_value=1, max_value=32, value=8, step=1, help="Higher values increase thoroughness of Monte Carlo sampling (8 is standard, 16-32 for deep search).", key=f"{key_prefix}_ex")
            n_poses = st.slider("Number of Poses (Modes)", min_value=1, max_value=20, value=5, step=1, key=f"{key_prefix}_nposes")
        with col_p2:
            energy_range = st.number_input("Energy Range (kcal/mol)", min_value=1.0, max_value=10.0, value=3.0, step=0.5, help="Maximum energy difference between best and worst pose.", key=f"{key_prefix}_er")
            seed = st.number_input("Random Seed (0 = Random)", min_value=0, max_value=999999, value=42, step=1, help="Fixed seed ensures exact numerical reproducibility.", key=f"{key_prefix}_seed")
        with col_p3:
            import os
            max_cpus = os.cpu_count() or 4
            cpu_threads = st.slider("CPU Threads", min_value=1, max_value=min(max_cpus, 16), value=min(max_cpus, 8), step=1, key=f"{key_prefix}_cpu")
            sf_name = st.selectbox("Scoring Function", ["vina", "vinardo"], index=0, key=f"{key_prefix}_sf")

        show_box_wireframe = st.checkbox("Show 3D Grid Search Box in Viewer", value=True, key=f"{key_prefix}_showbox")

    grid_box = {
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z,
        'size_x': size_x,
        'size_y': size_y,
        'size_z': size_z
    }

    return {
        'grid_box': grid_box,
        'exhaustiveness': exhaustiveness,
        'n_poses': n_poses,
        'energy_range': energy_range,
        'seed': seed,
        'cpu': cpu_threads,
        'sf_name': sf_name,
        'engine_mode': engine_mode,
        'gpu_threads': gpu_threads,
        'vina_gpu_bin': vina_gpu_bin,
        'show_box': show_box_wireframe
    }


def render_success_alert(message):
    """Render a themed success alert matching the Fazena design language."""
    ICON_CHECK = '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/><path d="m9 11 3 3L22 4"/></svg>'
    st.markdown(f"""
    <div class="fz-alert fz-alert-success">
        {ICON_CHECK} {message}
    </div>
    """, unsafe_allow_html=True)
