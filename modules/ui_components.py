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

    /* ── Design Tokens ── */
    :root {
        --fz-bg-primary:    #1D1D1D;
        --fz-bg-secondary:  #242424;
        --fz-bg-card:       #282828;
        --fz-bg-card-hover: #2F2F2F;
        --fz-accent:        #DA752D;
        --fz-accent-dim:    #c4672a;
        --fz-accent-glow:   rgba(218, 117, 45, 0.25);
        --fz-accent-glow-s: rgba(218, 117, 45, 0.10);
        --fz-text:          #FEFEFE;
        --fz-text-muted:    rgba(254, 254, 254, 0.60);
        --fz-text-dim:      rgba(254, 254, 254, 0.40);
        --fz-border:        rgba(255, 255, 255, 0.07);
        --fz-border-hover:  rgba(255, 255, 255, 0.13);
        --fz-radius-card:   16px;
        --fz-radius-btn:    6px;
        --fz-radius-input:  8px;
        --fz-font-heading:  'Audiowide', sans-serif;
        --fz-font-body:     'Inter', sans-serif;
        --fz-font-mono:     'JetBrains Mono', monospace;
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
       BUTTONS — Orange Accent
       ═══════════════════════════════════ */
    .stButton > button {
        background-color: transparent !important;
        color: var(--fz-text) !important;
        font-family: var(--fz-font-heading) !important;
        border: 1px solid var(--fz-accent) !important;
        border-radius: var(--fz-radius-btn) !important;
        padding: 0.65rem 2rem !important;
        transition: all 0.35s cubic-bezier(0.165, 0.84, 0.44, 1) !important;
        text-transform: uppercase;
        letter-spacing: 2px;
        font-size: 0.8rem !important;
        position: relative;
        overflow: hidden;
    }

    .stButton > button:hover {
        background-color: var(--fz-accent) !important;
        color: #000000 !important;
        box-shadow: 0 2px 12px var(--fz-accent-glow), inset 0 0 25px rgba(234, 148, 20, 0.3) !important;
        transform: translateY(-1px);
        border-color: var(--fz-accent) !important;
    }

    .stButton > button:active {
        transform: translateY(0);
    }

    /* Download Buttons */
    .stDownloadButton > button {
        background-color: transparent !important;
        color: var(--fz-text) !important;
        font-family: var(--fz-font-heading) !important;
        border: 1px solid rgba(234, 116, 20, 0.6) !important;
        border-radius: var(--fz-radius-btn) !important;
        padding: 0.65rem 2rem !important;
        transition: all 0.35s cubic-bezier(0.165, 0.84, 0.44, 1) !important;
        text-transform: uppercase;
        letter-spacing: 1.5px;
        font-size: 0.75rem !important;
    }

    .stDownloadButton > button:hover {
        background-color: var(--fz-accent) !important;
        color: #000000 !important;
        border-color: var(--fz-accent) !important;
        box-shadow: 0 2px 12px var(--fz-accent-glow) !important;
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
    """Render the Fazena footer branding with creator attribution."""
    st.markdown("""
    <div class="fazena-footer-brand">
        Developed by <a href="https://www.fazena.org" target="_blank">FAZENA</a> &bull;
        <span style="opacity: 0.85; font-size: 0.72rem; letter-spacing: 1px;">Yehor Franchuk, MrEgorin, FAZENA</span>
    </div>
    """, unsafe_allow_html=True)


def render_success_alert(message):
    """Render a themed success alert matching the Fazena design language."""
    ICON_CHECK = '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"/><path d="m9 11 3 3L22 4"/></svg>'
    st.markdown(f"""
    <div class="fz-alert fz-alert-success">
        {ICON_CHECK} {message}
    </div>
    """, unsafe_allow_html=True)
