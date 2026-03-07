import streamlit as st

def inject_custom_css():
    st.markdown("""
    <style>
    /* Google Fonts */
    @import url('https://fonts.googleapis.com/css2?family=Michroma&family=Inter:wght@300;400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');

    /* Global Variables - Fazena.org Inspired */
    :root {
        --primary-bg: #000000;
        --secondary-bg: #0a0a0a;
        --accent-orange: #ea7414;
        --accent-glow: rgba(234, 116, 20, 0.3);
        --text-main: #ffffff;
        --text-muted: rgba(255, 255, 255, 0.6);
        --card-bg: rgba(255, 255, 255, 0.03);
        --glass-border: rgba(255, 255, 255, 0.1);
    }

    /* Overall App Background with Vignette */
    .stApp {
        background: radial-gradient(circle at center, #0a0a0a 0%, #000000 100%);
        color: var(--text-main);
        font-family: 'Inter', sans-serif;
    }

    /* Sidebar Customization */
    section[data-testid="stSidebar"] {
        background-color: var(--primary-bg) !important;
        border-right: 1px solid var(--glass-border);
    }
    
    section[data-testid="stSidebar"] .stMarkdown p {
        color: var(--text-muted);
    }

    /* Typography Overrides */
    h1, h2, h3, h4 {
        font-family: 'Michroma', sans-serif !important;
        letter-spacing: 1px !important;
        text-transform: uppercase;
        color: var(--text-main) !important;
    }
    
    h1 {
        font-size: 2.5rem !important;
        margin-bottom: 1.5rem !important;
        border-left: 4px solid var(--accent-orange);
        padding-left: 1rem !important;
    }

    /* Cards / Containers with Glassmorphism */
    div[data-testid="stVerticalBlock"] > div:has(div.stMarkdown) {
        /* This targeting is tricky in Streamlit, applying to general containers */
    }
    
    .fazena-card {
        background: var(--card-bg);
        border: 1px solid var(--glass-border);
        border-radius: 12px;
        padding: 2rem;
        backdrop-filter: blur(12px);
        margin-bottom: 1rem;
        transition: all 0.3s ease;
    }
    
    .fazena-card:hover {
        border-color: var(--accent-orange);
        box-shadow: 0 0 20px var(--accent-glow);
    }

    /* Buttons - Fazena Orange */
    .stButton > button {
        background-color: transparent !important;
        color: var(--text-main) !important;
        font-family: 'Michroma', sans-serif;
        border: 1px solid var(--accent-orange) !important;
        border-radius: 4px !important;
        padding: 0.6rem 2rem !important;
        transition: all 0.4s cubic-bezier(0.165, 0.84, 0.44, 1) !important;
        text-transform: uppercase;
        letter-spacing: 2px;
        font-size: 0.8rem;
    }
    
    .stButton > button:hover {
        background-color: var(--accent-orange) !important;
        color: #000000 !important;
        box-shadow: 0 0 25px var(--accent-glow) !important;
        transform: scale(1.02);
    }

    /* Dataframes and Tables */
    .stDataFrame {
        border: 1px solid var(--glass-border);
        border-radius: 8px;
        background: var(--card-bg);
    }
    
    /* Sequence styling */
    .rna-sequence-container {
        font-family: 'JetBrains Mono', monospace;
        background: #050505;
        padding: 1.5rem;
        border-radius: 4px;
        border: 1px solid var(--glass-border);
        color: var(--accent-orange);
        word-break: break-all;
        line-height: 1.8;
        font-size: 1rem;
        margin-top: 1rem;
        box-shadow: inset 0 0 10px rgba(0,0,0,1);
    }
    
    /* Input Fields */
    .stTextInput > div > div > input, .stTextArea > div > textarea {
        background-color: #050505 !important;
        color: var(--text-main) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 4px !important;
        padding: 10px !important;
    }
    
    .stTextInput > div > div > input:focus, .stTextArea > div > textarea:focus {
        border-color: var(--accent-orange) !important;
        box-shadow: 0 0 10px var(--accent-glow) !important;
    }

    /* Tabs */
    .stTabs [data-baseweb="tab-list"] {
        background-color: transparent !important;
        gap: 20px;
    }

    .stTabs [data-baseweb="tab"] {
        font-family: 'Michroma', sans-serif !important;
        color: var(--text-muted) !important;
        background-color: transparent !important;
    }

    .stTabs [aria-selected="true"] {
        color: var(--accent-orange) !important;
        border-bottom-color: var(--accent-orange) !important;
    }
    </style>
    """, unsafe_allow_html=True)

def render_header(title, subtitle):
    st.markdown(f"<h1>{title}</h1>", unsafe_allow_html=True)
    st.markdown(f"<p style='color: var(--text-muted); font-size: 1.1rem; margin-bottom: 2rem;'>{subtitle}</p>", unsafe_allow_html=True)
