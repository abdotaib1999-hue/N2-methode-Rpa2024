import streamlit as st
import numpy as np
import pandas as pd
import math

# ==========================================
# PHASE 1: USER INPUTS (Sidebar)
# ==========================================
st.set_page_config(page_title="Seismic Vulnerability - N2 Method", layout="wide")
st.title("Détermination de la Performance Sismique (Méthode N2)")

st.sidebar.header("1. Paramètres de la Structure")
n_etages = st.sidebar.number_input("Nombre d'étages (n)", min_value=1, value=5, step=1)

# Dynamic table for mass and mode shapes
st.sidebar.subheader("Masses et Modes Propres")
default_data = {
    "Étage": [f"Niveau {i+1}" for i in range(n_etages)],
    "Masse (t)": [200.0] * n_etages,
    "Mode propre (φ)": np.linspace(0.2, 1.0, n_etages) # Linear default shape
}
df_structure = pd.DataFrame(default_data)
edited_df = st.sidebar.data_editor(df_structure, use_container_width=True, hide_index=True)

st.sidebar.header("2. Paramètres de Bilinéarisation")
Fy_star = st.sidebar.number_input("Force de fluage Fy* (kN)", value=1523.6, format="%.2f")
Sdy_star_mm = st.sidebar.number_input("Déplacement de fluage Sdy* (mm)", value=78.49, format="%.2f")
Sdy_star = Sdy_star_mm / 1000.0 # Convert to meters

st.sidebar.header("3. Paramètres Sismiques & Site")
A = st.sidebar.number_input("Coefficient d'accélération (A)", value=0.30)
S = st.sidebar.number_input("Paramètre de site (S)", value=1.30)
T1 = st.sidebar.number_input("T1 (s)", value=0.15)
T2 = st.sidebar.number_input("T2 (Tc) (s)", value=0.60)
T3 = st.sidebar.number_input("T3 (Td) (s)", value=2.00)
xi = st.sidebar.number_input("Amortissement ξ (%)", value=5.0)
eta = st.sidebar.number_input("Facteur de correction (η)", value=1.0)
k = st.sidebar.number_input("Exposant de sismicité (k)", value=2.7)


# ==========================================
# PHASE 2: PROGRAM CALCULATE (SDOF Properties)
# ==========================================
st.header("Propriétés du Système Équivalent (SDOF)")

# Extract arrays from editor
masses = edited_df["Masse (t)"].values
phis = edited_df["Mode propre (φ)"].values

# Calculations
m_star = np.sum(masses * phis)
M1 = np.sum(masses * (phis**2))
gamma = m_star / M1 if M1 != 0 else 0

Say = (Fy_star / m_star) / 9.81 if m_star != 0 else 0
K_eff = Say / Sdy_star if Sdy_star != 0 else 0
T_star = 2 * math.pi * math.sqrt((Sdy_star * m_star) / Fy_star) if Fy_star != 0 else 0

# Display Phase 2 Results
col1, col2, col3 = st.columns(3)
col1.metric("Masse Modale Effective (m*)", f"{m_star:.3f} t")
col2.metric("Facteur de Participation (Γ)", f"{gamma:.3f}")
col3.metric("Période Effective (T*)", f"{T_star:.3f} s")

col4, col5, col6 = st.columns(3)
col4.metric("Accélération de fluage (Say)", f"{Say:.3f} g")
col5.metric("Déplacement de fluage (Sdy*)", f"{Sdy_star:.4f} m")
col6.metric("Rigidité Effective (Keff)", f"{K_eff:.3f} g/m")

st.divider()

# ==========================================
# PHASE 3 & 4: SEISMIC DEMAND & TARGET DISPLACEMENT
# ==========================================
st.header("Évaluation de la Performance par État Limite")

# Define Limit States (Name, Tr*)
limit_states = [
    {"name": "Near Collapse (NC)", "Tr_star": 2475},
    {"name": "Significant Damage (SD)", "Tr_star": 475},
    {"name": "Damage Limitation (DL)", "Tr_star": 225},
    {"name": "RPA Standard", "Tr_star": 95},
    {"name": "FEMA 356", "Tr_star": 72}
]

# Function to calculate Sae(T*)
def calc_Sae(T, A, I, S, T1, T2, T3, eta):
    if T < T1:
        return A * I * S * (1 + (T / T1) * (2.5 * eta - 1))
    elif T1 <= T < T2:
        return A * I * S * (2.5 * eta)
    elif T2 <= T < T3:
        return A * I * S * (2.5 * eta) * (T2 / T) 
    else:
        return A * I * S * (2.5 * eta) * ((T2 * T3) / (T**2))

results = []

for ls in limit_states:
    Tr_star = ls["Tr_star"]
    
    # Calculate Importance Factor (I)
    I = (475 / Tr_star) ** (-1 / k)
    
    # Calculate Elastic Spectral Acceleration
    Sae_g = calc_Sae(T_star, A, I, S, T1, T2, T3, eta)
    
    # Phase 4: Target Displacement SDOF (dt*)
    dt_star = 0
    d_el_star = Sae_g * 9.81 * ((T_star / (2 * math.pi))**2) # using g = 9.81 m/s^2
    
    if T_star < T2:
        if Say >= Sae_g:
            dt_star = d_el_star
        else:
            qu = (Sae_g * m_star * 9.81) / Fy_star if Fy_star != 0 else 0
            # To avoid division by zero if qu is extremely small or exactly 0
            if qu > 0:
                dt_star = (1 + (qu - 1) * (T2 / T_star)) * (d_el_star / qu) if T_star != 0 else 0
            else:
                dt_star = d_el_star
    else:
        # Equal displacement rule
        dt_star = d_el_star
        
    # Final MDOF Target Displacement
    dt_mdof = gamma * dt_star
    
    # Check limit against elastic threshold
    elastic_limit_mdof = gamma * Sdy_star
    status = "Plastique (Dommages)" if dt_mdof > elastic_limit_mdof else "Élastique"
    
    results.append({
        "État Limite": ls["name"],
        "Période Retour (Tr*)": f"{Tr_star} ans",
        "Facteur (I)": f"{I:.3f}",
        "Sae(T*) [g]": f"{Sae_g:.3f}",
        "dt* SDOF (m)": f"{dt_star:.4f}",
        "dt MDOF (m)": f"{dt_mdof:.4f}",
        "État": status
    })

# Display Results in a clean DataFrame
df_results = pd.DataFrame(results)
st.dataframe(df_results, use_container_width=True)

# Add a visual warning if plastic limits are exceeded
if not df_results.empty:
    elastic_threshold = gamma * Sdy_star
    st.info(f"**Diagnostic :** Le déplacement élastique maximal (Toit) est de **{elastic_threshold:.4f} m**. Observez la colonne 'État' pour voir à quel moment la structure subit des dommages permanents.")
