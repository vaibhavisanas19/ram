import os
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
from Bio import PDB

# Streamlit webpage customization
st.set_page_config(page_title="Ramachandran Plot Generator", page_icon="🧬", layout="wide")
st.markdown(
    """
    <style>
        body {background-color: #f4f4f4;}
        .stApp {background-color: #e3f2fd;}
        h1 {color: #1565c0; text-align: center;}
        .stButton>button {background-color: #1976d2; color: white; border-radius: 10px;}
        .stFileUploader {text-align: center;}
    </style>
    """,
    unsafe_allow_html=True
)

# App Title
st.title("🧬 Ramachandran Plot Generator")
st.write("Upload a PDB file to analyze the phi and psi angles of protein structures.")

# File uploader
uploaded_file = st.file_uploader("Choose a PDB file", type=["pdb"], help="Upload a .pdb file")

if uploaded_file is not None:
    file_path = "uploaded.pdb"
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    # Initialize PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("uploaded", file_path)

    # PPBuilder to extract polypeptides
    ppb = PDB.PPBuilder()

    phi_angles = []
    psi_angles = []

    # Sidebar options
    st.sidebar.header("Options")
    show_residue_info = st.sidebar.checkbox("Show Residue Information", value=True)
    plot_color = st.sidebar.color_picker("Choose Plot Color", "#1976d2")

    # Loop through models and chains
    for model in structure:
        for chain in model:
            polypeptides = ppb.build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                if show_residue_info:
                    st.write(f"### Model {model.id}, Chain {chain.id}")
                    st.write(f"(Part {poly_index + 1} of {len(polypeptides)})")
                    st.write(f"Length: {len(poly)}")
                    st.write(f"From: {poly[0].resname} {poly[0].id[1]}")
                    st.write(f"To: {poly[-1].resname} {poly[-1].id[1]}")

                # Get phi-psi angles
                phi_psi_angles = poly.get_phi_psi_list()
                for res, angles in zip(poly, phi_psi_angles):
                    if angles[0] is not None and angles[1] is not None:
                        phi_angles.append(np.degrees(angles[0]))
                        psi_angles.append(np.degrees(angles[1]))
                    if show_residue_info:
                        st.write(f"Residue {res.resname} {res.id[1]}: Phi={angles[0]}, Psi={angles[1]}")

    # Plot Ramachandran plot using Streamlit
    st.write("## 📊 Ramachandran Plot")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(phi_angles, psi_angles, c=plot_color, marker='o', alpha=0.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("Phi (°)")
    ax.set_ylabel("Psi (°)")
    ax.set_title("Ramachandran Plot")
    ax.grid(True)
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axvline(0, color='black', linewidth=0.5)
    st.pyplot(fig)

    # Download button for processed data
    st.sidebar.download_button(
        label="Download Processed Data",
        data=f"Phi angles: {phi_angles}\nPsi angles: {psi_angles}",
        file_name="angles.txt",
        mime="text/plain"
    )
