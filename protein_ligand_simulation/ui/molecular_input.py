# molecular_input.py
import streamlit as st

def input_molecule(label="Upload Molecular File"):
    molecule = st.file_uploader(label, type=["sdf", "pdb", "mol2"])
    if molecule:
        st.write("Molecule loaded.")
        return molecule
    return None