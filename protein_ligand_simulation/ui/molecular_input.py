import streamlit as st

def input_molecule():
    molecule = st.file_uploader("Upload Molecular File", type=["sdf", "pdb", "mol2"])
    if molecule:
        st.write("Molecule loaded.")
        return molecule
    return None
