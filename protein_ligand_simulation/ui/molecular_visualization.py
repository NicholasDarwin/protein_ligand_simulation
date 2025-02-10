# ui/molecular_visualization.py
import streamlit as st
from visualization.py3dmol_viewer import visualize_structure

def display_simulation_results():
    st.title("Molecular Dynamics Simulation Results")
    pdb_file_path = st.text_input("Enter path to PDB file:", "output.pdb")
    
    if pdb_file_path:
        visualize_structure(pdb_file_path)
