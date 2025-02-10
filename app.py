# app.py
import streamlit as st
from ui.molecular_visualization import display_simulation_results
from calculations.molecular_dynamics import run_simulation

def main():
    st.title("Virtual Ligand Screening Pipeline")
    
    # Option to run molecular dynamics simulation
    if st.button("Run Molecular Dynamics Simulation"):
        pdb_file_path = st.text_input("Enter path to PDB file:", "input.pdb")
        if pdb_file_path:
            st.write("Running simulation...")
            run_simulation(pdb_file_path)
            st.write("Simulation completed!")
            display_simulation_results()

if __name__ == "__main__":
    main()
