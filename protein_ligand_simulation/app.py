import streamlit as st
from ui.molecular_visualization import display_simulation_results
from calculations.molecular_dynamics import run_simulation
from ui.molecular_input import input_molecule

def main():
    st.title("Virtual Ligand Screening Pipeline")
    
    # Add radio button to choose input method
    input_method = st.radio(
        "Choose input method:",
        ["Upload File", "Manual Input"]
    )
    
    # Initialize molecule variable
    molecule = None

    if input_method == "Upload File":
        # File upload logic
        molecule = input_molecule()
        if molecule:
            if st.button("Run Molecular Dynamics Simulation"):
                try:
                    st.write("Running simulation...")
                    run_simulation(molecule)
                    st.write("Simulation completed!")
                    display_simulation_results()
                except Exception as e:
                    st.error(f"Error during simulation: {str(e)}")
    
    else:
        # Manual input logic
        st.subheader("Manual Molecular Input")
        
        # Create input fields for molecular structure
        molecule_name = st.text_input("Molecule Name")
        
        # Add atomic structure inputs
        num_atoms = st.number_input("Number of Atoms", min_value=1, value=1)
        
        atoms_data = []
        for i in range(int(num_atoms)):
            st.write(f"Atom {i+1}")
            col1, col2, col3 = st.columns(3)
            with col1:
                atom_type = st.selectbox(f"Atom Type #{i+1}", 
                    ["C", "H", "O", "N", "P", "S"], key=f"type_{i}")
            with col2:
                x = st.number_input(f"X coordinate #{i+1}", key=f"x_{i}")
            with col3:
                y = st.number_input(f"Y coordinate #{i+1}", key=f"y_{i}")
            z = st.number_input(f"Z coordinate #{i+1}", key=f"z_{i}")
            
            atoms_data.append({
                "type": atom_type,
                "coordinates": (x, y, z)
            })

        # Create a mock molecule object after processing the atoms_data
        molecule = {
            "name": molecule_name,
            "atoms": atoms_data
        }
        
        if st.button("Process Molecule"):
            if molecule:
                st.write("Processing molecular structure...")
                st.write("Molecular structure processed!")
                st.write("Atoms data:", atoms_data)
                # Once the molecule is processed, run the simulation
                try:
                    st.write("Running simulation...")
                    run_simulation(molecule)
                    st.write("Simulation completed!")
                    display_simulation_results()
                except Exception as e:
                    st.error(f"Error during simulation: {str(e)}")
            else:
                st.error("Please input a valid molecular structure.")

if __name__ == "__main__":
    main()
