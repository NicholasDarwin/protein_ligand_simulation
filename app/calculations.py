# app/calculations.py

from app.pdb_fetcher import fetch_pdb
from app.simulation import run_gromacs_simulation

def fetch_pdb_for_protein(protein_name):
    """Fetch the PDB file for the given protein."""
    return fetch_pdb(protein_name)

def start_simulation(protein_name, simulation_name):
    """Start the GROMACS simulation for the given protein and save the output."""
    output_message = run_gromacs_simulation(protein_name)
    # Here, we simulate that the result is the output message from the GROMACS simulation.
    return f"Simulation: {simulation_name}\nProtein: {protein_name}\n\n{output_message}"

