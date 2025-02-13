# app/calculations.py

from app.pdb_fetcher import fetch_pdb
from app.simulation import run_gromacs_simulation

def fetch_pdb_for_protein(protein_name):
    """Fetch the PDB file for the given protein."""
    return fetch_pdb(protein_name)

def start_simulation(protein_name, ligand_name):
    # Call the simulation function with both parameters
    run_gromacs_simulation(protein_name, ligand_name)
    return "Simulation started successfully"