# app/calculations.py

from app.pdb_fetcher import fetch_pdb
from app.simulation import run_gromacs_simulation

def fetch_pdb_for_protein(protein_name):
    """Fetch the PDB file for the given protein."""
    return fetch_pdb(protein_name)

def start_simulation(protein_name, ligand_name):
    """Start the GROMACS simulation for the given protein and ligand."""
    output_message = run_gromacs_simulation(protein_name, ligand_name)
    return f"Protein: {protein_name}\nLigand: {ligand_name}\n\n{output_message}"