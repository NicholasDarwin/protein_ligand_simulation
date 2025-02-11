# Protein-Ligand MD Simulation

This app simulates molecular dynamics (MD) of protein-ligand complexes using GROMACS.

## Features
- Fetch protein and ligand data from the RCSB Protein Data Bank (PDB).
- Set up and run GROMACS MD simulations.
- Run simulations and analyze results.

## Installation

1. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

2. Install GROMACS and ensure it's accessible from your command line.

3. Run the app:
    ```bash
    python run.py
    ```

## Usage
- Enter the name of the protein and ligand.
- Click "Fetch PDB" to download the PDB file.
- Click "Run MD Simulation" to start the simulation.
