import numpy as np

def free_energy_perturbation(ligand1, ligand2, receptor):
    kB = 0.0019872041  # Boltzmann constant in kcal/mol/K
    T = 300  # Temperature in Kelvin
    beta = 1 / (kB * T)
    
    # Calculate energy differences between states
    dE = calculate_energy_difference(ligand1, ligand2, receptor)
    
    # Calculate free energy using Zwanzig equation
    dG = -kB * T * np.log(np.mean(np.exp(-beta * dE)))
    return dG