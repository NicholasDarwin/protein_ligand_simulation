import numpy as np
from scipy.linalg import eigh

def normal_mode_analysis(ligand, receptor):
    # Build Hessian matrix
    hessian = build_hessian_matrix(ligand, receptor)
    
    # Calculate mass-weighted Hessian
    masses = get_atomic_masses(ligand, receptor)
    mass_weighted_hessian = weight_hessian_by_mass(hessian, masses)
    
    # Calculate eigenvalues and eigenvectors
    eigenvals, eigenvecs = eigh(mass_weighted_hessian)
    
    # Calculate frequencies and entropy
    frequencies = np.sqrt(np.abs(eigenvals))
    entropy = calculate_vibrational_entropy(frequencies)
    
    return entropy