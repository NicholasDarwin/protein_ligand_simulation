import numpy as np

def radius_of_gyration(molecule):
    # Get atomic coordinates and masses
    coords = np.array(molecule.GetConformer().GetPositions())
    masses = np.array([atom.GetMass() for atom in molecule.GetAtoms()])
    
    # Calculate center of mass
    total_mass = np.sum(masses)
    com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    
    # Calculate radius of gyration
    rg = np.sqrt(np.sum(masses * np.sum((coords - com)**2, axis=1)) / total_mass)
    return rg