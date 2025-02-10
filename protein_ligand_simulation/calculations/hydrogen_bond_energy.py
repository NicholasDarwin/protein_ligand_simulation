import numpy as np

def hydrogen_bond_energy(donor, acceptor):
    # Calculate distance between donor and acceptor
    r = np.linalg.norm(donor - acceptor)
    
    # Parameters for hydrogen bonding
    epsilon = 5.0  # kcal/mol (typical H-bond strength)
    r0 = 2.8      # Angstroms (typical H-bond length)
    
    # Angular dependent term could be added for more accuracy
    energy = epsilon * ((r0/r)**12 - 2*(r0/r)**6)
    return energy