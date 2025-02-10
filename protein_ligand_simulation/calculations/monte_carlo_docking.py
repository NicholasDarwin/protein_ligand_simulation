import numpy as np
from scipy.spatial.transform import Rotation

def monte_carlo_docking(ligand, receptor, temperature, n_steps=1000):
    current_pose = ligand.GetConformer()
    current_energy = calculate_energy(ligand, receptor)
    best_pose = current_pose
    best_energy = current_energy
    
    kT = 0.0019872041 * temperature  # Boltzmann constant * temperature
    
    for step in range(n_steps):
        # Generate new pose through random rotation and translation
        new_pose = perturb_ligand(current_pose)
        new_energy = calculate_energy(new_pose, receptor)
        
        # Metropolis criterion
        if new_energy < current_energy or \
           np.random.random() < np.exp(-(new_energy - current_energy)/kT):
            current_pose = new_pose
            current_energy = new_energy
            
            if current_energy < best_energy:
                best_pose = current_pose
                best_energy = current_energy
    
    return best_pose, best_energy

def perturb_ligand(pose):
    # Random rotation
    rotation = Rotation.random()
    # Random translation
    translation = np.random.uniform(-2, 2, 3)
    
    new_pose = pose.copy()
    new_coords = rotation.apply(pose.GetPositions()) + translation
    new_pose.SetPositions(new_coords)
    return new_pose