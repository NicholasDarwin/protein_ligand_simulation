import numpy as np

def calculate_rmsd(coords1, coords2):
    if coords1.shape != coords2.shape:
        raise ValueError("Input coordinates must have the same shape")
    
    # Calculate squared differences
    diff_sq = (coords1 - coords2)**2
    
    # Calculate mean squared difference
    mean_sq_diff = np.mean(np.sum(diff_sq, axis=1))
    
    # Return RMSD
    return np.sqrt(mean_sq_diff)