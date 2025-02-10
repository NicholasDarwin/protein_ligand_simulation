import numpy as np
from scipy.spatial import ConvexHull

def solvent_accessible_surface_area(molecule):
    # Convert molecule coordinates to numpy array
    coords = np.array(molecule.GetConformer().GetPositions())
    
    # Add probe radius (typically 1.4Å for water)
    probe_radius = 1.4
    
    # Calculate convex hull
    hull = ConvexHull(coords)
    
    # Calculate surface area
    return hull.area + 4 * np.pi * probe_radius**2 * len(coords)