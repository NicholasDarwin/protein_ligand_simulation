# calculations/hydration_free_energy.py
import numpy as np
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import AllChem

def hydration_free_energy(ligand):
    """
    Calculate the hydration free energy of a ligand using an implicit solvent model.
    
    Parameters:
    -----------
    ligand : RDKit Mol object
        The ligand molecule for which to calculate hydration free energy
        
    Returns:
    --------
    float
        The calculated hydration free energy in kcal/mol
    """
    # Prepare the ligand by adding hydrogens and computing charges
    ligand = prepare_ligand(ligand)
    
    # Get atomic coordinates and properties
    conf = ligand.GetConformer()
    coords = conf.GetPositions()
    
    # Calculate polar contribution
    polar_contribution = calculate_polar_contribution(ligand)
    
    # Calculate non-polar contribution using SASA
    nonpolar_contribution = calculate_nonpolar_contribution(ligand)
    
    # Calculate cavity formation energy
    cavity_energy = calculate_cavity_formation(ligand)
    
    # Sum all contributions
    total_energy = polar_contribution + nonpolar_contribution + cavity_energy
    
    return total_energy

def prepare_ligand(ligand):
    """Prepare ligand by adding hydrogens and computing charges"""
    # Add hydrogens if they're not present
    if ligand.GetNumAtoms() != ligand.GetNumHeavyAtoms() + Chem.AddHs(ligand).GetNumAtoms() - ligand.GetNumHeavyAtoms():
        ligand = Chem.AddHs(ligand)
    
    # Generate 3D coordinates if not present
    if ligand.GetNumConformers() == 0:
        AllChem.EmbedMolecule(ligand, randomSeed=42)
    
    # Compute Gasteiger charges
    AllChem.ComputeGasteigerCharges(ligand)
    
    return ligand

def calculate_polar_contribution(ligand):
    """Calculate polar contribution using GB/SA model"""
    # Parameters
    epsilon_water = 78.5  # dielectric constant of water
    epsilon_vacuum = 1.0
    
    # Get Gasteiger charges
    charges = []
    for atom in ligand.GetAtoms():
        charge = atom.GetDoubleProp('_GasteigerCharge') if atom.HasProp('_GasteigerCharge') else 0.0
        charges.append(charge)
    
    # Calculate Born radii
    born_radii = [calculate_born_radius(atom) for atom in ligand.GetAtoms()]
    
    # Calculate polar solvation energy
    energy = 0.0
    for i, (qi, ri) in enumerate(zip(charges, born_radii)):
        energy += (qi**2)/(2*ri) * (1/epsilon_vacuum - 1/epsilon_water)
    
    return energy

def calculate_born_radius(atom):
    """Calculate Born radius for an atom based on its type"""
    # Default Born radii based on atom types (in Angstroms)
    born_radii = {
        'C': 1.7,
        'N': 1.55,
        'O': 1.52,
        'H': 1.2,
        'S': 1.8,
        'P': 1.8,
        'F': 1.47,
        'Cl': 1.75,
        'Br': 1.85,
        'I': 1.98
    }
    
    symbol = atom.GetSymbol()
    return born_radii.get(symbol, 1.7)  # Default to carbon radius if unknown

def calculate_nonpolar_contribution(ligand):
    """Calculate non-polar contribution using SASA"""
    # Surface tension parameter
    gamma = 0.0072  # kcal/(mol·Å²)
    
    # Calculate SASA using RDKit's built-in SASA calculation
    AllChem.ComputeMolVolume(ligand)  # This also computes SASA
    sasa = sum([atom.GetDoubleProp('_SASA') if atom.HasProp('_SASA') else 0.0 
                for atom in ligand.GetAtoms()])
    
    return gamma * sasa

def calculate_cavity_formation(ligand):
    """Calculate cavity formation energy"""
    # Volume-dependent term
    AllChem.ComputeMolVolume(ligand)
    volume = ligand.GetDoubleProp('_VMol') if ligand.HasProp('_VMol') else 0.0
    beta = 0.033  # kcal/(mol·Å³)
    
    return beta * volume