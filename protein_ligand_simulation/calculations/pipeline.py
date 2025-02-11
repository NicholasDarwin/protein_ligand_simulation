# calculations/pipeline.py
from rdkit import Chem
from rdkit.Chem import AllChem
from .mm_pbsa import mm_pbsa
from .docking_score import calculate_docking_score
from .hydration_free_energy import hydration_free_energy
from .monte_carlo_docking import monte_carlo_docking
from .molecular_dynamics import run_simulation

class CalculationPipeline:
    def __init__(self, ligand, receptor=None):
        self.ligand = self.prepare_molecule(ligand)
        self.receptor = self.prepare_molecule(receptor) if receptor else None
        self.results = {}
    
    def prepare_molecule(self, mol):
        """Prepare molecule for calculations"""
        if mol is None:
            return None
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Optimize geometry
        AllChem.MMFFOptimizeMolecule(mol)
        
        return mol
        
    def run_all_calculations(self):
        """Run all available calculations"""
        # Calculate docking score if receptor is available
        if self.receptor:
            self.results['docking_score'] = calculate_docking_score(self.ligand, self.receptor)
            self.results['mm_pbsa'] = mm_pbsa(self.ligand, self.receptor)
            
        # Calculate hydration free energy
        self.results['hydration_energy'] = hydration_free_energy(self.ligand)
        
        # Run molecular dynamics simulation
        try:
            simulation = run_simulation(self.ligand)
            self.results['md_simulation'] = simulation
        except Exception as e:
            self.results['md_simulation_error'] = str(e)
            
        return self.results