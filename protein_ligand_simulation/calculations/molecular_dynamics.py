# molecular_dynamics.py
import numpy as np
from .lennard_jones import lennard_jones_potential, lennard_jones_force
from .coulombs_law import coulombs_law
from .hydrogen_bond_energy import hydrogen_bond_energy
from .sasa import solvent_accessible_surface_area
from .rmsd import calculate_rmsd

class MolecularDynamics:
    def __init__(self, temperature=300, time_step=0.002, steps=1000):
        """
        Initialize molecular dynamics simulation
        
        Parameters:
        -----------
        temperature : float
            Temperature in Kelvin
        time_step : float
            Time step in picoseconds
        steps : int
            Number of simulation steps
        """
        self.temperature = temperature
        self.time_step = time_step
        self.steps = steps
        self.kb = 0.0019872041  # Boltzmann constant in kcal/mol/K
        
    def run_simulation(self, molecule):
        """
        Run molecular dynamics simulation
        
        Parameters:
        -----------
        molecule : dict
            Dictionary containing molecular information:
            - positions: np.array of atomic positions
            - velocities: np.array of atomic velocities
            - masses: np.array of atomic masses
            - charges: np.array of atomic charges
            - bonds: list of bonded atom pairs
        """
        # Initialize arrays for storing trajectory data
        trajectory = []
        energies = []
        
        # Get initial positions and velocities
        positions = molecule['positions']
        velocities = self._initialize_velocities(molecule['masses'])
        
        # Main simulation loop
        for step in range(self.steps):
            # Calculate forces
            forces = self._calculate_forces(molecule)
            
            # Update positions and velocities using Velocity Verlet algorithm
            positions, velocities = self._velocity_verlet_step(
                positions, velocities, forces, molecule['masses']
            )
            
            # Temperature control (simple velocity scaling)
            self._temperature_control(velocities, molecule['masses'])
            
            # Calculate energy
            energy = self._calculate_energy(molecule)
            
            # Store trajectory and energy
            trajectory.append(positions.copy())
            energies.append(energy)
            
            # Update molecule positions
            molecule['positions'] = positions
            
            # Print progress
            if step % 100 == 0:
                print(f"Step {step}/{self.steps}, Energy: {energy:.2f} kcal/mol")
        
        return trajectory, energies
    
    def _initialize_velocities(self, masses):
        """Initialize random velocities based on Maxwell-Boltzmann distribution"""
        velocities = np.random.normal(
            0, np.sqrt(self.kb * self.temperature / masses[:, np.newaxis]), 
            size=(len(masses), 3)
        )
        return velocities
    
    def _calculate_forces(self, molecule):
        """Calculate forces on each atom"""
        forces = np.zeros_like(molecule['positions'])
        
        # Non-bonded interactions
        for i in range(len(molecule['positions'])):
            for j in range(i + 1, len(molecule['positions'])):
                # Calculate distance
                r_ij = molecule['positions'][j] - molecule['positions'][i]
                distance = np.linalg.norm(r_ij)
                
                # Lennard-Jones force
                epsilon = 0.1  # kcal/mol
                sigma = 3.5    # Angstroms
                lj_force = lennard_jones_force(distance, epsilon, sigma)
                
                # Coulomb force
                coulomb_force = coulombs_law(
                    molecule['charges'][i], 
                    molecule['charges'][j], 
                    distance
                )
                
                # Total force
                total_force = (lj_force + coulomb_force) * r_ij / distance
                forces[i] += total_force
                forces[j] -= total_force
        
        # Bonded interactions (simple harmonic potential)
        k_bond = 100.0  # kcal/mol/A^2
        r0 = 1.5       # Equilibrium bond length
        for bond in molecule['bonds']:
            i, j = bond
            r_ij = molecule['positions'][j] - molecule['positions'][i]
            distance = np.linalg.norm(r_ij)
            force = -k_bond * (distance - r0) * r_ij / distance
            forces[i] += force
            forces[j] -= force
            
        return forces
    
    def _velocity_verlet_step(self, positions, velocities, forces, masses):
        """Implement Velocity Verlet integration"""
        # Update positions
        new_positions = positions + velocities * self.time_step + \
                       0.5 * forces / masses[:, np.newaxis] * self.time_step**2
                       
        # Calculate new forces
        new_forces = self._calculate_forces({'positions': new_positions})
        
        # Update velocities
        new_velocities = velocities + \
                        0.5 * (forces + new_forces) / masses[:, np.newaxis] * self.time_step
                        
        return new_positions, new_velocities
    
    def _temperature_control(self, velocities, masses):
        """Simple velocity scaling for temperature control"""
        kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
        current_temp = 2 * kinetic_energy / (3 * len(masses) * self.kb)
        
        scaling_factor = np.sqrt(self.temperature / current_temp)
        velocities *= scaling_factor
    
    def _calculate_energy(self, molecule):
        """Calculate total energy of the system"""
        # Kinetic energy
        kinetic = 0.5 * np.sum(molecule['masses'][:, np.newaxis] * 
                              molecule['velocities']**2)
        
        # Potential energy (non-bonded)
        potential = 0.0
        for i in range(len(molecule['positions'])):
            for j in range(i + 1, len(molecule['positions'])):
                r_ij = molecule['positions'][j] - molecule['positions'][i]
                distance = np.linalg.norm(r_ij)
                
                # Lennard-Jones
                potential += lennard_jones_potential(
                    distance, epsilon=0.1, sigma=3.5
                )
                
                # Coulomb
                potential += coulombs_law(
                    molecule['charges'][i], 
                    molecule['charges'][j], 
                    distance
                )
        
        return kinetic + potential