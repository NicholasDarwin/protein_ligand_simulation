# calculations/molecular_dynamics.py
import openmm
from openmm import app
from openmm.unit import *

def run_simulation(pdb_file_path, force_field="amber14", temperature=300*kelvin, time_step=2*femtoseconds, steps=50000):
    # Load the structure and force field
    pdb = app.PDBFile(pdb_file_path)
    forcefield = app.ForceField(force_field)
    modeller = app.Modeller(pdb.topology, pdb.positions)

    # Create the system
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff)

    # Integrator and simulation setup
    integrator = openmm.LangevinIntegrator(temperature, 1/picosecond, time_step)
    simulation = app.Simulation(modeller.topology, system, integrator)

    # Set up the output
    simulation.reporters.append(app.DCDReporter('output.dcd', 500))  # Trajectory file
    simulation.reporters.append(app.StateDataReporter('output.txt', 1000, step=True, potentialEnergy=True, temperature=True))

    # Run the simulation
    simulation.step(steps)
    return simulation
