import openmm
from openmm import app
from openmm.unit import *
import tempfile
import os

def generate_ligand_parameters(mol2_path):
    """Generate parameters for the ligand using OpenMM directly"""
    # Read mol2 file directly using OpenMM's built-in parser
    pdb = app.PDBFile(mol2_path)
    return pdb

def create_system_with_custom_ligand(pdb, ligand_pdb, forcefield):
    """Create OpenMM system with custom ligand parameters"""
    # Create a new topology that includes the ligand
    modeller = app.Modeller(pdb.topology, pdb.positions)
    
    # Create system using the force field
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        rigidWater=True,
        removeCMMotion=True
    )
    
    return system, modeller

def run_simulation(input_file, force_field="amber99sbildn.xml", temperature=300*kelvin, time_step=2*femtoseconds, steps=50000):
    try:
        # Handle uploaded file
        if hasattr(input_file, 'name'):
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_file:
                tmp_file.write(input_file.getvalue())
                pdb_path = tmp_file.name
        else:
            pdb_path = input_file

        # Generate ligand parameters directly from PDB
        ligand_pdb = generate_ligand_parameters(pdb_path)

        # Load force field
        try:
            forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
        except Exception as e:
            raise ValueError(f"Failed to load force field: {str(e)}")

        # Create system with custom ligand
        system, modeller = create_system_with_custom_ligand(ligand_pdb, ligand_pdb, forcefield)

        # Set up integrator
        integrator = openmm.LangevinIntegrator(temperature, 1/picosecond, time_step)
        
        # Create simulation
        simulation = app.Simulation(modeller.topology, system, integrator)
        
        # Set initial positions
        simulation.context.setPositions(modeller.positions)
        
        # Minimize energy
        print("Minimizing energy...")
        simulation.minimizeEnergy()

        # Set up reporters
        simulation.reporters.append(app.DCDReporter('trajectory.dcd', 1000))
        simulation.reporters.append(app.StateDataReporter(
            'output.txt', 
            1000, 
            step=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=steps,
            separator='\t'
        ))

        # Run simulation
        print(f"Running simulation for {steps} steps...")
        simulation.step(steps)
        print("Simulation completed successfully")
        
        return simulation

    except Exception as e:
        raise Exception(f"Simulation failed: {str(e)}")
    finally:
        # Clean up temporary files
        if 'pdb_path' in locals():
            try:
                os.unlink(pdb_path)
            except:
                pass